"""Run official SMR/HEIDI on all available probes and preserve unfiltered output."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT=Path(__file__).resolve().parents[2]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def symbol_mapping(path):
    mapping={}
    with gzip.open(path,'rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['#tax_id']!='9606':continue
            for ref in row['dbXrefs'].split('|'):
                if ref.startswith('Ensembl:'):
                    mapping.setdefault(ref.split(':',1)[1],set()).add(row['Symbol'])
    return {k:next(iter(v)) for k,v in mapping.items() if len(v)==1}


def main(config_path,tissue,output=None):
    config=json.loads(config_path.read_text());s=config['datasets'][tissue];settings=config['parameters']
    out=Path(output).resolve() if output else ROOT/'research/results/smr'/tissue
    out.mkdir(parents=True,exist_ok=True)
    prefix=out/'smr'
    exe=ROOT/config['binary'];gwas=ROOT/config['gwas'];ld=ROOT/config['LD_prefix'];eqtl=ROOT/s['prefix']
    inputs=[exe,gwas,config_path,Path(__file__),ROOT/config['gene_info'],ROOT/'research/references/author_main_table1.tsv']
    inputs += [Path(str(ld)+suffix) for suffix in ['.bed','.bim','.fam']]
    inputs += [Path(str(eqtl)+suffix) for suffix in ['.besd','.epi','.esi']]
    for path in inputs:
        if not path.is_file():raise FileNotFoundError(path)
    cmd=[str(exe),'--bfile',str(ld),'--gwas-summary',str(gwas),'--beqtl-summary',str(eqtl),
         '--peqtl-smr',str(settings['peqtl_smr']),'--heidi-skip-thresh',str(settings['heidi_skip_threshold']),
         '--peqtl-heidi',str(settings['peqtl_heidi']),'--ld-upper-limit',str(settings['heidi_ld_upper']),
         '--ld-lower-limit',str(settings['heidi_ld_lower']),'--cis-wind',str(settings['cis_window_kb']),
         '--diff-freq',str(settings['frequency_difference']),'--diff-freq-prop',str(settings['frequency_exclusion_fraction']),
         '--thread-num',str(settings['threads']),'--out',str(prefix)]
    record=dict(status='running',tissue=tissue,started_utc=datetime.now(timezone.utc).isoformat(),command=cmd,config=config,
                inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs})
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
    with (out/'execution.log').open('w') as log:
        process=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    record['exit_code']=process.returncode;record['finished_utc']=datetime.now(timezone.utc).isoformat()
    raw=Path(str(prefix)+'.smr')
    if process.returncode or not raw.is_file():
        record['status']='failed';(out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
        raise RuntimeError(f'SMR failed for {tissue}; inspect {out / "execution.log"}')
    record['status']='calculation_complete_postprocessing_pending'
    record['raw_smr_sha256']=digest(raw)
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
    d=pd.read_csv(raw,sep='\t',na_values=['NA','nan']).rename(columns={'probeID':'ProbeID'})
    required=['ProbeID','b_SMR','se_SMR','p_SMR','p_HEIDI']
    if not all(k in d for k in required):raise ValueError(f'Unexpected SMR columns: {d.columns.tolist()}')
    mapping=symbol_mapping(ROOT/config['gene_info'])
    d['gene_id']=d.ProbeID.str.replace(r'\.\d+$','',regex=True)
    d['symbol']=d.gene_id.map(mapping)
    d['gene_mapping_status']=np.where(d.symbol.notna(),'unique_current_NCBI_symbol','unmapped_or_ambiguous')
    for key in ['p_SMR','p_HEIDI']:
        d[key]=pd.to_numeric(d[key],errors='coerce').where(lambda x:x.between(0,1))
    d['SMR_FDR_all_tested_probes']=np.nan
    valid=d.p_SMR.notna()
    if valid.any():d.loc[valid,'SMR_FDR_all_tested_probes']=false_discovery_control(d.loc[valid,'p_SMR'].to_numpy())
    d['heidi_status']=np.where(d.p_HEIDI.isna(),'not_computed','computed')
    d['nominal_SMR_and_HEIDI_rule']=pd.Series(pd.NA,index=d.index,dtype='boolean')
    available=d.p_SMR.notna() & d.p_HEIDI.notna()
    d.loc[available,'nominal_SMR_and_HEIDI_rule']=(d.loc[available,'p_SMR']<.05)&(d.loc[available,'p_HEIDI']>.01)
    d['tissue']=tissue;d['candidate_scope']='all_eligible_genes_not_formal_protein_candidates'
    d.to_csv(out/'all_probe_results.tsv',sep='\t',index=False,na_rep='NA')
    d[d.nominal_SMR_and_HEIDI_rule.fillna(False)].to_csv(out/'nominal_supported_probes.tsv',sep='\t',index=False,na_rep='NA')
    author=pd.read_csv(ROOT/'research/references/author_main_table1.tsv',sep='\t')
    targets=author[['target']].copy();targets['symbol']=targets.target.replace({'IL23':'IL23A','IL12':'IL12A','DDX58':'RIGI'})
    comparison=targets.merge(d,on='symbol',how='left');comparison['evidence_origin']='author_reference_comparison_only'
    comparison.to_csv(out/'author_target_comparison.tsv',sep='\t',index=False,na_rep='NA')
    record.update(status='SMR_HEIDI_all_available_probes_completed',probes=len(d),heidi_computed=int(d.p_HEIDI.notna().sum()),
                  nominal_supported=int(d.nominal_SMR_and_HEIDI_rule.fillna(False).sum()),
                  FDR_05_probes=int((d.SMR_FDR_all_tested_probes<.05).sum()),formal_candidates_used=False,
                  note='Nominal SMR p<.05 and HEIDI p>.01 is the paper rule; genome-wide BH is reported separately; missing HEIDI stays missing',
                  outputs={p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name!='run.json'})
    (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['status','tissue','probes','heidi_computed','nominal_supported','FDR_05_probes']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,default=ROOT/'research/config/smr.json')
    parser.add_argument('--tissue',required=True);parser.add_argument('--output');a=parser.parse_args()
    main(a.config,a.tissue,a.output)
