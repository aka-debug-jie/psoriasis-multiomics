"""Export full cis eQTL records for frozen follow-up and author-comparison genes."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def registry():
 mapping={}
 annotation=ROOT/'research/data/metadata/Homo_sapiens.gene_info.gz'
 with gzip.open(annotation,'rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['#tax_id']!='9606':continue
   ids={x.split(':',1)[1] for x in r['dbXrefs'].split('|') if x.startswith('Ensembl:')}
   if ids:mapping.setdefault(r['Symbol'],set()).update(ids)
 genes={}
 followup=pd.read_csv(ROOT/'research/results/available_followup/followup_proteins.tsv',sep='\t')
 for _,r in followup.iterrows():genes[r.gene_id]=dict(gene_id=r.gene_id,symbol=r.mapped_gene,available_followup=True,author_reference=False)
 author=pd.read_csv(ROOT/'research/references/author_main_table1.tsv',sep='\t')
 for gene in author.target:
  symbol={'IL23':'IL23A','IL12':'IL12A','DDX58':'RIGI'}.get(gene,gene);ids=mapping.get(symbol,set())
  if len(ids)!=1:raise ValueError(f'Author-reference gene mapping not unique: {symbol}: {ids}')
  gid=next(iter(ids))
  if gid not in genes:genes[gid]=dict(gene_id=gid,symbol=symbol,available_followup=False,author_reference=True)
  else:genes[gid]['author_reference']=True
 return pd.DataFrame(genes.values()).sort_values('gene_id')

def main(tissues):
 settings=json.loads((ROOT/'research/config/smr.json').read_text())
 reg=registry()
 root=ROOT/'research/data/processed/coloc';root.mkdir(parents=True,exist_ok=True)
 reg.to_csv(root/'gene_registry.tsv',sep='\t',index=False)
 wanted=set();tissue_records=[]
 for tissue in tissues:
  data_dir=root/tissue;data_dir.mkdir(parents=True,exist_ok=True)
  smr=ROOT/'research/results/smr'/tissue/'all_probe_results.tsv'
  results=pd.read_csv(smr,sep='\t')
  indices=reg.merge(results[['ProbeID','gene_id','ProbeChr','Probe_bp','topSNP','topSNP_chr','topSNP_bp','p_SMR','p_HEIDI']],on='gene_id',how='left')
  indices['index_status']=indices.ProbeID.notna().map({True:'SMR_top_eQTL_available',False:'no_eligible_SMR_index'})
  indices.to_csv(data_dir/'index_registry.tsv',sep='\t',index=False,na_rep='NA')
  probes=data_dir/'probes.txt';probes.write_text('\n'.join(sorted(indices.ProbeID.dropna().unique()))+'\n')
  prefix=data_dir/'full_cis';target=Path(str(prefix)+'.txt')
  eqtl=ROOT/settings['datasets'][tissue]['prefix']
  cmd=[str(ROOT/settings['binary']),'--beqtl-summary',str(eqtl),'--extract-probe',str(probes),'--query','1','--out',str(prefix)]
  with (data_dir/'export.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
  count=0
  for chunk in pd.read_csv(target,sep='\t',chunksize=100000,usecols=['SNP']):
   count+=len(chunk);wanted.update(chunk.SNP)
  tissue_records.append(dict(tissue=tissue,rows=count,command=cmd,query_threshold=1,
                             source_prefix=str(eqtl.relative_to(ROOT)),SMR_results_sha256=digest(smr),
                             exported_sha256=digest(target),index_sha256=digest(data_dir/'index_registry.tsv')))
  print(tissue,'exported',count,'full cis pairs',flush=True)
 gwas=ROOT/'research/data/processed/smr_gwas/mvp_full.ma';subset=root/'mvp_regional.tsv'
 kept=0
 with gwas.open() as src,subset.open('w') as out:
  out.write(src.readline())
  for line in src:
   if line.split('\t',1)[0] in wanted:out.write(line);kept+=1
 record=dict(status='full_available_cis_queries_exported',tissues=tissue_records,requested_variant_ids=len(wanted),MVP_retained=kept,
             query_does_not_filter_by_significance=True,region_filter_stage='coloc script uses each SMR index, +/-500kb and +/-1Mb',
             source_coverage='Full published cis pairs only; locus windows can extend beyond source cis boundary',
             source_MVP_sha256=digest(gwas),regional_MVP_sha256=digest(subset),registry_sha256=digest(root/'gene_registry.tsv'),
             script_sha256=digest(Path(__file__)))
 (root/'run.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps(dict(requested=len(wanted),matched=kept)))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--tissues',nargs='+',default=['eQTLGen_Blood']);a=p.parse_args();main(a.tissues)
