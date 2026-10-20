"""Export full eQTL regions for the version-matched substituted target registry.

Same extraction algorithm as prepare_coloc_all_genes.py; preserves the previous branch.
"""
from pathlib import Path
import json
import subprocess
import numpy as np
import pandas as pd
from prepare_coloc_inputs import ROOT,registry,digest


def main():
 settings=json.loads((ROOT/'research/config/smr.json').read_text())
 reg=pd.read_csv(ROOT/'research/data/processed/coloc_substitute_registry/gene_registry.tsv',sep='\t');root=ROOT/'research/data/processed/coloc_mrna_substitute';root.mkdir(parents=True,exist_ok=True)
 reg.to_csv(root/'gene_registry.tsv',sep='\t',index=False)
 wanted=set();records=[]
 for tissue,source in settings['datasets'].items():
  out=root/tissue;out.mkdir(parents=True,exist_ok=True)
  prefix=ROOT/source['prefix']
  epi=pd.read_csv(str(prefix)+'.epi',sep='\t',header=None,names=['ProbeChr','ProbeID','genetic_distance','Probe_bp','source_gene','orientation'])
  epi['gene_id']=epi.ProbeID.str.replace(r'\.\d+$','',regex=True)
  targets=reg.merge(epi,on='gene_id',how='left')
  probes=out/'probes.txt';probes.write_text('\n'.join(sorted(targets.ProbeID.dropna().unique()))+'\n')
  exported=out/'full_cis.txt'
  cmd=[str(ROOT/settings['binary']),'--beqtl-summary',str(prefix),'--extract-probe',str(probes),'--query','1','--out',str(out/'full_cis')]
  with (out/'export.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
  q=pd.read_csv(exported,sep='\t')
  wanted.update(q.SNP)
  indices=[]
  for _,row in targets.iterrows():
   record=row.to_dict();record['topSNP']=None;record['topSNP_chr']=None;record['topSNP_bp']=None
   record['index_selection']='maximum_source_abs_beta_over_SE_before_outcome_match'
   if pd.isna(row.ProbeID):
    record['index_status']='probe_absent_from_source';indices.append(record);continue
   g=q[q.Probe==row.ProbeID].drop_duplicates().copy()
   bad=set(g.loc[g.SNP.duplicated(keep=False),'SNP'])
   g=g[~g.SNP.isin(bad)]
   g=g[np.isfinite(g.b)&np.isfinite(g.SE)&(g.SE>0)]
   record['source_conflicting_index_variants_removed']=len(bad)
   if g.empty:
    record['index_status']='no_valid_source_cis_pair';indices.append(record);continue
   g['absolute_z']=abs(g.b/g.SE)
   top=g.sort_values(['absolute_z','SNP'],ascending=[False,True]).iloc[0]
   record.update(topSNP=top.SNP,topSNP_chr=top.Chr,topSNP_bp=top.BP,index_absolute_z=top.absolute_z,
                 index_source_p=top.p,index_status='source_index_frozen_without_association_P_threshold')
   indices.append(record)
  idx=out/'index_registry.tsv';pd.DataFrame(indices).to_csv(idx,sep='\t',index=False,na_rep='NA')
  records.append(dict(tissue=tissue,rows=len(q),command=cmd,query_threshold=1,
                      eqtl_BESD_sha256=digest(Path(str(prefix)+'.besd')),
                      eqtl_ESI_sha256=digest(Path(str(prefix)+'.esi')),eqtl_EPI_sha256=digest(Path(str(prefix)+'.epi')),
                      exported_sha256=digest(exported),index_sha256=digest(idx)))
  print(tissue,'pairs',len(q),'indices',sum(r['topSNP'] is not None for r in indices),flush=True)
 gwas=ROOT/'research/data/processed/smr_gwas/mvp_full.ma';subset=root/'mvp_regional.tsv';kept=0
 with gwas.open() as src,subset.open('w') as out:
  out.write(src.readline())
  for line in src:
   if line.split('\t',1)[0] in wanted:out.write(line);kept+=1
 report=dict(status='all_registered_genes_full_available_cis_exported',tissues=records,requested_variant_ids=len(wanted),MVP_retained=kept,
             index_selection='max abs(b/SE) on valid, unambiguous source variants; SNP lexical tie break; before MVP matching; never substitute index',
             query_does_not_filter_by_significance=True,source_coverage='published full cis pairs; no extension beyond source support',
             source_MVP_sha256=digest(gwas),regional_MVP_sha256=digest(subset),registry_sha256=digest(root/'gene_registry.tsv'),
             script_sha256=digest(Path(__file__)),registry_helper_sha256=digest(ROOT/'research/scripts/prepare_coloc_inputs.py'))
 (root/'run.json').write_text(json.dumps(report,indent=2)+'\n')
 print('MVP shared',kept)

if __name__=='__main__':main()
