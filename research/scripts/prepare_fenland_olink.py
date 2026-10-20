"""Normalize the authorized Olink substitute, retaining raw units and estimated scales."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]

def digest(path,algorithm='sha256'):
 h=hashlib.new(algorithm)
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def main():
 meta_path=ROOT/'research/data/metadata/Fenland_Olink_6787142.json';meta=json.loads(meta_path.read_text())
 source_dir=ROOT/'research/data/raw/Fenland_Olink'
 for item in meta['files']:
  p=source_dir/item['key'];assert digest(p,'md5')==item['checksum'].split(':')[1],p
 target_path=source_dir/'Olink.protein.targets.Fenland.20220107.txt'
 targets=pd.read_csv(target_path,sep='\t')
 targets['id']=targets.pheno.astype(str)+'.'+targets.chr.astype(str)+'.'+targets.region_start.astype(str)+'.'+targets.region_end.astype(str)
 assert not targets.id.duplicated().any()
 multiplicity=targets.groupby('pheno').Uniprot_ID.nunique()
 targets['assay_multi_protein']=targets.pheno.map(multiplicity)>1
 out=ROOT/'research/data/processed/Fenland_Olink';out.mkdir(parents=True,exist_ok=True)
 targets.to_csv(out/'target_registry.tsv',sep='\t',index=False)
 by_id=targets.set_index('id');top={};moments={};rows=0;valid_rows=0;sample_sizes=set();regions_seen=set()
 requested=pd.read_csv(ROOT/'research/data/processed/coloc_all_genes/gene_registry.tsv',sep='\t')
 selected=set(targets.loc[targets.ensembl_gene_id.isin(requested.gene_id),'id'])
 selected_chunks=[]
 raw=source_dir/'Summary.statistics.cis.Olink.proteins.Fenland.20220107.txt.gz'
 for chunk in pd.read_csv(raw,sep='\t',chunksize=100000):
  rows+=len(chunk);regions_seen.update(chunk.id);sample_sizes.update(chunk.TotalSampleSize.dropna().unique().tolist())
  if not set(chunk.id).issubset(by_id.index):raise ValueError('Source region missing from matching target registry')
  good=np.isfinite(chunk.Effect)&np.isfinite(chunk.StdErr)&(chunk.StdErr>0)&chunk.Freq1.between(.001,.999)&(chunk['info']>=.4)&(chunk.TotalSampleSize>2)
  valid=chunk[good].copy();valid_rows+=len(valid)
  # Same no-intercept regression used by coloc::sdY.est, with reported per-SNP N.
  inverse_var=1/valid.StdErr**2
  n_genotype_var=2*valid.TotalSampleSize*valid.Freq1*(1-valid.Freq1)
  valid['scale_xy']=inverse_var*n_genotype_var
  valid['scale_xx']=inverse_var**2
  for identifier,g in valid.groupby('id',sort=False):
   agg=moments.setdefault(identifier,[0.,0.,0]);agg[0]+=g.scale_xy.sum();agg[1]+=g.scale_xx.sum();agg[2]+=len(g)
  eligible=valid[valid.rsid.astype(str).str.fullmatch(r'rs[0-9]+') & valid.Allele1.isin(list('ACGT')) & valid.Allele2.isin(list('ACGT')) & (valid.Allele1!=valid.Allele2)].copy()
  eligible['z_abs']=abs(eligible.Effect/eligible.StdErr)
  for identifier,g in eligible.groupby('id',sort=False):
   item=g.sort_values(['z_abs','rsid'],ascending=[False,True]).iloc[0].drop(labels=['scale_xy','scale_xx']).to_dict()
   if identifier not in top or (-item['z_abs'],item['rsid'])<(-top[identifier]['z_abs'],top[identifier]['rsid']):top[identifier]=item
  keep=chunk[chunk.id.isin(selected)]
  if len(keep):selected_chunks.append(keep)
  if rows%1000000==0:print('Fenland rows',rows,flush=True)
 scales=pd.DataFrame([dict(id=k,sdY_estimated=np.sqrt(v[0]/v[1]),scale_variants=v[2],scale_estimator='coloc_sdY_est_no_intercept_regression_reported_N') for k,v in moments.items()])
 scales=targets.merge(scales,on='id',how='left');scales.to_csv(out/'assay_scale_diagnostics.tsv',sep='\t',index=False)
 index=pd.DataFrame(top.values()).merge(scales,on='id',suffixes=('','_target'),validate='one_to_one')
 index['beta_estimated_SD']=index.Effect/index.sdY_estimated
 index['se_estimated_SD']=index.StdErr/index.sdY_estimated
 index['F_statistic']=(index.Effect/index.StdErr)**2
 index['source_published_threshold_pass']=index.log10p>-np.log10(4.5e-11)
 index['genomewide_5e8_threshold_pass']=index.log10p>8-np.log10(5)
 index['main_MR_assay_eligible']=~index.assay_multi_protein
 index['cohort']='Fenland_Olink_substitute_485'
 index['source_role']='authorized_substitute_not_Science2021_SomaScan'
 index.to_csv(out/'source_indices.tsv',sep='\t',index=False,na_rep='NA')
 regional=pd.concat(selected_chunks,ignore_index=True)
 regional.to_csv(out/'reference_and_followup_regions.tsv.gz',sep='\t',index=False)
 inputs=[raw,target_path,meta_path,ROOT/'research/data/processed/coloc_all_genes/gene_registry.tsv']
 record=dict(status='authorized_Fenland_Olink_input_prepared',source_rows=rows,valid_scale_rows=valid_rows,
             registered_regions=len(targets),regions_seen=len(regions_seen),unique_phenotypes=targets.pheno.nunique(),
             source_indices=len(index),multi_protein_regions=int(targets.assay_multi_protein.sum()),
             published_4_5e11_instruments=int((index.source_published_threshold_pass & index.main_MR_assay_eligible).sum()),
             genomewide_5e8_instruments=int((index.genomewide_5e8_threshold_pass & index.main_MR_assay_eligible).sum()),
             reported_sample_sizes=sorted(sample_sizes),selected_regional_pairs=len(regional),
             scale_summary=index.sdY_estimated.describe().to_dict(),
             scale_boundary='Original beta units not explicitly confirmed; raw beta/SE retained; optional estimated-SD rescaling is a reconstruction, not claimed original unit',
             inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs},script_sha256=digest(Path(__file__)),
             outputs={p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name!='run.json'},
             finished_utc=datetime.now(timezone.utc).isoformat())
 (out/'run.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps({k:v for k,v in record.items() if k not in ('inputs','outputs')},indent=2))

if __name__=='__main__':main()
