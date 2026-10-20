"""Select a declared two-source follow-up set; quarantine gene-label conflicts."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
import requests
import pandas as pd
from audit_mr_selection import gene_locations

ROOT=Path(__file__).resolve().parents[2]

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def main():
 out=ROOT/'research/results/available_followup';out.mkdir(parents=True,exist_ok=True)
 meta=ROOT/'research/data/metadata/uniprot_followup';meta.mkdir(parents=True,exist_ok=True)
 paths=[ROOT/'research/results/mr_available'/n for n in ['available_meta.tsv','finngen_available_meta.tsv','stuart_available_meta.tsv']]
 x,f,s=[pd.read_csv(p,sep='\t') for p in paths]
 for label,d in [('finngen',f),('stuart',s)]:
  x=x.merge(d[['uniprot','p','beta']].rename(columns={'p':'p_'+label,'beta':'beta_'+label}),on='uniprot',how='left')
 x=x[(x.p<1.89e-5)&((x.p_finngen<.05)|(x.p_stuart<.05))].copy()
 # Do not let a current UniProt name override conflicting assay identities.
 # A single published source gene must map uniquely to the primary assembly.
 mapping=[]
 for _,row in x.iterrows():
  gene=str(row.gene)
  conflict=';' in gene or ',' in gene
  mapping.append(dict(uniprot=row.uniprot,mapped_gene=None if conflict else gene,
                      mapping_status='cross_cohort_gene_conflict' if conflict else 'original_source_gene_label'))
 x=x.merge(pd.DataFrame(mapping),on='uniprot',validate='one_to_one')
 gtf=ROOT/'research/data/metadata/gencode.v19.annotation.gtf.gz'
 locations=gene_locations(gtf)
 locations=locations[locations.gene_chr.isin(['chr'+str(i) for i in range(1,23)]+['chrX','chrY'])]
 rows=[]
 for _,row in x.iterrows():
  g=locations[locations.gene==row.mapped_gene]
  data=row.to_dict();data['gene_id']=None
  if len(g)==1:
   item=g.iloc[0];data.update(item.drop(labels=['gene']).to_dict());data['gene_id']=item.gene_id.split('.')[0]
   flags=[bool(item[k]) for k in ['mhc_overlap','mhc_contained','mhc_tss']]
   data['eligibility']='MHC_excluded' if all(flags) else 'MHC_boundary_unresolved' if any(flags) else 'available_followup'
  else:data['eligibility']='gene_location_unresolved'
  data['replication_direction_concordant']=any(pd.notna(data['p_'+label]) and data['p_'+label]<.05 and data['beta']*data['beta_'+label]>0 for label in ['finngen','stuart'])
  data['scope']='two_available_sources_no_Steiger_not_formal_three_source_discovery'
  rows.append(data)
 d=pd.DataFrame(rows);d.to_csv(out/'all_nominally_replicated_proteins.tsv',sep='\t',index=False,na_rep='NA')
 chosen=d[d.eligibility=='available_followup'];chosen.to_csv(out/'followup_proteins.tsv',sep='\t',index=False,na_rep='NA')
 record=dict(status='partial_source_followup_set_frozen',date=datetime.now(timezone.utc).isoformat(),
             discovery_p=1.89e-5,external_p=.05,requires_one_external_dataset=True,requires_same_direction=False,
             gene_mapping='Single original source gene label then unique primary-chromosome GENCODE v19 interval; conflicting labels quarantined',
             MHC_rule='25-33Mb; overlap/containment/TSS must agree, otherwise unresolved',
             before_MHC=len(d),followup_proteins=len(chosen),eligibility_counts=d.eligibility.value_counts().to_dict(),
             formal_candidate_list=False,steiger_applied=False,
             inputs={str(p.relative_to(ROOT)):digest(p) for p in paths+[gtf]},
             script_sha256=digest(Path(__file__)),outputs={p.name:digest(p) for p in out.glob('*.tsv')})
 (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(chosen[['uniprot','mapped_gene','gene_id']].to_string(index=False));print(record['eligibility_counts'])

if __name__=='__main__':main()
