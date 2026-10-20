"""GO/KEGG over-representation of the declared available-source follow-up set."""
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
import csv,gzip,json,hashlib,argparse
import numpy as np
import pandas as pd
from scipy.stats import hypergeom,false_discovery_control
from audit_mr_selection import gene_locations

ROOT=Path(__file__).resolve().parents[2]

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def read_ontology(path):
 terms={};alt={};current=None;header=[]
 with path.open() as f:
  for line in f:
   line=line.rstrip()
   if line=='[Term]':
    current={'parents':set(),'alt':[],'obsolete':False};continue
   if line.startswith('['):current=None;continue
   if current is None:
    if line.startswith(('format-version:','data-version:','date:')):header.append(line)
    continue
   if line.startswith('id: '):
    current['id']=line[4:];terms[current['id']]=current
   elif line.startswith('name: '):current['name']=line[6:]
   elif line.startswith('namespace: '):current['namespace']=line[11:]
   elif line.startswith('is_a: '):current['parents'].add(line.split()[1])
   elif line.startswith('relationship: part_of '):current['parents'].add(line.split()[2])
   elif line.startswith('alt_id: '):current['alt'].append(line.split()[1])
   elif line=='is_obsolete: true':current['obsolete']=True
 for key,val in terms.items():
  if not val['obsolete']:
   for other in val['alt']:alt[other]=key
 terms={k:v for k,v in terms.items() if not v['obsolete']}
 return terms,alt,header

def enrich(sets,labels,background,query,min_size=10,max_size=500):
 annotated=set().union(*sets.values()) & background
 selected=query & annotated
 rows=[]
 for term,members in sorted(sets.items()):
  hitset=members & annotated;hits=hitset & selected
  if min_size<=len(hitset)<=max_size:
   p=hypergeom.sf(len(hits)-1,len(annotated),len(hitset),len(selected))
   rows.append(dict(term=term,description=labels.get(term,term),overlap=len(hits),query_annotated=len(selected),
                    term_background=len(hitset),background_annotated=len(annotated),p_value=float(p),gene_ids=';'.join(sorted(hits))))
 d=pd.DataFrame(rows)
 if len(d):d['fdr_bh']=false_discovery_control(d.p_value.to_numpy());d=d.sort_values(['fdr_bh','p_value','term'])
 return d,dict(background_annotated=len(annotated),query_annotated=len(selected),query_unannotated=sorted(query-annotated),tested_terms=len(d))

def main(scope='available'):
 raw=ROOT/'research/data/raw/functional';out=ROOT/('research/results/functional_'+scope);out.mkdir(parents=True,exist_ok=True)
 first=['research/results/mr_available/available_meta.tsv','research/results/available_followup/followup_proteins.tsv'] if scope=='available' else ['research/results/mr_fenland_substitute/all_proteins_comparison.tsv','research/results/mr_fenland_substitute/followup_proteins.tsv']
 paths=[ROOT/first[0],ROOT/first[1],
        ROOT/'research/data/metadata/Homo_sapiens.gene_info.gz',ROOT/'research/data/metadata/gencode.v19.annotation.gtf.gz',
        raw/'goa_human.gaf.gz',raw/'go-basic.obo',raw/'kegg_hsa_pathway_gene.tsv',raw/'kegg_hsa_pathway_names.tsv']
 original=pd.read_csv(paths[0],sep='\t');follow=pd.read_csv(paths[1],sep='\t')
 if scope=='substitute':
  identity_path=ROOT/'research/data/processed/substitute_mr/protein_registry.tsv';paths.append(identity_path)
  identities=pd.read_csv(identity_path,sep='\t').set_index('uniprot')
  original['gene']=original.uniprot.map(identities.gene_symbol)
 symbols=defaultdict(set);aliases=defaultdict(set)
 with gzip.open(paths[2],'rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['#tax_id']!='9606':continue
   symbols[r['Symbol']].add(r['GeneID'])
   for syn in r['Synonyms'].split('|'):
    if syn!='-':aliases[syn].add(r['GeneID'])
 locations=gene_locations(paths[3]);locations=locations[locations.gene_chr.isin(['chr'+str(i) for i in range(1,23)]+['chrX','chrY'])]
 locs={name:g for name,g in locations.groupby('gene')};mapping=[];uid_gene={}
 for _,r in original.iterrows():
  gene=str(r.gene);ids=symbols.get(gene,set()) or aliases.get(gene,set());loc=locs.get(gene)
  reason='eligible'
  if ';' in gene or ',' in gene:reason='conflicting_source_gene'
  elif loc is None or len(loc)!=1:reason='nonunique_gene_location'
  elif bool(loc.iloc[0].mhc_overlap):reason='MHC_excluded'
  elif len(ids)!=1:reason='nonunique_Entrez_mapping'
  mapped=next(iter(ids)) if len(ids)==1 else None
  mapping.append(dict(uniprot=r.uniprot,source_gene=gene,GeneID=mapped,background_status=reason))
  if reason=='eligible':uid_gene[r.uniprot]=mapped
 pd.DataFrame(mapping).to_csv(out/'background_mapping.tsv',sep='\t',index=False,na_rep='NA')
 background=set(uid_gene.values());query={uid_gene[u] for u in follow.uniprot if u in uid_gene}
 assert len(query)>0 and set(follow.uniprot).issubset(uid_gene)
 terms,alt,ontology_header=read_ontology(paths[5])
 @lru_cache(None)
 def ancestors(term):
  if term not in terms:return frozenset()
  parents=terms[term]['parents'];return frozenset({term}.union(*(ancestors(p) for p in parents)))
 go=defaultdict(set);gaf_header=[];unknown_terms=set();positive_annotations=0
 with gzip.open(paths[4],'rt') as f:
  for line in f:
   if line.startswith('!'):gaf_header.append(line.strip());continue
   r=line.rstrip('\n').split('\t')
   if r[0]!='UniProtKB' or r[1] not in uid_gene or r[12].split('|')[0]!='taxon:9606':continue
   if 'NOT' in r[3].split('|') or r[6]=='ND':continue
   term=alt.get(r[4],r[4])
   if term not in terms:unknown_terms.add(term);continue
   positive_annotations+=1
   for parent in ancestors(term):go[parent].add(uid_gene[r[1]])
 summaries={};outputs=[]
 for namespace in ['biological_process','molecular_function','cellular_component']:
  sets={k:v for k,v in go.items() if terms[k]['namespace']==namespace}
  table,summary=enrich(sets,{k:v['name'] for k,v in terms.items()},background,query)
  name='GO_'+namespace+'.tsv';table.to_csv(out/name,sep='\t',index=False);outputs.append(name)
  summary['significant_terms']=int((table.fdr_bh<.05).sum()) if len(table) else 0;summaries[namespace]=summary
 kegg=defaultdict(set)
 with paths[6].open() as f:
  for line in f:
   gene,pathway=line.rstrip().split('\t');kegg[pathway.removeprefix('path:')].add(gene.removeprefix('hsa:'))
 names={}
 with paths[7].open() as f:
  for line in f:
   key,name=line.rstrip().split('\t',1);names[key.removeprefix('path:')]=name
 table,summary=enrich(kegg,names,background,query);table.to_csv(out/'KEGG.tsv',sep='\t',index=False);outputs.append('KEGG.tsv')
 summary['significant_terms']=int((table.fdr_bh<.05).sum()) if len(table) else 0;summaries['KEGG']=summary
 record=dict(status='GO_KEGG_computed',formal_candidate_list_used=False,source_proteins=len(follow),query_genes=len(query),background_genes=len(background),
             selection=scope+' source follow-up; source gene conflicts/MHC/nonunique genes excluded from background',
             test='one-sided hypergeometric; BH separately within each GO namespace and KEGG; all eligible terms including zero overlaps',
             term_size=[10,500],GO_propagation=['is_a','part_of'],NOT_and_ND_excluded=True,
             summaries=summaries,gaf_header=gaf_header,ontology_header=ontology_header,unresolved_GO_terms=sorted(unknown_terms),
             inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},script_sha256=sha(Path(__file__)),
             outputs={p.name:sha(p) for p in out.glob('*.tsv')})
 (out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(summaries,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--scope',choices=['available','substitute'],default='available');main(p.parse_args().scope)
