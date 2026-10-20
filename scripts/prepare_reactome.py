"""Prepare a fixed, redistributable Reactome teaching annotation; preserve GO inputs."""
from pathlib import Path
import argparse,csv,collections,hashlib,json,urllib.request
p=argparse.ArgumentParser();p.add_argument('--course',type=Path,default=Path('course'));p.add_argument('--cache',type=Path,required=True);a=p.parse_args()
a.cache.mkdir(parents=True,exist_ok=True)
url='https://download.reactome.org/92/NCBI2Reactome_All_Levels.txt'
source=a.cache/'NCBI2Reactome_All_Levels.txt'
if not source.exists():urllib.request.urlretrieve(url,source)
version='92'
def read(name):
 with(a.course/'data/functional'/name).open(encoding='utf-8')as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows,columns):
 with(a.course/'data/functional'/name).open('w',encoding='utf-8',newline='')as f:
  w=csv.DictWriter(f,columns,delimiter='\t');w.writeheader();w.writerows(rows)
background={r['GeneID']for r in read('opportunity_background.tsv')}
terms=collections.defaultdict(set);descriptions={};annotated=set()
with source.open(encoding='utf-8')as f:
 for row in csv.reader(f,delimiter='\t'):
  if len(row)<6 or row[5]!='Homo sapiens' or row[0]not in background:continue
  terms[row[1]].add(row[0]);descriptions[row[1]]=row[3];annotated.add(row[0])
old=read('term_gene.tsv');columns=list(old[0]);retained=[r for r in old if r['family'] in ('BP','MF','CC')]
for term,members in sorted(terms.items()):
 if 10<=len(members)<=500:
  for gene in sorted(members):retained.append({'family':'Reactome','term':term,'GeneID':gene,'description':descriptions[term]})
write('term_gene.tsv',retained,columns)
bg=read('annotated_background.tsv');write('annotated_background.tsv',[r for r in bg if r['family']in ('BP','MF','CC')]+[{'family':'Reactome','GeneID':g}for g in sorted(annotated)],list(bg[0]))
meta={'source':url,'version':version,'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'license':'CC0-1.0','license_url':'https://reactome.org/license','species':'Homo sapiens','selection':'All levels; opportunity background; 10<=term size<=500; full annotated background stored separately','annotated_background':len(annotated),'tested_terms':sum(10<=len(v)<=500 for v in terms.values())}
(a.course/'data/functional/reactome_source.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2))
