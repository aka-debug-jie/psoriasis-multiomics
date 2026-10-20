"""Extract a fixed MVP EUR teaching panel from GWAS Catalog harmonised files."""
from pathlib import Path
import argparse,csv,json,math,zlib,concurrent.futures,os,hashlib,time,gzip
import requests,pandas as pd,pysam
os.environ.setdefault("CURL_CA_BUNDLE",requests.certs.where())
p=argparse.ArgumentParser();p.add_argument('--course',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);p.add_argument('--panel',type=Path,default=Path(__file__).resolve().parents[1]/'config/clinical_public_panel.tsv');a=p.parse_args()
a.cache.mkdir(parents=True,exist_ok=True)
family=pd.read_csv(a.course/'data/clinical/family_17.tsv',sep='\t')
x=pd.read_csv(a.course/'data/genetics/exposures.tsv',sep='\t');x=x[x.uniprot.isin(family.uniprot)]
positions=sorted({(str(int(r.chr)),int(r.pos))for r in x.itertuples()})
ids=set(x.rsid.astype(str))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def get(url,**kwargs):
 for attempt in range(2):
  try:
   response=requests.get(url,timeout=60,**kwargs);response.raise_for_status();return response
  except requests.RequestException:
   if attempt:raise
   time.sleep(2)
def work(row):
 accession=row['source_id'];folder=a.cache/accession;folder.mkdir(exist_ok=True)
 done=folder/'complete.json'
 if done.exists():return json.loads(done.read_text())
 number=int(accession[4:]);lo=((number-1)//1000)*1000+1;hi=lo+999
 base=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{lo:08}-GCST{hi:08}/{accession}/'
 api=get('https://www.ebi.ac.uk/gwas/rest/api/v2/studies/'+accession).json()
 (folder/'catalog.json').write_text(json.dumps(api,indent=2))
 meta=get(base+accession+'.tsv.gz-meta.yaml').text;(folder/'metadata.yaml').write_text(meta)
 # Keep exact official metadata and the directory listing for licence exception review.
 (folder/'directory.html').write_text(get(base).text)
 url=base+'harmonised/'+accession+'.h.tsv.gz'
 rows=[];seen=set();excluded=[]
 def accept(r):
  rs=r.get('rsid','');matches=ids & set(rs.replace(';',',').split(','))
  if not matches:return
  try:
   odds=float(r['odds_ratio']);low=float(r['ci_lower']);high=float(r['ci_upper'])
   if not(0<low<high and odds>0):raise ValueError('Invalid OR/CI')
   beta=math.log(odds);se=(math.log(high)-math.log(low))/(2*1.959963984540054)
   if not(math.isfinite(beta)and math.isfinite(se)and se>0):raise ValueError('Nonfinite effect')
   if int(r['num_cases'])!=int(row['num_cases'])or int(r['num_controls'])!=int(row['num_controls']):raise ValueError('Sample count mismatch')
  except (ValueError,KeyError)as e:
   excluded.append({'rsid':rs,'reason':str(e)});return
  for match in matches:rows.append({'rsid':match,'EA':r['effect_allele'],'OA':r['other_allele'],'beta':beta,'se':se,'p_reported':r['p_value']})
 index=folder/(accession+'.tbi')
 indexed=requests.get(url+'.tbi',timeout=60)
 if indexed.status_code==404:
  # Some Catalog records have not been harmonised. Stream the original public file;
  # retain exact source records and never substitute a similarly named phenotype.
  url=base+accession+'.tsv.gz'
  print(accession,'streaming original (no harmonised index)',flush=True)
  response=get(url,stream=True)
  try:
   with gzip.GzipFile(fileobj=response.raw)as stream:
    for chunk in pd.read_csv(stream,sep='\t',dtype=str,chunksize=250000):
     selected=chunk[chunk.rsid.isin(ids)]
     for r in selected.to_dict('records'):accept(r)
  finally:response.close()
 else:
  indexed.raise_for_status();index.write_bytes(indexed.content)
  response=get(url,headers={'Range':'bytes=0-16383'})
  if response.status_code!=206:raise ValueError('Server did not honor header range')
  header=zlib.decompressobj(31).decompress(response.content).decode().splitlines()[0].lstrip('#').split('\t')
  with pysam.TabixFile(url,index=str(index))as tab:
   for chrom,pos in positions:
    if chrom not in tab.contigs:continue
    for line in tab.fetch(chrom,pos-1,pos):
     if line in seen:continue
     seen.add(line);accept(dict(zip(header,line.split('\t'))))
 df=pd.DataFrame(rows,columns=['rsid','EA','OA','beta','se','p_reported']).drop_duplicates()
 dup=df.rsid.duplicated(keep=False)
 for rs in df.loc[dup,'rsid']:excluded.append({'rsid':rs,'reason':'conflicting source records'})
 df=df[~dup]
 if df.empty:raise ValueError(accession+' has no usable variants')
 out=folder/(accession+'.tsv');df.to_csv(out,sep='\t',index=False,na_rep='NA')
 (folder/'excluded.json').write_text(json.dumps(excluded,indent=2))
 result={**row,'url':url,'matched_variants':len(df),'sha256':digest(out),'metadata_sha256':digest(folder/'metadata.yaml'),'terms':'https://www.ebi.ac.uk/about/terms-of-use/','se_method':'log-CI width / (2*qnorm(.975)); original P retained separately'}
 done.write_text(json.dumps(result,indent=2));print(accession,len(df),'variants',flush=True);return result
def resilient(row):
 for attempt in range(2):
  try:return work(row)
  except (OSError,ValueError,requests.RequestException)as exc:
   print(row['source_id'],'retrieval attempt',attempt+1,'failed:',type(exc).__name__,flush=True)
   if attempt:return {'source_id':row['source_id'],'error':str(exc)}
   time.sleep(2)
with a.panel.open()as f:panel=list(csv.DictReader(f,delimiter='\t'))
with concurrent.futures.ThreadPoolExecutor(max_workers=3)as pool:results=list(pool.map(resilient,panel))
failed=[r for r in results if 'error'in r]
(a.cache/'last_failures.json').write_text(json.dumps(failed,indent=2))
if failed:raise RuntimeError('Some sources failed; completed inputs remain cached: '+','.join(r['source_id']for r in failed))
# Replace only the staged public teaching inputs, after all source extractions pass.
clinical=a.course/'data/clinical'
for file in clinical.glob('*.tsv'):
 if file.name!='family_17.tsv':file.unlink()
endpoints=[]
for row in results:
 accession=row['source_id'];(clinical/(accession+'.tsv')).write_bytes((a.cache/accession/(accession+'.tsv')).read_bytes())
 endpoints.append({'phenocode':accession,'phenotype':row['phenotype'],'num_cases':row['num_cases'],'num_controls':row['num_controls'],'path_https':row['url'],'scope':'public_teaching_MVP_EUR_PheCode_panel','comparison_note':row['comparison_note']})
with(clinical/'endpoints.tsv').open('w',newline='')as f:
 w=csv.DictWriter(f,list(endpoints[0]),delimiter='\t');w.writeheader();w.writerows(endpoints)
(clinical/'sources.json').write_text(json.dumps(results,indent=2))
print('All 22 public teaching outcomes prepared',flush=True)
