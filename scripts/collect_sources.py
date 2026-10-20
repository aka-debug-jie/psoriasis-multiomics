"""Collect exact-version source offers and archives; never infer a release approval."""
from pathlib import Path
import argparse,json,hashlib,urllib.request,urllib.error,urllib.parse,tarfile,concurrent.futures,subprocess,shutil,re,threading
import yaml
from collect_component_metadata import collect,sha256_file,safe_part
p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--sources',type=Path,required=True);a=p.parse_args()
a.sources.mkdir(parents=True,exist_ok=True);manifest=collect(a.runtime,a.output)
ROOT=Path(__file__).resolve().parents[1]
lock=json.loads((ROOT/'course/renv.lock').read_text())['Packages']
_locks={};_locks_guard=threading.Lock()
def acquire(url,destination,expected=None):
 with _locks_guard:lock=_locks.setdefault(str(destination),threading.Lock())
 with lock:return _acquire(url,destination,expected)
def _acquire(url,destination,expected=None):
 destination.parent.mkdir(parents=True,exist_ok=True)
 if not destination.exists():
  request=urllib.request.Request(url,headers={'User-Agent':'psoriasis-multiomics-source-archive/4.0'})
  with urllib.request.urlopen(request,timeout=90)as response,destination.with_suffix(destination.suffix+'.part').open('wb')as f:shutil.copyfileobj(response,f)
  destination.with_suffix(destination.suffix+'.part').replace(destination)
 digest=sha256_file(destination)
 if expected and digest.lower()!=expected.lower():raise ValueError('Source checksum mismatch: '+destination.name)
 return {'url':url,'path':destination.relative_to(a.sources).as_posix(),'sha256':digest}
def first(urls,destination):
 errors=[]
 for url in urls:
  try:return acquire(url,destination)
  except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError)as e:errors.append(str(e))
 raise RuntimeError('; '.join(errors))
base_r={'base','compiler','datasets','graphics','grDevices','grid','methods','parallel','splines','stats','stats4','tcltk','tools','utils'}
R_SOURCE='https://cran.r-project.org/src/base/R-4/R-4.4.2.tar.gz'
def process(c):
 sources=[]
 try:
  name,version=c['name'],c['version']
  if c['ecosystem']=='r':
   if name in base_r:
    c['source_offer']={'url':R_SOURCE,'component':'R-4.4.2'}
   else:
    r=lock.get(name,{})
    if r.get('RemoteSha')and r.get('RemoteUsername')and r.get('RemoteRepo'):
     url=f"https://github.com/{r['RemoteUsername']}/{r['RemoteRepo']}/archive/{r['RemoteSha']}.tar.gz"
     sources=[acquire(url,a.sources/'r'/(name+'-'+version+'.tar.gz'))]
    else:
     filename=name+'_'+version+'.tar.gz'
     urls=['https://cran.r-project.org/src/contrib/'+filename,'https://cran.r-project.org/src/contrib/Archive/'+name+'/'+filename]
     if r.get('Source')=='Bioconductor' or c.get('source_metadata',{}).get('BiocViews'):
      urls=['https://bioconductor.org/packages/3.20/bioc/src/contrib/'+filename,'https://bioconductor.org/packages/3.20/bioc/src/contrib/Archive/'+name+'/'+filename]+urls
     sources=[first(urls,a.sources/'r'/filename)]
  elif c['ecosystem']=='python':
   # conda-pack is installed by conda; its wheel version may be a placeholder.
   same=[x for x in manifest['components']if x['ecosystem']=='conda'and x['name'].replace('-','_')==name.replace('-','_')]
   if same and version=='0.0.0':
    c['source_offer']={'component':'conda:'+same[0]['name'],'version':same[0]['version'],'reason':'Conda artifact and recipe identify this distribution; wheel metadata has a placeholder version.'}
   else:
    with urllib.request.urlopen('https://pypi.org/pypi/'+name+'/'+version+'/json',timeout=60)as response:info=json.load(response)
    distributions=[x for x in info['urls']if x['packagetype']=='sdist']
    if distributions:
     item=distributions[0];sources=[acquire(item['url'],a.sources/'python'/item['filename'],item['digests']['sha256'])]
    else:c['source_offer']={'url':info['info'].get('project_url'),'status':'No sdist; review installed wheel and upstream licence'}
  elif c['ecosystem']=='conda':
   directory=a.output/'conda'/(safe_part(name)+'-'+safe_part(version)+'-'+safe_part(c['build']))/'info/recipe'
   recipes=list(directory.glob('meta.yaml'))+list(directory.glob('recipe.yaml'))
   if recipes:
    try:recipe=yaml.safe_load(recipes[0].read_text())
    except yaml.YAMLError:recipe={}
    entries=recipe.get('source',[])if isinstance(recipe,dict)else[]
    if isinstance(entries,dict):entries=[entries]
    for i,item in enumerate(entries):
     if not isinstance(item,dict):continue
     urls=item.get('url',[])
     if isinstance(urls,str):urls=[urls]
     urls=[u for u in urls if isinstance(u,str)and u.startswith(('https://','http://'))and '{{'not in u and '${{'not in u]
     if urls:
      filename=Path(urllib.parse.urlsplit(urls[0]).path).name or 'source.tar.gz'
      entry=first(urls,a.sources/'conda'/(safe_part(name)+'-'+safe_part(version))/filename)
      if item.get('sha256')and isinstance(item['sha256'],str)and re.fullmatch('[0-9a-fA-F]{64}',item['sha256']):
       if entry['sha256']!=item['sha256']:raise ValueError('Conda source checksum mismatch')
      sources.append(entry)
     elif item.get('git_url')and item.get('git_rev'):
      c.setdefault('source_offers',[]).append({'git_url':item['git_url'],'git_rev':str(item['git_rev']),'recipe':recipes[0].relative_to(a.output).as_posix()})
   if not sources and 'source_offers'not in c:c['source_offer']={'artifact_url':c['artifact_url'],'status':'Review recipe, metapackage or binary redistribution terms'}
  elif c['ecosystem']=='jupyter':
   c['source_offer']={'component':'python:jupyterlab','note':'Frontend dependency notices accompany the JupyterLab source archive.'}
  c['sources']=sources;c['source_status']='collected_for_review'
 except Exception as exc:
  c['source_status']='collection_incomplete';c['collection_error']=str(exc)
 return c
# R core is shared by its standard-library components.
core=acquire(R_SOURCE,a.sources/'R-4.4.2.tar.gz')
for path in [a.runtime/'R/COPYING',a.runtime/'lib/R/COPYING']:
 if path.exists():
  (a.output/'r-core').mkdir(exist_ok=True);shutil.copy2(path,a.output/'r-core/COPYING')
with concurrent.futures.ThreadPoolExecutor(max_workers=6)as pool:
 manifest['components']=list(pool.map(process,manifest['components']))
manifest['r_core_source']=core
manifest['review_status']='pending_manual_review'
(a.output/'components.json').write_text(json.dumps(manifest,indent=2)+'\n')
(a.output/'SOURCE_ACCESS.md').write_text('''# Corresponding source and notices

This release includes component-specific licences and build recipes. The matching
source archives are offered alongside the student packages as source archives in
the same GitHub release. `components.json` lists the upstream source URLs, versions,
checksums, and any unresolved collection items. An archived file is not itself a
licence approval; the release review must resolve each required item before publication.

Rtsne includes software developed by Laurens van der Maaten. See its retained
licence for the complete attribution and conditions. Original course licensing
does not replace licences of R, Python, Jupyter, R packages, or external tools.
''')
print('Collected component sources; unresolved:',sum(c['source_status']=='collection_incomplete'for c in manifest['components']),flush=True)
