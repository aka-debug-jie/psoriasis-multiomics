"""Create a platform bundle. A candidate is not released until validation passes."""
from pathlib import Path
import argparse,subprocess,shutil,json,hashlib,csv,os,zipfile,platform
p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);p.add_argument('--language',choices=['zh-CN','en'],default='zh-CN');p.add_argument('--course',type=Path);p.add_argument('--platform',choices=['windows-x64','macos-intel','macos-arm64'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--smr',type=Path,required=True);p.add_argument('--plink',type=Path,required=True);p.add_argument('--runtime-archive',type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];bundle=a.output.resolve();res=bundle/'resources';res.mkdir(parents=True,exist_ok=True)
shutil.copytree(a.course or root/'course',res/'course',dirs_exist_ok=True)
templates=root/('translations/en/templates'if a.language=='en'else'templates')
for name,destination in [('view_course.html',bundle/'View course.html'),('help.html',res/'help.html')]:
 text=(templates/name).read_text().replace('课程资源','resources').replace('查看课程.html','View course.html').replace('运行课程.bat','Run course.bat').replace('运行课程.command','Run course.command')
 destination.write_text(text,encoding='utf-8')
shutil.copy2(root/'runtime/launch.py',res/'launch.py')
win=a.platform.startswith('windows');entry='Run course.bat' if win else 'Run course.command'
shutil.copy2(root/'runtime'/('start.bat'if win else 'start.command'),bundle/entry);(bundle/entry).chmod(0o755)
if win:shutil.copy2(root/'runtime/bootstrap.ps1',res/'bootstrap.ps1')
if win:(bundle/entry).write_bytes((bundle/entry).read_bytes().replace(b'\r\n',b'\n').replace(b'\n',b'\r\n'))
osname='windows'if win else'macos';arch='arm64'if a.platform.endswith('arm64')else'x86_64'
tools=res/'course/tools'/a.platform;tools.mkdir(parents=True,exist_ok=True)
rows=[]
for name,src in [('smr',a.smr),('plink',a.plink)]:
 dst=tools/(name+('.exe'if win else''));shutil.copy2(src,dst);dst.chmod(0o755)
 rows.append({'tool':name,'platform':osname,'arch':arch,'version':'1.4.3'if name=='smr'else'1.9.0','url':'','sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'origin':'native_build_with_documented_patches'if name=='smr'else'official_binary','path':str(dst.relative_to(res/'course')).replace('\\','/')})
with(res/'course/tools/manifest.tsv').open('w',newline='',encoding='utf-8')as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
# Include locally linked DLLs required by the Windows build.
if win:
 for dll in a.smr.parent.glob('*.dll'):shutil.copy2(dll,tools/dll.name)
seed=res/'runtime'/('seed.zip'if win else'seed.tar.gz');seed.parent.mkdir(exist_ok=True)
(res/'runtime/platform.txt').write_text(a.platform+'\n',encoding='utf-8')
if a.runtime_archive:shutil.copy2(a.runtime_archive,seed)
else:subprocess.run(['conda-pack','-p',str(a.runtime.resolve()),'-o',str(seed),'--format','zip'if win else'tar.gz','--force'],check=True)
manifest={'language':a.language,'version':'4.0.0','platform':a.platform,'source_commit':os.environ.get('GITHUB_SHA','local'),'runtime_sha256':hashlib.sha256(seed.read_bytes()).hexdigest(),'R':'4.4.2','status':'candidate_requires_native_validation','built_system':platform.platform()}
(res/'runtime/language.txt').write_text(a.language+'\n',encoding='utf-8')
(res/'runtime/sha256.txt').write_text(manifest['runtime_sha256']+'\n')
expanded=sum(f.stat().st_size for f in a.runtime.rglob('*')if f.is_file()and not f.is_symlink())
(res/'runtime/min_free_kb.txt').write_text(str(int((expanded*1.2+512*1024**2)/1024))+'\n')
(res/'bundle.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('Candidate bundle:',bundle)
