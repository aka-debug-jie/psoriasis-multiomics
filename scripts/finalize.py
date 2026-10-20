from pathlib import Path
import argparse,json,hashlib,zipfile
p=argparse.ArgumentParser();p.add_argument('--platform',required=True);p.add_argument('--language',choices=['zh-CN','en'],default='zh-CN');a=p.parse_args()
root=Path(__file__).resolve().parents[1];bundle=root/'.artifacts'/('pmo-'+a.language+'-'+a.platform)
res=bundle/'resources';m=json.loads((res/'bundle.json').read_text());v=json.loads((res/'native_validation.json').read_text())
assert m['status']=='native_execution_passed'and m['relocation_passed']
assert len(v['lessons'])==13 and all(x['status']=='passed'for x in v['lessons'])
# A build with unverified offline behavior is retained privately, never labeled ready for students.
assert v['network_isolation']=='passed','Offline isolation must pass before student release'
assert v['ordinary_user'] and v['offline_first_initialization']
for path in (res/'course/notebooks').glob('[0-9][0-9]_*.ipynb'):
 nb=json.loads(path.read_text(encoding='utf-8'))
 assert all(not c.get('outputs') and c.get('execution_count') is None for c in nb['cells'] if c['cell_type']=='code'),path
 assert path.read_bytes()==(res/'course/notebook_templates'/path.name).read_bytes(),path
files=sorted(p for p in bundle.rglob('*')if p.is_file() and p.name!='native_validation.json')
for file in files:
 relative=file.relative_to(bundle)
 assert not any(part in {'work','.git','.ipynb_checkpoints','__pycache__'}for part in relative.parts),relative
 assert file.name not in {'server.json','.env','.course-root'} and file.suffix not in {'.log','.key','.pem'},relative
assert {p.name for p in bundle.iterdir()}=={'resources','View course.html','Run course.bat'if a.platform.startswith('windows')else'Run course.command'}
if a.platform.startswith('windows'):assert (res/'bootstrap.ps1').read_bytes().startswith(b'\xef\xbb\xbf')
m['status']='ready_for_students'
out=root/'.artifacts/release';out.mkdir(exist_ok=True)
archive=out/(bundle.name+'.zip');temporary=archive.with_suffix('.zip.tmp')
with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED,compresslevel=6)as z:
 for f in files:
  if f==res/'bundle.json':z.writestr((Path(bundle.name)/f.relative_to(bundle)).as_posix(),json.dumps(m,ensure_ascii=False,indent=2))
  else:z.write(f,Path(bundle.name)/f.relative_to(bundle))
with zipfile.ZipFile(temporary)as z:assert z.testzip()is None
temporary.replace(archive)
(out/(a.language+'-'+a.platform+'-manifest.json')).write_text(json.dumps({'archive':archive.name,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'validation':v},ensure_ascii=False,indent=2),encoding='utf-8')
(res/'bundle.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
(res/'native_validation.json').unlink()
print('Student package ready:',archive)
