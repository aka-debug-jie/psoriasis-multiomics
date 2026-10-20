"""Execute the unpacked package with bundled software and compare to reference."""
from pathlib import Path
import argparse,subprocess,sys,os,shutil,json,importlib.util,csv,hashlib,platform as host_platform
p=argparse.ArgumentParser();p.add_argument('--platform',required=True);p.add_argument('--language',choices=['zh-CN','en'],default='zh-CN');p.add_argument('--inside',type=Path);p.add_argument('--smoke-only',action='store_true');a=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];win=sys.platform=='win32'
if a.inside is None:
 original=ROOT/'.artifacts'/('pmo-'+a.language+'-'+a.platform)
 testroot=ROOT/'.build'/'中文 空格路径'/original.name
 if testroot.exists():shutil.rmtree(testroot)
 shutil.copytree(original,testroot)
 subprocess.run([sys.executable,str(ROOT/'scripts/offline_native.py'),'--platform',a.platform,'--bundle',str(testroot)],check=True)
 # The tested files and outputs become the release candidate; no runtime cache is included.
 shutil.copytree(testroot/'resources/course/notebooks',original/'resources/course/notebooks',dirs_exist_ok=True)
 subprocess.run([sys.executable,str(ROOT/'scripts/prepare_student_notebooks.py'),str(original/'resources/course')],check=True)
 shutil.copy2(testroot/'resources/work/native_validation.json',original/'resources/native_validation.json')
 # Moving after first use must recreate only the disposable environment.
 marker=testroot/'resources/work/results/student_saved.txt';marker.write_text('keep my results')
 moved=testroot.with_name(testroot.name+' 已移动');testroot.rename(moved)
 subprocess.run([sys.executable,str(ROOT/'scripts/offline_native.py'),'--platform',a.platform,'--bundle',str(moved),'--smoke-only'],check=True)
 assert (moved/'resources/work/results/student_saved.txt').read_text()=='keep my results'
 shutil.rmtree(moved)
 manifest=json.loads((original/'resources/bundle.json').read_text());manifest.update(status='native_execution_passed',relocation_passed=True)
 (original/'resources/bundle.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
 sys.exit(0)

bundle=a.inside;res=bundle/'resources';course=res/'course'
with(course/'data/manifest.tsv').open(encoding='utf-8')as f:
 inputs=list(csv.DictReader(f,delimiter='\t'))
for item in inputs:
 path=course/item['path']
 assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256'], 'Input identity mismatch: '+item['path']
spec=importlib.util.spec_from_file_location('course_launch',res/'launch.py');launch=importlib.util.module_from_spec(spec);spec.loader.exec_module(launch)
env,rbin=launch.make_environment();os.environ.clear();os.environ.update(env)
import nbformat
from nbclient import NotebookClient
# Check the enforced boundary from inside the restricted account before computation.
assert os.environ.get('COURSE_OFFLINE_ENFORCED')=='1', 'Native offline enforcement missing'
import socket,ctypes
if win:
 assert not ctypes.windll.shell32.IsUserAnAdmin(), 'Test must run without administrator rights'
else:
 assert os.geteuid()!=0 and os.geteuid()!=501, 'Test must run as the isolated ordinary user'
for address in [('1.1.1.1',443),('2606:4700:4700::1111',443)]:
 try:
  connection=socket.create_connection(address,timeout=3)
 except OSError:pass
 else:
  connection.close();raise AssertionError('External network unexpectedly reachable')

from service_smoke import check_service
print('Checking the classroom service and R kernel',flush=True)
service=check_service(bundle)
print('Classroom service and R kernel passed',flush=True)
if a.smoke_only:
 print('Moved classroom server and R kernel passed');sys.exit(0)
statuses=[]
for path in sorted((course/'notebooks').glob('*.ipynb')):
 print('Executing '+path.name,flush=True)
 nb=nbformat.read(path,as_version=4)
 NotebookClient(nb,timeout=1200,kernel_name='ir-psoriasis',resources={'metadata':{'path':str(course/'notebooks')}},allow_errors=False).execute()
 nbformat.write(nb,path);statuses.append({'lesson':path.stem[:2],'status':'passed'})
 print('Passed '+path.name,flush=True)
subprocess.run([str(rbin),'--vanilla','--slave','-e','source("tests/check_results.R");check_course_results(Sys.getenv("TEACHING_NOTEBOOK_OUTPUT"))'],cwd=course,env=env,check=True)
subprocess.run([sys.executable,str(ROOT/'scripts/export_notebooks.py'),str(course),'--language',json.loads((res/'bundle.json').read_text()).get('language','zh-CN')],check=True,env=env)
# Compare deterministic tables to the shipped reference on the same data.
subprocess.run([str(rbin),'--vanilla','--slave','-f',str(ROOT/'scripts/compare_results.R')],cwd=course,env=env,check=True)
result={'platform':a.platform,'language':json.loads((res/'bundle.json').read_text()).get('language','zh-CN'),'tested_os':host_platform.platform(),'source_commit':json.loads((res/'bundle.json').read_text())['source_commit'],'R':'4.4.2','verified_input_files':len(inputs),'lessons':statuses,'result_checks':21,'numeric_comparison':'passed','network_isolation':'passed','ordinary_user':True,'offline_first_initialization':True,'service':service}
for name in ['numeric_comparison.tsv','single_cell_comparison.tsv']:
 with(res/'work/results'/name).open(encoding='utf-8')as f:result[name]=list(csv.DictReader(f,delimiter='\t'))
(res/'work/native_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('Native execution and reference comparisons passed')
