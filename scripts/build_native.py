"""Build on the target OS; the output is still a candidate until test_bundle passes."""
from pathlib import Path
import argparse,os,subprocess,sys,shutil,urllib.request,zipfile,json,re,tarfile,hashlib
p=argparse.ArgumentParser();p.add_argument('--platform',required=True);p.add_argument('--smr-diagnostic',action='store_true');p.add_argument('--trace-heidi',action='store_true');a=p.parse_args()
assert not a.trace_heidi or a.smr_diagnostic
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
BUILD=ROOT/'.build';BUILD.mkdir(exist_ok=True);RUNTIME=BUILD/('runtime' if sys.platform=='win32' else 'runtime-'+('p'*120))
conda=os.environ.get('CONDA_EXE')
if not conda:
 candidate=Path(os.environ.get('CONDA',sys.prefix))/('Scripts/conda.exe'if sys.platform=='win32'else'bin/conda')
 conda=str(candidate)if candidate.exists()else shutil.which('conda')
win=sys.platform=='win32'
expected={'windows-x64':('win32','AMD64'),'macos-intel':('darwin','x86_64'),'macos-arm64':('darwin','arm64')}
import platform
assert sys.platform==expected[a.platform][0] and platform.machine().lower()==expected[a.platform][1].lower(),(sys.platform,platform.machine())
def run(argv,**kw):print('Run:',argv[0],flush=True);return subprocess.run([str(x)for x in argv],check=True,**kw)
from runtime_cache import restore,save
runtime_archive=restore(ROOT,a.platform,RUNTIME)
if runtime_archive is None and not(a.smr_diagnostic and RUNTIME.is_dir()):
    packages=['python=3.12','pip','conda-pack','zlib','libuv','zeromq','libpng','pkg-config']
    if not win:
     packages+=['c-compiler','cxx-compiler','make','cmake','llvm-openmp']
     if not a.smr_diagnostic:packages+=['r-base=4.4.2','fortran-compiler','libxml2','libcurl','openssl']
    run([conda,'create','-y','-p',RUNTIME,'--override-channels','-c','conda-forge',*packages])
    python=RUNTIME/('python.exe'if win else'bin/python')
    if not a.smr_diagnostic:run([python,'-m','pip','install','-r','requirements-build.txt'])
    if win and not a.smr_diagnostic:
     rscript=shutil.which('Rscript')or shutil.which('Rscript.exe')
     if not rscript:raise RuntimeError('Build runner requires R 4.4.2 from setup-r')
     rhome=subprocess.check_output([rscript,'--vanilla','-e','cat(R.home())'],text=True).strip()
     shutil.copytree(rhome,RUNTIME/'R',dirs_exist_ok=True)
     rscript=RUNTIME/'R/bin/Rscript.exe'
    elif not win:rscript=RUNTIME/'bin/Rscript'
    if not a.smr_diagnostic:run([conda,'run','--no-capture-output','-p',RUNTIME,rscript,'--vanilla','scripts/install_r_packages.R',RUNTIME])
    runtime_archive=None
source=BUILD/'smr-source';commit='eedb1a62ba4c66ab0b09acd0205876955a803f18'
if not source.exists():
 source.mkdir();run(['git','-C',source,'init']);run(['git','-C',source,'fetch','--depth','1','https://github.com/JianYang-Lab/SMR.git',commit]);run(['git','-C',source,'checkout','--detach','FETCH_HEAD'])
assert subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()==commit
if win and not a.trace_heidi:run([sys.executable,'scripts/patch_smr_platform.py',source])
if a.trace_heidi:run([sys.executable,'scripts/trace_smr_heidi.py',source])
smrbuild=BUILD/'smr-build';install=BUILD/'smr'
cmake=['cmake','-S',source,'-B',smrbuild,'-G','Ninja'if win else'Unix Makefiles','-DCMAKE_BUILD_TYPE=Release',f'-DCMAKE_INSTALL_PREFIX={install}','-DSMR_BLAS_BACKEND=none','-DBUILD_WITH_MIMALLOC=OFF','-DSMR_HAS_MARCH_X86_64_V3=FALSE','-DSMR_HAS_MARCH_ARMV85A=FALSE']
if win:
 ucrt=Path(os.environ['MSYS2_ROOT'])/'ucrt64/bin';env=os.environ.copy();env['PATH']=str(ucrt)+os.pathsep+env['PATH']
 assert (ucrt/'gcc.exe').is_file(),'Missing MSYS2 UCRT compiler'
 cmake[0]=str(ucrt/'cmake.exe')
 cmake +=[f'-DCMAKE_C_COMPILER={ucrt}/gcc.exe',f'-DCMAKE_CXX_COMPILER={ucrt}/g++.exe']
 run(cmake,env=env);run(['cmake','--build',smrbuild,'--parallel','2'],env=env);run(['cmake','--install',smrbuild],env=env)
 smr=install/'bin/smr.exe'
 # Resolve non-system import DLLs recursively, keeping them beside the executable.
 candidates={p.name.lower():p for p in ucrt.glob('*.dll')}
 candidates.update({p.name.lower():p for p in smrbuild.rglob('*.dll')})
 queue=[smr];seen=set()
 while queue:
  binary=queue.pop();text=subprocess.check_output([str(ucrt/'objdump.exe'),'-p',str(binary)],text=True)
  for name in re.findall(r'DLL Name:\s*(\S+)',text):
   key=name.lower()
   if key in seen:continue
   seen.add(key)
   if key in candidates:
    dst=smr.parent/candidates[key].name;shutil.copy2(candidates[key],dst);queue.append(dst)
else:
 cmake+=['-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0','-DCMAKE_AR=/usr/bin/ar','-DCMAKE_RANLIB=/usr/bin/ranlib']
 run([conda,'run','--no-capture-output','-p',RUNTIME,*cmake]);run([conda,'run','--no-capture-output','-p',RUNTIME,'cmake','--build',smrbuild,'--parallel','2']);run(['cmake','--install',smrbuild])
 smr=install/'bin/smr'
 # Bundle non-system dylibs so SMR does not depend on Homebrew or the build prefix.
 queue=[smr];seen=set()
 while queue:
  binary=queue.pop();deps=subprocess.check_output(['otool','-L',str(binary)],text=True).splitlines()[1:]
  for line in deps:
   dep=line.strip().split(' (')[0]
   if dep.startswith(('/usr/lib/','/System/','@loader_path/')):continue
   path=Path(dep)
   if dep.startswith('@rpath/'):path=RUNTIME/'lib'/Path(dep).name
   if not path.exists():raise RuntimeError('Unresolved SMR dylib: '+dep)
   target=smr.parent/path.name
   if path.name not in seen:
    seen.add(path.name);shutil.copy2(path,target);target.chmod(0o755);queue.append(target)
   run(['install_name_tool','-change',dep,'@loader_path/'+path.name,binary])
 for item in [smr,*smr.parent.glob('*.dylib')]:run(['codesign','--force','--sign','-',item])
if a.smr_diagnostic:
 run([sys.executable,'scripts/diagnose_smr.py','--smr',smr,'--platform',a.platform,'--label','trace'if a.trace_heidi else'baseline'])
 sys.exit(0)
url='https://s3.amazonaws.com/plink1-assets/'+('plink_win64_20260927.zip'if win else'plink_mac_20260927.zip')
archive=BUILD/'plink.zip';urllib.request.urlretrieve(url,archive)
with zipfile.ZipFile(archive)as z:z.extractall(BUILD/'plink')
plink=next((BUILD/'plink').rglob('plink.exe'if win else'plink'));plink.chmod(0o755)
expected_hash='efc55270aba7216e4fcebdf781ee04e9fa556c2fe5be1f71c836ba73930a6b6b'if win else'13a7fc9b908d0feb968d0835108b7fa5416fe27f44cad9bb7a0887676a644c40'
assert hashlib.sha256(plink.read_bytes()).hexdigest()==expected_hash,'PLINK source identity changed'
help_result=subprocess.run([str(smr)],capture_output=True,text=True)
assert 'Version 1.4.3'in help_result.stdout+help_result.stderr,'Unexpected SMR version'
# Capture installed components, exact source archives, dependency trees and build patches.
notices=BUILD/('notices-'+a.platform);sources=BUILD/('sources-'+a.platform)
run([conda,'run','--no-capture-output','-p',RUNTIME,'python','scripts/collect_sources.py','--runtime',RUNTIME,'--output',notices,'--sources',sources])
sources.mkdir(exist_ok=True)
plink_commit='bceae81647e473639d6a4a99be8d23821ae3566f'
plink_source_url='https://github.com/chrchang/plink-ng/archive/'+plink_commit+'.tar.gz'
plink_source=sources/'PLINK-1.9.0-source.tar.gz';urllib.request.urlretrieve(plink_source_url,plink_source)
with tarfile.open(sources/'SMR-1.4.3-complete-source.tar.gz','w:gz')as archive:
 for base,label in [(source,'SMR-1.4.3'),(smrbuild/'_deps','dependencies')]:
  if not base.exists():continue
  for f in base.rglob('*'):
   relative=f.relative_to(base)
   if '.git'in relative.parts or any(part.endswith('-build')for part in relative.parts):continue
   if f.is_file()and not f.is_symlink():archive.add(f,arcname=Path(label)/relative)
shutil.copy2(ROOT/'scripts/patch_smr_platform.py',sources/'SMR_windows_port.py')
(notices/'external_tools.json').write_text(json.dumps({'SMR':{'version':'1.4.3','commit':commit,'build_arguments':[str(x)for x in cmake],'complete_source_archive':'SMR-1.4.3-complete-source.tar.gz'},'PLINK':{'version':'1.9.0','commit':plink_commit,'source_url':plink_source_url,'source_sha256':hashlib.sha256(plink_source.read_bytes()).hexdigest(),'binary_sha256':expected_hash}},indent=2))
if win:
 # Record the exact installed MSYS2 package owning each distributed DLL.
 from collect_msys_sources import collect_msys_sources
 collect_msys_sources(Path(os.environ['MSYS2_ROOT']),smr.parent,notices,sources,smrbuild/'_deps')
if runtime_archive is None:runtime_archive=save(ROOT,a.platform,RUNTIME,conda)
for language in ['zh-CN','en']:
 course=BUILD/'editions'/language
 run([conda,'run','--no-capture-output','-p',RUNTIME,'python','scripts/prepare_language.py','--language',language,'--output',course])
 out=ROOT/'.artifacts'/('pmo-'+language+'-'+a.platform)
 run([conda,'run','--no-capture-output','-p',RUNTIME,'python','scripts/assemble.py','--runtime',RUNTIME,'--platform',a.platform,'--language',language,'--course',course,'--output',out,'--smr',smr,'--plink',plink,'--runtime-archive',runtime_archive])
 if not win:
  for dll in smr.parent.glob('*.dylib'):shutil.copy2(dll,out/'resources/course/tools'/a.platform/dll.name)
 shutil.copytree(notices,out/'resources/licenses',dirs_exist_ok=True)
 for name in ['LICENSE','LICENSE-CONTENT.md','THIRD_PARTY_NOTICES.md','DATA_LICENSES.tsv']:
  shutil.copy2(ROOT/name,out/'resources/licenses'/name)
# Sources are separate download assets; they never become student work files.
source_release=ROOT/'.artifacts/sources';source_release.mkdir(parents=True,exist_ok=True)
source_index=[]
for folder in sorted(sources.iterdir()):
 members=[folder]if folder.is_file()else sorted(f for f in folder.rglob('*')if f.is_file())
 groups=[];group=[];total=0
 for f in members:
  if group and total+f.stat().st_size>1500*1024**2:groups.append(group);group=[];total=0
  group.append(f);total+=f.stat().st_size
 if group:groups.append(group)
 for part,group in enumerate(groups,1):
  name=source_release/('sources-'+a.platform+'-'+folder.name+'-'+str(part)+'.zip')
  with zipfile.ZipFile(name,'w',zipfile.ZIP_STORED)as archive:
   for f in group:archive.write(f,f.relative_to(sources).as_posix())
  assert name.stat().st_size<2*1024**3,'An individual source is too large: '+str(name)
  source_index.append({'asset':name.name,'bytes':name.stat().st_size,'sha256':hashlib.sha256(name.read_bytes()).hexdigest(),'members':[f.relative_to(sources).as_posix()for f in group]})
(source_release/('sources-'+a.platform+'-index.json')).write_text(json.dumps(source_index,indent=2))
shutil.rmtree(RUNTIME)  # Tests must not resolve libraries from the build prefix.
print('Both language candidates built:',a.platform,flush=True)
