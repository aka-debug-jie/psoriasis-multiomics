"""Pack the freshly built native runtime and provide a quiet GitHub CLI wrapper."""
from pathlib import Path
import os,subprocess,json,hashlib,sys,shutil

def gh(*args,check=True):
 cli=shutil.which('gh')
 if not cli and sys.platform=='win32':cli=str(Path(os.environ['ProgramFiles'])/'GitHub CLI/gh.exe')
 return subprocess.run([cli,*map(str,args)],check=check,capture_output=True,text=True)

def restore(root,platform,prefix):
 # Public builds never read historical private runtime caches.
 return None

def save(root,platform,prefix,conda):
 ext='.zip'if platform.startswith('windows')else'.tar.gz'
 folder=root/'.build/runtime-archives';folder.mkdir(parents=True,exist_ok=True);archive=folder/(platform+ext)
 subprocess.run([str(conda),'run','--no-capture-output','-p',str(prefix),'conda-pack','-p',str(prefix),'-o',str(archive),'--format','zip'if ext=='.zip'else'tar.gz','--force'],check=True)
 return archive
