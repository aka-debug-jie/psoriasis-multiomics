"""Record owners and source-only packages for DLLs copied from the MSYS2 toolchain."""
from pathlib import Path
import hashlib,json,re,shlex,subprocess,urllib.request,urllib.parse

def collect_msys_sources(msys_root,dll_directory,notices,sources,dependency_build):
 bash=msys_root/'usr/bin/bash.exe';records=[];packages={}
 def command(args):return subprocess.check_output([str(bash),'-lc',' '.join(shlex.quote(str(x))for x in args)],text=True).strip()
 for dll in sorted(dll_directory.glob('*.dll')):
  digest=hashlib.sha256(dll.read_bytes()).hexdigest()
  built=[p for p in dependency_build.rglob(dll.name)if p.is_file()and hashlib.sha256(p.read_bytes()).hexdigest()==digest]
  if built:
   records.append({'dll':dll.name,'sha256':digest,'owner':'SMR build dependency',
                   'build_artifact':built[0].relative_to(dependency_build).as_posix(),
                   'source_archive':'SMR-1.4.3-complete-source.tar.gz'})
   continue
  installed_file=msys_root/'ucrt64/bin'/dll.name
  assert installed_file.is_file()and hashlib.sha256(installed_file.read_bytes()).hexdigest()==digest, 'Unresolved DLL origin: '+dll.name
  installed='/ucrt64/bin/'+dll.name
  owner=command(['pacman','-Qqo',installed])
  identity=command(['pacman','-Q',owner]);name,version=identity.split(None,1)
  info=command(['pacman','-Qi',owner])
  if name not in packages:
   url='https://packages.msys2.org/packages/'+name
   with urllib.request.urlopen(url,timeout=60)as response:page=response.read().decode()
   links=[urllib.parse.urljoin(url,s)for s in re.findall(r'href=["\']([^"\']+\.src\.tar\.(?:zst|gz|xz))["\']',page)]
   links=[s for s in links if version in urllib.parse.unquote(s)]
   item={'package':name,'version':version,'official_page':url,'source_status':'source_package_unresolved','package_info':info}
   if links:
    source_url=links[0];dest=sources/'msys2'/Path(urllib.parse.urlsplit(source_url).path).name;dest.parent.mkdir(exist_ok=True)
    if not dest.exists():urllib.request.urlretrieve(source_url,dest)
    item.update(source_status='source_package_captured_review_upstream_sources',source_url=source_url,source_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),source_archive=dest.relative_to(sources).as_posix())
   packages[name]=item
  records.append({'dll':dll.name,'sha256':digest,'owner':name,'version':version})
 (notices/'msys2-components.json').write_text(json.dumps({'libraries':records,'packages':list(packages.values()),'review_status':'pending_manual_review'},indent=2)+'\n')
