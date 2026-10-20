"""Fetch the reviewed teaching snapshot at build time, never during student startup."""
from pathlib import Path,PurePosixPath
import hashlib,json,subprocess,zipfile,os
ROOT=Path(__file__).resolve().parents[1]
meta=json.loads((ROOT/'config/data_asset.json').read_text())
folder=ROOT/'.build/inputs';folder.mkdir(parents=True,exist_ok=True)
archive=folder/meta['asset']
if not archive.exists():subprocess.run(['gh','release','download','v4.0.0','--repo','aka-debug-jie/psoriasis-multiomics','--pattern',meta['asset'],'--dir',str(folder)],check=True)
assert hashlib.sha256(archive.read_bytes()).hexdigest()==meta['sha256'],'Data asset checksum mismatch'
with zipfile.ZipFile(archive)as z:
 for member in z.infolist():
  path=PurePosixPath(member.filename)
  assert path.parts and path.parts[0]=='data'and '..'not in path.parts and not path.is_absolute(),member.filename
  assert (member.external_attr>>16)&0o170000!=0o120000,'Unexpected symbolic link'
 z.extractall(ROOT/'course')
print('Reviewed course input snapshot verified and extracted.')
