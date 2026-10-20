"""Verify the reviewed static site asset before a Pages deployment."""
from pathlib import Path,PurePosixPath
import sys,json,hashlib,zipfile
archive=Path(sys.argv[1]);target=Path(sys.argv[2]);root=Path(__file__).resolve().parents[1]
meta=json.loads((root/'config/site_asset.json').read_text())
assert hashlib.sha256(archive.read_bytes()).hexdigest()==meta['sha256'],'Site asset checksum mismatch'
with zipfile.ZipFile(archive)as z:
 for item in z.infolist():
  path=PurePosixPath(item.filename)
  assert not path.is_absolute()and '..'not in path.parts
  assert (item.external_attr>>16)&0o170000!=0o120000
 z.extractall(target)
print('Verified static site extracted.')
