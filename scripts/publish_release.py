"""Upload verified assets to a draft. Publishing the release is a separate step."""
from pathlib import Path
import subprocess,sys,os,json,argparse
from runtime_cache import gh
p=argparse.ArgumentParser();p.add_argument('--platform',required=True);p.add_argument('--language',choices=['zh-CN','en'],required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1];os.chdir(root)
subprocess.run([sys.executable,'scripts/finalize.py','--platform',a.platform,'--language',a.language],check=True)
tag='v4.0.0';view=gh('release','view',tag,'--json','isDraft',check=False)
assert view.returncode==0 and json.loads(view.stdout)['isDraft'],'Create the reviewed draft release first'
files=[root/'.artifacts/release'/('pmo-'+a.language+'-'+a.platform+'.zip'),root/'.artifacts/release'/(a.language+'-'+a.platform+'-manifest.json')]
gh('release','upload',tag,*files,'--clobber')
