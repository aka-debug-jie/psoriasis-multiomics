"""Execute selected course chapters in an already configured R kernel."""
from pathlib import Path
import argparse,json,time,os
import nbformat
from nbclient import NotebookClient
p=argparse.ArgumentParser();p.add_argument('course',type=Path);p.add_argument('--chapters',default=','.join(f'{i:02}' for i in range(13)));p.add_argument('--kernel',default='ir-psoriasis');p.add_argument('--report',type=Path,required=True);a=p.parse_args()
chapters=set(a.chapters.split(','));report=[]
os.environ['TEACHING_LANG']=(a.course/'language.txt').read_text().strip()
for file in sorted((a.course/'notebooks').glob('*.ipynb')):
 if file.name[:2]not in chapters:continue
 start=time.monotonic();print('Execute',file.name,flush=True)
 nb=nbformat.read(file,as_version=4)
 NotebookClient(nb,timeout=1200,kernel_name=a.kernel,resources={'metadata':{'path':str(file.parent.resolve())}},allow_errors=False).execute()
 nbformat.write(nb,file)
 report.append({'chapter':file.stem,'seconds':round(time.monotonic()-start,2),'code_cells':sum(c.cell_type=='code'for c in nb.cells),'status':'passed'})
 a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2))
 print('Passed',file.name,flush=True)
