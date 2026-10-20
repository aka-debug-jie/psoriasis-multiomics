"""Execute selected R notebooks in fresh kernels and report actual outcomes."""
import argparse,json,os,time,traceback
from pathlib import Path
import nbformat
from nbclient import NotebookClient
p=argparse.ArgumentParser();p.add_argument('--lessons',nargs='+',default=[f'{i:02d}'for i in range(1,13)]);p.add_argument('--kernel',default='ir-psoriasis');p.add_argument('--output',required=True);p.add_argument('--save',action='store_true');args=p.parse_args()
root=Path(__file__).resolve().parents[1];out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
os.environ['JUPYTER_PATH']=str(root/'.jupyter')+os.pathsep+os.environ.get('JUPYTER_PATH','')
os.environ['TEACHING_NOTEBOOK_OUTPUT']=str(out/'results');os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
status=[]
for item in args.lessons:
 cid=f'{int(item):02d}';start=time.monotonic()
 try:
  path=next((root/'notebooks').glob(cid+'_*.ipynb'))
  nb=nbformat.read(path,as_version=4)
  NotebookClient(nb,timeout=1200,kernel_name=args.kernel,resources={'metadata':{'path':str(root/'notebooks')}},allow_errors=False,store_widget_state=False).execute()
  if args.save:nbformat.write(nb,path)
  status.append(dict(lesson=cid,frontend='jupyter',status='passed',seconds=round(time.monotonic()-start,2)))
 except Exception as e:
  (out/(cid+'.log')).write_text(traceback.format_exc())
  text=str(e);blocked=any(k in text for k in ['本课需要匹配系统','本课缺少','未找到','适配 R'])
  status.append(dict(lesson=cid,frontend='jupyter',status='unavailable'if blocked else 'failed',seconds=round(time.monotonic()-start,2)))
 (out/'status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
print(json.dumps(status,ensure_ascii=False,indent=2))
if any(s['status']!='passed'for s in status):raise SystemExit(1)
