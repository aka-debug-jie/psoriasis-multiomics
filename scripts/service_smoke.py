"""Exercise the actual bundled Jupyter server and its R kernel over loopback."""
from pathlib import Path
import json,subprocess,time,urllib.request,uuid,sys,os,datetime,re

def check_service(bundle):
 res=Path(bundle)/'resources';state_file=res/'work/server.json'
 entry=Path(bundle)/('Run course.bat'if sys.platform=='win32'else'Run course.command')
 cmd=['cmd.exe','/c',str(entry),'--no-browser']if sys.platform=='win32'else['bash',str(entry),'--no-browser']
 with(res/'work/service_test.log').open('w',encoding='utf-8')as log:
  process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
  state=None
  def api(path,method='GET',body=None):
   data=None if body is None else json.dumps(body).encode()
   req=urllib.request.Request(state['base']+path,data=data,method=method,headers={'Authorization':'token '+state['token'],'Content-Type':'application/json'})
   with urllib.request.urlopen(req,timeout=10)as response:
    raw=response.read();return json.loads(raw)if raw else None
  try:
   for _ in range(720):
    if process.poll()is not None:raise RuntimeError('Student launcher exited; see service_test.log')
    if state_file.exists():
     state=json.loads(state_file.read_text());break
    time.sleep(.25)
   assert state is not None,'Server did not become ready'
   assert api('/api/status')['kernels']==0
   assert api('/api/contents/course/notebooks/00_start.ipynb?content=0')['type']=='notebook'
   # Reopening must attach to the existing classroom, not spawn another server.
   subprocess.run(cmd,check=True,timeout=30,stdout=log,stderr=log)
   assert json.loads(state_file.read_text())['pid']==state['pid']
   kernel=api('/api/kernels','POST',{'name':'ir-psoriasis'})
   import websocket
   ws=websocket.create_connection(state['base'].replace('http://','ws://')+'/api/kernels/'+kernel['id']+'/channels?token='+state['token'],timeout=30)
   message={'header':{'msg_id':uuid.uuid4().hex,'username':'student','session':uuid.uuid4().hex,'msg_type':'execute_request','version':'5.3','date':datetime.datetime.now(datetime.timezone.utc).isoformat()},'parent_header':{},'metadata':{},'channel':'shell','content':{'code':'1 + 1','silent':False,'store_history':False,'user_expressions':{},'allow_stdin':False,'stop_on_error':True}}
   ws.send(json.dumps(message));found=False;replied=False;deadline=time.monotonic()+45
   ws.settimeout(5)
   while not(found and replied):
    if time.monotonic()>deadline:raise RuntimeError('R kernel response timed out')
    try:received=json.loads(ws.recv())
    except websocket.WebSocketTimeoutException:continue
    kind=received.get('msg_type')or received.get('header',{}).get('msg_type')
    if kind in ('execute_result','display_data'):found='2'in received['content']['data'].get('text/plain','')
    if kind=='execute_reply':
     assert received['content']['status']=='ok';replied=True
    if kind=='error':raise RuntimeError(str(received['content']))
   ws.close();assert found,'R kernel did not return 2'
   api('/api/kernels/'+kernel['id'],'DELETE')
  finally:
   if state is None:
    for name in ('service_test.log','server.log'):
     path=res/'work'/name
     if path.exists():print(name+'\n'+re.sub(r'(token=)[^\s&]+',r'\1[redacted]','\n'.join(path.read_text(errors='replace').splitlines()[-80:])),flush=True)
   if state:
    try:api('/api/shutdown','POST')
    except OSError:pass
   try:process.wait(timeout=25)
   except subprocess.TimeoutExpired:process.terminate();process.wait(timeout=10)
  assert process.returncode==0
 return {'actual_launcher':True,'contents_api':True,'R_kernel_result':2,'repeat_launch_reused_server':True,'shutdown':True}
