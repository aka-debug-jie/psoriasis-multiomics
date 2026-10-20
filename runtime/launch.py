"""Start the bundled, local-only R classroom. Uses the bundled Python."""
from pathlib import Path
import argparse,json,os,secrets,socket,subprocess,sys,time,urllib.request,urllib.error,webbrowser,tempfile,re,shlex

RES=Path(__file__).resolve().parent
WORK=RES/'work';RUNTIME=Path(os.environ.get('COURSE_RUNTIME',str(WORK/'runtime')));COURSE=RES/'course'
META=json.loads((RES/'bundle.json').read_text(encoding='utf-8')) if (RES/'bundle.json').exists() else {}
LANG=META.get('language','zh-CN')
def tr(zh,en):return en if LANG=='en' else zh

def request(url,token,method='GET'):
    req=urllib.request.Request(url,headers={'Authorization':'token '+token},method=method)
    with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req,timeout=2) as response:return response.read()

def fix_r_shell_paths():
    # Conda's R wrapper stores unquoted build paths. Quote the unpacked location.
    home=RUNTIME/'lib/R'
    values={'R_HOME_DIR':home,'R_SHARE_DIR':home/'share','R_INCLUDE_DIR':home/'include','R_DOC_DIR':home/'doc'}
    seen=set()
    for path in [RUNTIME/'bin/R',home/'bin/R']:
        if not path.exists():continue
        resolved=path.resolve()
        if resolved in seen or not resolved.is_relative_to(RUNTIME.resolve()):continue
        seen.add(resolved)
        source=resolved.read_text(encoding='utf-8')
        if not source.startswith('#!'):continue
        updated=source
        for key,value in values.items():
            updated=re.sub('^'+key+'=.*$',lambda match:key+'='+shlex.quote(str(value)),updated,flags=re.M)
        if updated!=source:resolved.write_text(updated,encoding='utf-8')

def make_environment():
    WORK.mkdir(exist_ok=True)
    env=os.environ.copy()
    for name in ('PYTHONPATH','PYTHONHOME','R_LIBS','R_LIBS_USER','R_LIBS_SITE','R_HOME','RENV_PROJECT','CONDA_PREFIX','TEACHING_SMR','TEACHING_PLINK','JUPYTER_CONFIG_PATH','GH_TOKEN','GITHUB_TOKEN','ACTIONS_RUNTIME_TOKEN','ACTIONS_ID_TOKEN_REQUEST_TOKEN'):
        env.pop(name,None)
    win=sys.platform=='win32'
    if win:env.update(PYTHONUTF8='1',PYTHONIOENCODING='utf-8')
    if not win:fix_r_shell_paths()
    binary=RUNTIME/'R/bin/x64/Rterm.exe' if win else RUNTIME/'bin/R'
    if not binary.is_file():raise RuntimeError(tr('内置 R 文件缺失，请重新解压完整课程包。','Bundled R is missing. Extract a fresh complete course package.'))
    lib=RUNTIME/'course-library'
    if not lib.is_dir():raise RuntimeError(tr('课程包库缺失，请重新解压完整课程包。','The bundled R library is missing. Extract a fresh complete package.'))
    config=WORK/'jupyter';(config/'kernels/ir-psoriasis').mkdir(parents=True,exist_ok=True)
    for sub in ('tmp','results','jupyter-runtime','jupyter-config'):(WORK/sub).mkdir(exist_ok=True)
    env.update(TEACHING_LANG=LANG,TEACHING_BUNDLED='1',TEACHING_PROFILE='interactive',TEACHING_FRONTEND='jupyter',
        TEACHING_LIBRARY=str(lib),TEACHING_NOTEBOOK_OUTPUT='../work/results',
        JUPYTER_PATH=str(config),JUPYTER_DATA_DIR=str(config),JUPYTER_CONFIG_DIR=str(WORK/'jupyter-config'),
        JUPYTER_RUNTIME_DIR=str(WORK/'jupyter-runtime'),IPYTHONDIR=str(WORK/'ipython'),
        PYTHONNOUSERSITE='1',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',
        TMPDIR=tempfile.gettempdir() if win else '/tmp',TEMP=tempfile.gettempdir() if win else '/tmp',TMP=tempfile.gettempdir() if win else '/tmp')
    # Never consult a student's globally installed R or Python.
    system_path=os.path.join(os.environ.get('SystemRoot','C:\\Windows'),'System32') if win else '/usr/bin:/bin:/usr/sbin:/sbin'
    bins=[RUNTIME,RUNTIME/'Scripts',RUNTIME/'Library/bin',RUNTIME/'R/bin/x64']if win else[RUNTIME/'bin']
    env['PATH']=os.pathsep.join(str(p)for p in bins)+os.pathsep+system_path
    expression='.libPaths(c('+json.dumps(str(lib).replace('\\','/'))+',.Library)); IRkernel::main()'
    spec={'argv':[str(binary),'--vanilla','--slave','-e',expression,'--args','{connection_file}'],
          'display_name':tr('R · 多组学课程','R · Multi-omics course'),'language':'R'}
    (config/'kernels/ir-psoriasis/kernel.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
    return env,binary

def startup_lock():
    handle=(WORK/'startup.lock').open('a+b')
    if handle.tell()==0:handle.write(b'0');handle.flush()
    handle.seek(0)
    if sys.platform=='win32':
        import msvcrt
        deadline=time.monotonic()+150
        while True:
            try:msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1);break
            except OSError:
                if time.monotonic()>deadline:handle.close();raise RuntimeError(tr('另一窗口正在启动，请稍候重试。','Another window is starting the course. Please retry shortly.'))
                time.sleep(.25)
    else:
        import fcntl
        fcntl.flock(handle,fcntl.LOCK_EX)
    return handle

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--no-browser',action='store_true');args=p.parse_args()
    WORK.mkdir(exist_ok=True)
    lock=startup_lock()
    try:env,rbin=make_environment()
    except BaseException:lock.close();raise
    if args.check:
        lock.close()
        subprocess.run([str(rbin),'--vanilla','--slave','-e',
          '.libPaths(c(Sys.getenv("TEACHING_LIBRARY"),.Library)); source("R/environment.R"); s<-course_package_status(frontend="both");stopifnot(all(s$available));source("R/tool_support.R");tool_path("smr");tool_path("plink");print(s)'],cwd=COURSE,env=env,check=True)
        print(tr('课程环境与外部工具已就绪。','The course environment and external tools are ready.'));return
    state_path=WORK/'server.json'
    if state_path.exists():
        try:
            state=json.loads(state_path.read_text())
            body=json.loads(request(state['base']+'/api/status',state['token']))
            if 'started' in body and state.get('root')==str(RES):
                if not args.no_browser:webbrowser.open(state['url'])
                lock.close()
                print(tr('课程已经运行，已打开现有课堂。','The course is already running. Opened the existing session.'));return
        except PermissionError:lock.close();raise
        except (OSError,ValueError,KeyError):pass
    token=secrets.token_hex(24)
    with socket.socket()as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    base=f'http://127.0.0.1:{port}';url=base+'/lab/tree/course/notebooks/00_start.ipynb?token='+token
    env['JUPYTER_TOKEN']=token
    # Keep configuration, logs, and student outputs together inside course resources.
    config={'ServerApp':{'root_dir':str(RES),'ip':'127.0.0.1','port':port,'port_retries':0,'open_browser':False},
            'LanguageServerManager':{'autodetect':False},
            'LabApp':{'extension_manager':'readonly','news_url':None,'check_for_updates_class':'jupyterlab.handlers.announcements.NeverCheckForUpdate'}}
    if sys.platform=='win32':
        config['ServerApp']['cookie_secret_file']=str(WORK/'jupyter-runtime'/('cookie-'+secrets.token_hex(12)+'.key'))
    config_file=WORK/'jupyter-config/jupyter_config.json';config_file.write_text(json.dumps(config),encoding='utf-8')
    with (WORK/'server.log').open('a',encoding='utf-8')as log:
        process=subprocess.Popen([sys.executable,'-m','jupyterlab','--no-browser','--config',str(config_file)],env=env,cwd=RES,stdout=log,stderr=subprocess.STDOUT)
        try:
            for _ in range(480):
                if process.poll() is not None:raise RuntimeError(tr('课堂未启动成功，请查看 resources/work/server.log。','The server did not start. See resources/work/server.log.'))
                try:request(base+'/api/status',token);break
                except OSError:time.sleep(.25)
            else:raise RuntimeError(tr('课堂启动超时，请查看 resources/work/server.log。','Server startup timed out. See resources/work/server.log.'))
            temp_state=WORK/('server-'+secrets.token_hex(6)+'.tmp')
            temp_state.write_text(json.dumps({'base':base,'url':url,'token':token,'root':str(RES),'pid':process.pid}),encoding='utf-8')
            if os.name!='nt':temp_state.chmod(0o600)
            os.replace(temp_state,state_path)
            lock.close()
            if not args.no_browser:webbrowser.open(url)
            print(tr('课堂已打开。学习结束后，在 Jupyter 的 File 菜单选择 Shut Down。','The course is open. Choose File → Shut Down in Jupyter when finished.'),flush=True)
            print(tr('也可以在此窗口按 Ctrl+C 结束课堂。','You can also press Ctrl+C in this window to stop the course.'),flush=True)
            process.wait()
        except KeyboardInterrupt:
            try:request(base+'/api/shutdown',token,'POST')
            except OSError:pass
        finally:
            lock.close()
            if process.poll()is None:
                try:process.wait(timeout=15)
                except subprocess.TimeoutExpired:process.terminate();process.wait(timeout=10)
            cleanup_lock=startup_lock()
            try:
                try:current=json.loads(state_path.read_text())
                except (FileNotFoundError,ValueError):current={}
                if current.get('pid')==process.pid and current.get('token')==token:state_path.unlink(missing_ok=True)
            finally:cleanup_lock.close()

if __name__=='__main__':
    try:main()
    except Exception as exc:print(str(exc).replace(str(RES),'%COURSE%').replace(str(Path.home()),'%USERPROFILE%'),file=sys.stderr);sys.exit(1)
