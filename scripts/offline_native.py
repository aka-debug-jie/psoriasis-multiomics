"""Run the unpacked course as an ordinary local user, with outbound networking blocked.
Only called on disposable native CI machines; never shipped in the student folder.
"""
from pathlib import Path
import argparse,subprocess,sys,os,secrets,json,base64,shutil,zipfile
p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,required=True);p.add_argument('--platform',required=True);p.add_argument('--smoke-only',action='store_true');p.add_argument('--probe-only',action='store_true');a=p.parse_args()
root=Path(__file__).resolve().parents[1];bundle=a.bundle.resolve();res=bundle/'resources';rt=res/'work/runtime';user='course_student'
def run(args,**kw):return subprocess.run([str(x)for x in args],check=True,**kw)
if sys.platform=='win32':
    user='course_'+secrets.token_hex(4)
    if a.probe_only:
        res.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(bundle/'probe.zip','w')as z:z.writestr('示例.txt','offline-ok')
        long_fixture=Path('\\\\?\\'+str(bundle))/('path-'+'a'*70)/('b'*70)/('c'*70)/'probe.txt'
        long_fixture.parent.mkdir(parents=True,exist_ok=True);long_fixture.write_text('long path permission probe')
    password=secrets.token_urlsafe(24)+'Aa1!';rule='CourseOffline-'+secrets.token_hex(5)
    python=rt/'python.exe'
    ps=lambda x:"'"+str(x).replace("'","''")+"'"
    # The ordinary user performs first extraction as well as computation.
    initialize=(f"Add-Type -AssemblyName System.IO.Compression.FileSystem; [IO.Compression.ZipFile]::ExtractToDirectory({ps(bundle/'probe.zip')},{ps(bundle/'解压 检查')}); if([IO.File]::ReadAllText({ps(bundle/'解压 检查/示例.txt')}) -ne 'offline-ok'){{throw 'Unicode ZIP extraction failed'}}; Remove-Item -LiteralPath {ps(bundle/('path-'+'a'*70))} -Recurse -Force")if a.probe_only else f'& {ps(bundle/"Run course.bat")} --check\nif($LASTEXITCODE -ne 0){{exit $LASTEXITCODE}}'
    compute='exit 0'if a.probe_only else f"$env:COURSE_RUNTIME=(Get-Content -LiteralPath {ps(res/'work/runtime_location.json')} -Raw -Encoding UTF8 | ConvertFrom-Json).path\n"+f'& (Join-Path $env:COURSE_RUNTIME "python.exe") -E -s -X utf8 {ps(root/"scripts/test_bundle.py")} --platform {ps(a.platform)} --inside {ps(bundle)} '+('--smoke-only'if a.smoke_only else'')+'\nexit $LASTEXITCODE'
    child=f'''
$ErrorActionPreference='Stop'
$env:COURSE_OFFLINE_ENFORCED='1'
$env:TEMP=$PSScriptRoot; $env:TMP=$PSScriptRoot; $env:TMPDIR=$PSScriptRoot
$identity=[Security.Principal.WindowsIdentity]::GetCurrent()
$principal=New-Object Security.Principal.WindowsPrincipal($identity)
if($principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)){{throw 'Test account must not be an administrator'}}
$env:PYTHONHOME='C:\\course-nonexistent-python'; $env:PYTHONPATH='C:\\course-nonexistent-site'
{initialize}
$probe=New-Object System.Net.Sockets.TcpClient
try {{
 $connection=$probe.ConnectAsync('1.1.1.1',443)
 if($connection.Wait(3000) -and $probe.Connected){{throw 'PowerShell can reach the external network'}}
}} catch [System.AggregateException] {{}} finally {{$probe.Dispose()}}
{compute}
'''
    child_encoded=base64.b64encode(child.encode('utf-16le')).decode()
    script=f'''
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$user={ps(user)}; $rule={ps(rule)}
$password=ConvertTo-SecureString {ps(password)} -AsPlainText -Force
try {{
 New-LocalUser -Name $user -Password $password -AccountNeverExpires | Out-Null
 Add-LocalGroupMember -SID 'S-1-5-32-545' -Member $user
 & icacls.exe {ps('\\\\?\\'+str(bundle))} /grant "$($user):(OI)(CI)F" /T /Q | Out-Null
 if($LASTEXITCODE -ne 0){{throw 'Cannot grant test folder access'}}
 $remote=@('0.0.0.0-126.255.255.255','128.0.0.0-255.255.255.255','::2-ffff:ffff:ffff:ffff:ffff:ffff:ffff:ffff')
 $sid=(Get-LocalUser -Name $user).SID.Value
 New-NetFirewallRule -DisplayName $rule -Direction Outbound -Action Block -LocalUser "D:(A;;CC;;;$sid)" -RemoteAddress $remote -Profile Any | Out-Null
 $temp=Join-Path $env:SystemRoot ('Temp\\'+$rule); New-Item -ItemType Directory -Path $temp -Force | Out-Null
 & icacls.exe $temp /grant "$($user):(OI)(CI)F" /T /Q | Out-Null
 $env:TEMP=$temp; $env:TMP=$temp; $env:TMPDIR=$temp
 $env:COURSE_OFFLINE_ENFORCED='1'
 $credential=New-Object System.Management.Automation.PSCredential("$env:COMPUTERNAME\\$user",$password)
 New-Item -ItemType Directory -Path {ps(res/'work')} -Force | Out-Null
 $childPath=Join-Path $temp 'run.ps1'
 [IO.File]::WriteAllText($childPath,[Text.Encoding]::Unicode.GetString([Convert]::FromBase64String('{child_encoded}')),(New-Object Text.UTF8Encoding($true)))
 $powershell=Join-Path $env:SystemRoot 'System32\\WindowsPowerShell\\v1.0\\powershell.exe'
 $arguments='-NoProfile -ExecutionPolicy Bypass -File "'+$childPath+'"'
 $process=Start-Process -FilePath $powershell -ArgumentList $arguments -WorkingDirectory {ps(bundle)} -Credential $credential -LoadUserProfile -PassThru -RedirectStandardOutput {ps(res/'work/offline_stdout.log')} -RedirectStandardError {ps(res/'work/offline_stderr.log')}
 $null=$process.Handle; $seen=0; $deadline=(Get-Date).AddMinutes(30)
 do {{
  $done=$process.WaitForExit(10000)
  $lines=@(Get-Content -LiteralPath {ps(res/'work/offline_stdout.log')} -Encoding UTF8)
  if($lines.Count -gt $seen){{$lines[$seen..($lines.Count-1)] | Write-Output; $seen=$lines.Count}}
  if(-not $done -and (Get-Date) -gt $deadline){{Stop-Process -Id $process.Id -Force;throw 'Windows course validation exceeded 30 minutes'}}
 }} while(-not $done)
 if($process.ExitCode -ne 0){{Get-Content -LiteralPath {ps(res/'work/offline_stdout.log')} -Encoding UTF8 -Tail 40;Get-Content -LiteralPath {ps(res/'work/offline_stderr.log')} -Encoding UTF8 -Tail 80;throw "Offline course failed: $($process.ExitCode)"}}
 Start-Sleep -Seconds 2
 $orphans=@(Get-Process -IncludeUserName -ErrorAction SilentlyContinue | Where-Object {{$_.UserName -eq "$env:COMPUTERNAME\\$user" -and $_.ProcessName -in @('R','Rterm','python','pythonw')}})
 if($orphans.Count){{$orphans | Stop-Process -Force;throw 'R or Jupyter processes remained after shutdown'}}
}} finally {{
 Get-NetFirewallRule -DisplayName $rule -ErrorAction SilentlyContinue | Remove-NetFirewallRule
 if($temp -and (Test-Path -LiteralPath $temp)){{Remove-Item -LiteralPath $temp -Recurse -Force -ErrorAction SilentlyContinue}}
 if($sid){{Get-CimInstance Win32_UserProfile -Filter "SID='$sid'" -ErrorAction SilentlyContinue | Remove-CimInstance -ErrorAction SilentlyContinue}}
 if(Get-LocalUser -Name $user -ErrorAction SilentlyContinue){{Remove-LocalUser -Name $user}}
}}
'''
    encoded=base64.b64encode(script.encode('utf-16le')).decode()
    status=subprocess.run(['powershell.exe','-NoProfile','-EncodedCommand',encoded],check=False).returncode
    if status:raise RuntimeError('Windows ordinary-user validation failed; see PowerShell output above')
    if a.probe_only:print('Windows ordinary-user offline probe passed',flush=True)
else:
    python=rt/'bin/python';uid=1501
    listing=subprocess.check_output(['dscl','.','-list','/Users','UniqueID'],text=True)
    used={int(s.split()[-1])for s in listing.splitlines()if s.split()[-1].lstrip('-').isdigit()}
    while uid in used:uid+=1
    anchor='com.apple/course-offline';policy=root/'.build/offline.pf'
    policy.write_text(f'block return out quick on ! lo0 proto {{ tcp udp }} from any to any user {uid}\n')
    created=False
    try:
        run(['sudo','dscl','.','-create','/Users/'+user]);created=True
        for key,value in [('UniqueID',str(uid)),('PrimaryGroupID','20'),('UserShell','/bin/bash'),('NFSHomeDirectory','/private/tmp/course_student')]:run(['sudo','dscl','.','-create','/Users/'+user,key,value])
        run(['sudo','mkdir','-p','/private/tmp/course_student']);run(['sudo','chown','-R',f'{uid}:20',bundle,'/private/tmp/course_student'])
        run(['sudo','pfctl','-f','/etc/pf.conf'])
        enabled=subprocess.run(['sudo','pfctl','-e'],capture_output=True,text=True)
        if enabled.returncode!=0 and 'already enabled'not in enabled.stderr.lower():raise RuntimeError(enabled.stderr)
        run(['sudo','pfctl','-a',anchor,'-f',policy])
        run(['sudo','-u',user,'-H','env','PYTHONHOME=/course-nonexistent-python','PYTHONPATH=/course-nonexistent-site','bash',bundle/'Run course.command','--check'])
        run(['sudo','-u',user,'-H','env','COURSE_OFFLINE_ENFORCED=1',python,'-E','-s',root/'scripts/test_bundle.py','--platform',a.platform,'--inside',bundle]+(['--smoke-only']if a.smoke_only else[]))
    finally:
        subprocess.run(['sudo','pfctl','-a',anchor,'-F','rules'],check=False)
        if created:subprocess.run(['sudo','dscl','.','-delete','/Users/'+user],check=False)
        subprocess.run(['sudo','chown','-R',f'{os.getuid()}:{os.getgid()}',bundle],check=False)
