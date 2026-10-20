param([Parameter(ValueFromRemainingArguments=$true)][string[]]$LaunchArgs)
$ErrorActionPreference = 'Stop'
$res = $PSScriptRoot
$english = $false
$lock = $null
$stage = $null
function Say($zh, $en) { if ($english) { Write-Host $en } else { Write-Host $zh } }
function Fail($code, $zh, $en) {
    if ($english) { throw "[$code] $en" } else { throw "[$code] $zh" }
}
function Writable($folder) {
    try {
        [IO.Directory]::CreateDirectory($folder) | Out-Null
        $probe = Join-Path $folder ('.write-' + [guid]::NewGuid().ToString('N'))
        [IO.File]::WriteAllText($probe, 'ok')
        [IO.File]::Delete($probe)
        return $true
    } catch { return $false }
}
try {
    $meta = Get-Content -LiteralPath (Join-Path $res 'bundle.json') -Raw -Encoding UTF8 | ConvertFrom-Json
    $english = $meta.language -eq 'en'
    if ($env:PROCESSOR_ARCHITECTURE -eq 'ARM64' -or $env:PROCESSOR_ARCHITEW6432 -eq 'ARM64' -or -not [Environment]::Is64BitOperatingSystem) {
        Fail 'ARCH' '此包适用于 Windows x64，请下载匹配电脑的版本。' 'This package needs Windows x64. Download the package for your computer.'
    }
    $work = Join-Path $res 'work'
    if (-not (Writable $work)) {
        Fail 'WRITE' '课程文件夹不可写，请将完整文件夹移到自己的本地目录。' 'The course folder is not writable. Move the complete folder to a local folder you own.'
    }
    $longestCourse = (Get-ChildItem -LiteralPath (Join-Path $res 'course') -File -Recurse | ForEach-Object { $_.FullName.Length } | Measure-Object -Maximum).Maximum
    if ($longestCourse -gt 240 -or (Join-Path $work 'snapshots\07_single_cell\20261020120000_practice_001\pseudobulk_candidates.tsv').Length -gt 240) {
        Fail 'PATH' '课程数据路径过长。请把完整文件夹移到 D:\course 或其他较短的可写目录，再打开入口。' 'Course data paths are too long. Move the complete folder to D:\course or another short writable folder, then open the launcher.'
    }
    $seed = Join-Path $res 'runtime\seed.zip'
    Say '正在检查内置环境，首次打开需要稍等片刻。' 'Checking the bundled environment. First launch may take a moment.'
    $digest = (Get-FileHash -LiteralPath $seed -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($digest -ne $meta.runtime_sha256) {
        Fail 'PACKAGE' '环境文件校验失败，请重新下载并完整解压课程包。' 'The runtime checksum does not match. Download and extract a fresh package.'
    }
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $archive = [IO.Compression.ZipFile]::OpenRead($seed)
    try {
        $longest = ($archive.Entries | ForEach-Object { $_.FullName.Length } | Measure-Object -Maximum).Maximum
        $expanded = ($archive.Entries | Measure-Object -Property Length -Sum).Sum
        $count = $archive.Entries.Count
    } finally { $archive.Dispose() }
    $key = $digest.Substring(0,16)
    $candidates = @((Join-Path $work 'runtime'))
    if ($env:LOCALAPPDATA) { $candidates += Join-Path $env:LOCALAPPDATA "PMO\$key\r" }
    $candidates += Join-Path ([IO.Path]::GetTempPath()) "PMO\$key\r"
    $dst = $null
    $required = [double]$expanded + $count * 4096.0 + 512MB
    foreach ($candidate in ($candidates | Select-Object -Unique)) {
        $parent = Split-Path $candidate
        if ([Math]::Max($candidate.Length, $parent.Length + 14) + 1 + $longest -gt 240) { continue }
        if (-not (Writable $parent)) { continue }
        $candidateReady = Join-Path $candidate '.pmo-ready.json'
        $candidateOwned = Join-Path $candidate '.pmo-owned.json'
        if (Test-Path -LiteralPath $candidateReady) {
            try { $r = Get-Content -LiteralPath $candidateReady -Raw -Encoding UTF8 | ConvertFrom-Json }
            catch { continue }
            if ($r.sha256 -ne $digest -or $r.path -ne $candidate -or $r.protocol -ne 1) { continue }
            $dst = $candidate; break
        }
        if (Test-Path -LiteralPath $candidate) {
            if (-not (Test-Path -LiteralPath $candidateOwned)) { continue }
            try { $owned = Get-Content -LiteralPath $candidateOwned -Raw -Encoding UTF8 | ConvertFrom-Json }
            catch { continue }
            if ($owned.sha256 -ne $digest -or $owned.path -ne $candidate) { continue }
            if ((Get-Item -LiteralPath $candidate).Attributes -band [IO.FileAttributes]::ReparsePoint) { continue }
        }
        $drive = [IO.DriveInfo]::new([IO.Path]::GetPathRoot($candidate))
        if ($drive.AvailableFreeSpace -lt $required) { continue }
        $dst = $candidate; break
    }
    if (-not $dst) {
        Fail 'CACHE' '没有可用的短环境目录，请把完整课程文件夹移到 D:\course 后重试。' 'No short writable runtime location is available. Move the complete course folder to D:\course and retry.'
    }
    $parent = Split-Path $dst
    $deadline = [DateTime]::UtcNow.AddMinutes(5)
    while (-not $lock) {
        try { $lock = [IO.File]::Open((Join-Path $parent '.pmo-init.lock'), 'OpenOrCreate', 'ReadWrite', 'None') }
        catch {
            if (($_.Exception.HResult -band 0xFFFF) -ne 32) { throw }
            if ([DateTime]::UtcNow -gt $deadline) { Fail 'BUSY' '另一窗口仍在准备环境，请等待其完成后重试。' 'Another window is preparing the runtime. Wait for it to finish and retry.' }
            Start-Sleep -Milliseconds 500
        }
    }
    $readyPath = Join-Path $dst '.pmo-ready.json'
    $ready = $false
    if (Test-Path -LiteralPath $readyPath) {
        try {
            $readyInfo = Get-Content -LiteralPath $readyPath -Raw -Encoding UTF8 | ConvertFrom-Json
            $ready = $readyInfo.sha256 -eq $digest -and $readyInfo.path -eq $dst -and $readyInfo.protocol -eq 1
        } catch { $ready = $false }
    }
    if (-not $ready) {
        $drive = [IO.DriveInfo]::new([IO.Path]::GetPathRoot($dst))
        if ($drive.AvailableFreeSpace -lt $required) { Fail 'SPACE' '环境所在磁盘空间不足，请至少释放 2 GB 后重试。' 'Not enough space for the runtime. Free at least 2 GB on its drive and retry.' }
        Say '正在准备内置环境，无需安装其他软件。' 'Preparing the bundled environment. No additional installation is needed.'
        $stage = Join-Path $parent ('tmp-' + [guid]::NewGuid().ToString('N').Substring(0,8))
        [IO.Compression.ZipFile]::ExtractToDirectory($seed, $stage)
        @{sha256=$digest; path=$dst; protocol=1} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $stage '.pmo-owned.json') -Encoding UTF8
        if (Test-Path -LiteralPath $dst) {
            if (Test-Path -LiteralPath $readyPath) { Fail 'CACHE_ID' '缓存版本不匹配，请使用新解压的课程文件夹。' 'The runtime cache has a different identity. Use a freshly extracted course folder.' }
            Remove-Item -LiteralPath $dst -Recurse -Force
        }
        Move-Item -LiteralPath $stage -Destination $dst
        $stage = $null
        $env:PATH = "$dst;$(Join-Path $dst 'Library\bin');$(Join-Path $dst 'Scripts');$env:SystemRoot\System32"
        & (Join-Path $dst 'python.exe') -E -s -X utf8 (Join-Path $dst 'Scripts\conda-unpack-script.py')
        if ($LASTEXITCODE -ne 0) { Fail 'UNPACK' '环境准备未完成，请关闭窗口后重新打开课程。' 'Runtime preparation failed. Close this window and open the course again.' }
        & (Join-Path $dst 'python.exe') -E -s -c "import sys, jupyterlab; assert sys.version_info[:2] == (3,12)"
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $dst 'R\bin\x64\Rterm.exe'))) { Fail 'RUNTIME' '内置软件不完整，请重新下载完整课程包。' 'Bundled software is incomplete. Download a fresh course package.' }
        $readyTemp = Join-Path $dst '.pmo-ready.tmp'
        @{sha256=$digest; path=$dst; protocol=1} | ConvertTo-Json | Set-Content -LiteralPath $readyTemp -Encoding UTF8
        Move-Item -LiteralPath $readyTemp -Destination $readyPath
    }
    $lock.Dispose(); $lock = $null
    @{path=$dst; sha256=$digest} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $work 'runtime_location.json') -Encoding UTF8
    $env:COURSE_RUNTIME = $dst
    $env:PATH = "$dst;$(Join-Path $dst 'Library\bin');$(Join-Path $dst 'Scripts');$env:SystemRoot\System32"
    & (Join-Path $dst 'python.exe') -E -s -X utf8 (Join-Path $res 'launch.py') @LaunchArgs
    exit $LASTEXITCODE
} catch {
    $message = $_.Exception.Message.Replace($res, '%COURSE%')
    if ($env:USERPROFILE) { $message = $message.Replace($env:USERPROFILE, '%USERPROFILE%') }
    Write-Host $message
    Say '课程未启动。请打开 resources/help.html，按错误编号排查。' 'The course did not start. Open resources/help.html and check the error code.'
    exit 1
} finally {
    if ($lock) { $lock.Dispose() }
    if ($stage -and (Test-Path -LiteralPath $stage)) { Remove-Item -LiteralPath $stage -Recurse -Force -ErrorAction SilentlyContinue }
}
