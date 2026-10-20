@echo off
chcp 65001 >nul
setlocal
set "COURSE_RESOURCES=%~dp0resources"
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy Bypass -File "%COURSE_RESOURCES%\bootstrap.ps1" %*
if not errorlevel 1 exit /b 0
if "%~1"=="--check" exit /b 1
pause
exit /b 1
