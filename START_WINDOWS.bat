@echo off
setlocal EnableExtensions DisableDelayedExpansion
set "PLASMA_PS=%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe"
if not exist "%PLASMA_PS%" goto no_powershell
if not exist "%~dp0START_WINDOWS.ps1" goto incomplete_archive
if /I "%~1"=="debug" goto debug
start "" "%PLASMA_PS%" -NoLogo -NoProfile -STA -WindowStyle Hidden -ExecutionPolicy Bypass -File "%~dp0START_WINDOWS.ps1"
if errorlevel 1 goto launch_failed
exit /b 0

:debug
"%PLASMA_PS%" -NoLogo -NoProfile -STA -ExecutionPolicy Bypass -File "%~dp0START_WINDOWS.ps1" -Diagnostic
set "PLASMA_BOOT_EXIT=%ERRORLEVEL%"
if not "%PLASMA_BOOT_EXIT%"=="0" pause
exit /b %PLASMA_BOOT_EXIT%

:no_powershell
echo Windows PowerShell was not found. Run: python run.py --gui
goto launch_failed
:incomplete_archive
echo Extract the entire archive first. START_WINDOWS.ps1 is missing.
:launch_failed
echo Startup failed. Open a command prompt here and run START_WINDOWS.bat debug
pause
exit /b 1
