param([switch]$Diagnostic, [switch]$ProbeOnly, [switch]$NoDialog)

# Run only in this process; do not modify registry, file associations or policy.
$ErrorActionPreference = 'Stop'
$Root = $PSScriptRoot
$Log = Join-Path $Root 'startup.log'
$SavedPython = Join-Path $Root 'python_path.txt'
$Probe = Join-Path $Root 'src\launch_probe.py'
$Candidates = New-Object 'System.Collections.Generic.List[string]'
$Seen = @{}
$Utf8 = New-Object System.Text.UTF8Encoding($false)

function Write-StartupLog([string]$Message) {
    [System.IO.File]::AppendAllText($Log, $Message + [Environment]::NewLine, $Utf8)
    if ($Diagnostic) { Write-Host $Message }
}

function Add-Python([string]$Value) {
    if ([string]::IsNullOrWhiteSpace($Value)) { return }
    $Path = [Environment]::ExpandEnvironmentVariables($Value.Trim().Trim('"'))
    if ([System.IO.Path]::GetFileName($Path) -ieq 'pythonw.exe') {
        $Path = Join-Path ([System.IO.Path]::GetDirectoryName($Path)) 'python.exe'
    }
    # Avoid launching the Microsoft Store through an execution-alias stub.
    if ($env:LOCALAPPDATA) {
        $Aliases = Join-Path $env:LOCALAPPDATA 'Microsoft\WindowsApps\'
        if ($Path.StartsWith($Aliases, [StringComparison]::OrdinalIgnoreCase)) { return }
    }
    if (-not [System.IO.File]::Exists($Path)) { return }
    $Key = $Path.ToLowerInvariant()
    if (-not $Seen.ContainsKey($Key)) {
        $Seen[$Key] = $true
        [void]$Candidates.Add($Path)
    }
}

function New-PythonProcess([string]$Executable, [string]$Arguments) {
    $Info = New-Object System.Diagnostics.ProcessStartInfo
    $Info.FileName = $Executable
    $Info.Arguments = $Arguments
    $Info.WorkingDirectory = $Root
    $Info.UseShellExecute = $false
    $Info.CreateNoWindow = $true
    $Info.RedirectStandardOutput = $true
    $Info.RedirectStandardError = $true
    $Info.StandardOutputEncoding = [Text.Encoding]::UTF8
    $Info.StandardErrorEncoding = [Text.Encoding]::UTF8
    [void]$Info.EnvironmentVariables.Remove('PYTHONHOME')
    [void]$Info.EnvironmentVariables.Remove('PYTHONPATH')
    $Info.EnvironmentVariables['PYTHONUTF8'] = '1'
    $Info.EnvironmentVariables['PYTHONIOENCODING'] = 'utf-8'
    $Process = New-Object System.Diagnostics.Process
    $Process.StartInfo = $Info
    return $Process
}

function Test-Python([string]$Executable) {
    $Process = $null
    try {
        $Process = New-PythonProcess $Executable ('-I "' + $Probe + '"')
        [void]$Process.Start()
        $OutputTask = $Process.StandardOutput.ReadToEndAsync()
        $ErrorTask = $Process.StandardError.ReadToEndAsync()
        if (-not $Process.WaitForExit(10000)) {
            $Process.Kill()
            Write-StartupLog ('Skipped (probe timeout): ' + $Executable)
            return $null
        }
        $Output = $OutputTask.Result
        $Errors = $ErrorTask.Result
        if ($Process.ExitCode -ne 0) {
            Write-StartupLog ('Skipped: ' + $Executable + ' | ' + $Output.Trim() + ' ' + $Errors.Trim())
            return $null
        }
        $Report = $Output | ConvertFrom-Json
        if (-not $Report.ok -or -not [System.IO.File]::Exists([string]$Report.executable)) {
            Write-StartupLog ('Invalid probe: ' + $Executable)
            return $null
        }
        Write-StartupLog ('Selected Python ' + $Report.version + ', 64-bit: ' + $Report.executable)
        return $Report
    } catch {
        Write-StartupLog ('Skipped: ' + $Executable + ' | ' + $_.Exception.Message)
        return $null
    } finally {
        if ($null -ne $Process) { $Process.Dispose() }
    }
}

try {
    [System.IO.File]::WriteAllText($Log, ('Plasma kinetics startup ' + (Get-Date).ToString('s') + [Environment]::NewLine), $Utf8)
    if (-not [System.IO.File]::Exists($Probe)) { throw 'Extract the entire archive: src\launch_probe.py is missing.' }
    Add-Python $env:PLASMA_PYTHON
    if ([System.IO.File]::Exists($SavedPython)) {
        foreach ($Line in Get-Content -LiteralPath $SavedPython -Encoding UTF8) {
            if ($Line.Trim() -and -not $Line.Trim().StartsWith('#')) { Add-Python $Line; break }
        }
    }
    Add-Python (Join-Path $Root '.venv\Scripts\python.exe')
    if ($env:VIRTUAL_ENV) { Add-Python (Join-Path $env:VIRTUAL_ENV 'Scripts\python.exe') }
    if ($env:CONDA_PREFIX) { Add-Python (Join-Path $env:CONDA_PREFIX 'python.exe') }

    # Standard installs are found even when PATH and .pyw associations are stale.
    $Bases = @()
    if ($env:LOCALAPPDATA) { $Bases += (Join-Path $env:LOCALAPPDATA 'Programs\Python') }
    if ($env:ProgramFiles) { $Bases += $env:ProgramFiles }
    if ($env:SystemDrive) { $Bases += ($env:SystemDrive + '\') }
    foreach ($Base in $Bases) {
        foreach ($Directory in Get-ChildItem -LiteralPath $Base -Directory -Filter 'Python3*' -ErrorAction SilentlyContinue | Sort-Object Name) {
            Add-Python (Join-Path $Directory.FullName 'python.exe')
        }
    }
    foreach ($RegistryRoot in @('HKCU:\Software\Python', 'HKLM:\Software\Python', 'HKLM:\Software\Wow6432Node\Python')) {
        foreach ($Key in Get-Item -Path ($RegistryRoot + '\*\*\InstallPath') -ErrorAction SilentlyContinue) {
            Add-Python ([string]$Key.GetValue('ExecutablePath'))
            $Directory = [string]$Key.GetValue('')
            if ($Directory) { Add-Python (Join-Path $Directory 'python.exe') }
        }
    }
    foreach ($Command in Get-Command python.exe, python3.exe -CommandType Application -All -ErrorAction SilentlyContinue) { Add-Python $Command.Source }
    if ($env:USERPROFILE) {
        foreach ($Folder in @('anaconda3', 'miniconda3')) { Add-Python (Join-Path (Join-Path $env:USERPROFILE $Folder) 'python.exe') }
    }
    # Listing runtimes does not start the launcher's broken default interpreter.
    foreach ($Launcher in Get-Command py.exe -CommandType Application -All -ErrorAction SilentlyContinue) {
        try {
            foreach ($Line in @(& $Launcher.Source -0p 2>$null)) {
                if ([string]$Line -match '(?i)((?:[a-z]:\\|\\\\).+python(?:w)?\.exe)\s*$') { Add-Python $Matches[1] }
            }
        } catch { Write-StartupLog ('Could not enumerate py.exe: ' + $_.Exception.Message) }
    }
    $Selected = $null
    foreach ($Candidate in $Candidates) {
        $Selected = Test-Python $Candidate
        if ($null -ne $Selected) { break }
    }
    if ($null -eq $Selected -and -not $ProbeOnly -and -not $NoDialog) {
        Add-Type -AssemblyName System.Windows.Forms
        $Dialog = New-Object System.Windows.Forms.OpenFileDialog
        $Dialog.Title = 'Choose installed Python 3.11+ (64-bit): python.exe'
        $Dialog.Filter = 'Python executable|python.exe;python3*.exe|Executable files|*.exe'
        try {
            if ($Dialog.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) { $Selected = Test-Python $Dialog.FileName }
        } finally { $Dialog.Dispose() }
    }
    if ($null -eq $Selected) { throw 'No working 64-bit CPython 3.11+ with Tkinter was found. Install Python from python.org or select its python.exe. See startup.log.' }
    if ($ProbeOnly) { $Selected | ConvertTo-Json -Compress; exit 0 }
    try { [System.IO.File]::WriteAllText($SavedPython, ([string]$Selected.executable + [Environment]::NewLine), $Utf8) } catch { Write-StartupLog 'The selected Python path could not be saved; continuing.' }
    $Run = Join-Path $Root 'run.py'
    if (-not [System.IO.File]::Exists($Run)) { throw 'Extract the entire archive: run.py is missing.' }
    $Gui = New-PythonProcess ([string]$Selected.executable) ('-u "' + $Run + '" --gui')
    try {
        [void]$Gui.Start()
        $OutputTask = $Gui.StandardOutput.ReadToEndAsync()
        $ErrorTask = $Gui.StandardError.ReadToEndAsync()
        $Gui.WaitForExit()
        if ($OutputTask.Result) { Write-StartupLog $OutputTask.Result }
        if ($ErrorTask.Result) { Write-StartupLog $ErrorTask.Result }
        if ($Gui.ExitCode -ne 0) { throw ('The program exited with code ' + $Gui.ExitCode + '. See startup.log.') }
        Write-StartupLog 'Program closed normally.'
    } finally { $Gui.Dispose() }
    exit 0
} catch {
    $Message = $_.Exception.Message + [Environment]::NewLine + $Log
    try { Write-StartupLog $Message } catch { Write-Host $Message }
    if (-not $NoDialog -and -not $ProbeOnly) {
        try {
            Add-Type -AssemblyName System.Windows.Forms
            [void][System.Windows.Forms.MessageBox]::Show($Message, 'Plasma kinetics - startup error', 'OK', 'Error')
        } catch { Write-Host $Message }
    }
    exit 1
}
