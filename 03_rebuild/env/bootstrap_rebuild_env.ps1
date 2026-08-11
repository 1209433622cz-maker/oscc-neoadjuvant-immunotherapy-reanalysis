param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
  [string]$Python = "C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe",
  [switch]$Optional
)

$ErrorActionPreference = "Stop"

$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$env:GSE200996_WORKSPACE = $Workspace

$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

if (-not (Test-Path -LiteralPath $Rscript)) {
  throw "Rscript not found: $Rscript"
}

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$installOut = Join-Path $logDir "env_install_R_$stamp.out.log"
$installErr = Join-Path $logDir "env_install_R_$stamp.err.log"
$checkOut = Join-Path $logDir "env_check_R_$stamp.out.log"
$checkErr = Join-Path $logDir "env_check_R_$stamp.err.log"

$installScript = Join-Path $Workspace "03_rebuild\env\install_R_packages.R"
$checkScript = Join-Path $Workspace "03_rebuild\env\check_rebuild_environment.R"
$pythonCheck = Join-Path $Workspace "03_rebuild\env\check_python_environment.py"

Write-Host "Workspace: $Workspace"
Write-Host "Rscript:   $Rscript"
Write-Host "Optional R packages: $Optional"

$rArgs = @($installScript)
if ($Optional) {
  $rArgs += "--optional"
}

$proc = Start-Process -FilePath $Rscript `
  -ArgumentList $rArgs `
  -WorkingDirectory $Workspace `
  -RedirectStandardOutput $installOut `
  -RedirectStandardError $installErr `
  -WindowStyle Hidden `
  -Wait `
  -PassThru

if ($proc.ExitCode -ne 0) {
  throw "R package installation failed with exit code $($proc.ExitCode). Logs: $installOut ; $installErr"
}

$proc = Start-Process -FilePath $Rscript `
  -ArgumentList @($checkScript) `
  -WorkingDirectory $Workspace `
  -RedirectStandardOutput $checkOut `
  -RedirectStandardError $checkErr `
  -WindowStyle Hidden `
  -Wait `
  -PassThru

if ($proc.ExitCode -ne 0) {
  throw "R environment check failed with exit code $($proc.ExitCode). Logs: $checkOut ; $checkErr"
}

if (Test-Path -LiteralPath $Python) {
  $pyReq = Join-Path $Workspace "03_rebuild\env\requirements-python.txt"
  $pyPipOut = Join-Path $logDir "env_install_python_$stamp.out.log"
  $pyPipErr = Join-Path $logDir "env_install_python_$stamp.err.log"
  if (Test-Path -LiteralPath $pyReq) {
    $proc = Start-Process -FilePath $Python `
      -ArgumentList @("-m", "pip", "install", "-r", $pyReq) `
      -WorkingDirectory $Workspace `
      -RedirectStandardOutput $pyPipOut `
      -RedirectStandardError $pyPipErr `
      -WindowStyle Hidden `
      -Wait `
      -PassThru
    if ($proc.ExitCode -ne 0) {
      Write-Warning "Python package installation reported issues. Logs: $pyPipOut ; $pyPipErr"
    }
  }

  $pyOut = Join-Path $logDir "env_check_python_$stamp.out.log"
  $pyErr = Join-Path $logDir "env_check_python_$stamp.err.log"
  $proc = Start-Process -FilePath $Python `
    -ArgumentList @($pythonCheck) `
    -WorkingDirectory $Workspace `
    -RedirectStandardOutput $pyOut `
    -RedirectStandardError $pyErr `
    -WindowStyle Hidden `
    -Wait `
    -PassThru
  if ($proc.ExitCode -ne 0) {
    Write-Warning "Python environment check reported missing packages. Logs: $pyOut ; $pyErr"
  }
} else {
  Write-Warning "Python not found at $Python; skipped Python check."
}

Write-Host "Environment bootstrap complete."
Write-Host "Latest R status:      $logDir\ENVIRONMENT_STATUS_LATEST.md"
Write-Host "Latest Python status: $logDir\PYTHON_ENVIRONMENT_STATUS_LATEST.md"
