param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Python = "D:\bioinfor\python.exe"
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
  $Python = (Get-Command python -ErrorAction Stop).Source
}

$script = Join-Path $Workspace "03_rebuild\analysis\120_validate_locked_family_gse195832_bulk.py"
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$out = Join-Path $logDir "120_gse195832_locked_family_$stamp.out.log"
$err = Join-Path $logDir "120_gse195832_locked_family_$stamp.err.log"

Write-Host "Python: $Python"
Write-Host "Script: $script"
Write-Host "Stdout: $out"
Write-Host "Stderr: $err"

$proc = Start-Process -FilePath $Python `
  -ArgumentList @($script) `
  -WorkingDirectory $Workspace `
  -RedirectStandardOutput $out `
  -RedirectStandardError $err `
  -WindowStyle Hidden `
  -Wait `
  -PassThru

if ($proc.ExitCode -ne 0) {
  throw "GSE195832 locked-family validation failed with exit code $($proc.ExitCode). See $err"
}
Write-Host "GSE195832 locked-family validation complete."
