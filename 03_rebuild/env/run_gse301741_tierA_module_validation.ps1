param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$env:GSE200996_WORKSPACE = $Workspace

$script = Join-Path $Workspace "03_rebuild\analysis\39_validate_gse301741_tierA_modules.R"
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$out = Join-Path $logDir "39_gse301741_tierA_module_validation_$stamp.out.log"
$err = Join-Path $logDir "39_gse301741_tierA_module_validation_$stamp.err.log"

$proc = Start-Process -FilePath $Rscript `
  -ArgumentList @($script) `
  -WorkingDirectory $Workspace `
  -RedirectStandardOutput $out `
  -RedirectStandardError $err `
  -WindowStyle Hidden `
  -Wait `
  -PassThru

Write-Host "Stdout: $out"
Write-Host "Stderr: $err"
if ($proc.ExitCode -ne 0) {
  throw "GSE301741 Tier A module validation failed with exit code $($proc.ExitCode)."
}

Write-Host "GSE301741 Tier A module validation complete."
