param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
  [int]$MinimumFreeMemoryGB = 18
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$env:GSE200996_WORKSPACE = $Workspace

$os = Get-CimInstance Win32_OperatingSystem
$freeGb = [math]::Round(($os.FreePhysicalMemory * 1024) / 1GB, 2)
$totalGb = [math]::Round(($os.TotalVisibleMemorySize * 1024) / 1GB, 2)
Write-Host "Free physical memory: $freeGb GB"
Write-Host "Total visible memory: $totalGb GB"
if ($freeGb -lt $MinimumFreeMemoryGB) {
  throw "Free memory is below $MinimumFreeMemoryGB GB. Defer the discovery cohort-sensitivity job."
}

$script = Join-Path $Workspace "03_rebuild\analysis\43_discovery_cohort_adjusted_pseudobulk_sensitivity.R"
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$out = Join-Path $logDir "43_discovery_cohort_sensitivity_$stamp.out.log"
$err = Join-Path $logDir "43_discovery_cohort_sensitivity_$stamp.err.log"

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
  throw "Discovery cohort-sensitivity analysis failed with exit code $($proc.ExitCode)."
}
Write-Host "Discovery cohort-sensitivity analysis complete."
