param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\runtimes\R-4.3.3\bin\Rscript.exe",
  [double]$MinimumFreeMemoryGB = 8
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$Rscript = (Resolve-Path -LiteralPath $Rscript).Path
$env:GSE200996_WORKSPACE = $Workspace
$env:R_LIBS_USER = "C:\Program Files\R\R-4.3.3\library"
$env:R_MAX_VSIZE = "24Gb"

$os = Get-CimInstance Win32_OperatingSystem
$freeGB = [math]::Round($os.FreePhysicalMemory / 1MB, 2)
Write-Host "Free physical memory: $freeGB GB"
if ($freeGB -lt $MinimumFreeMemoryGB) {
  throw "Free memory is below $MinimumFreeMemoryGB GB. Wait for other analyses to finish."
}

$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"

function Invoke-RStep {
  param([string]$Step, [string]$ScriptName)
  $script = Join-Path $Workspace "03_rebuild\analysis\$ScriptName"
  $stdout = Join-Path $logDir "${Step}_${stamp}.out.log"
  $stderr = Join-Path $logDir "${Step}_${stamp}.err.log"
  Write-Host "Starting $Step"
  Write-Host "Script: $script"
  Write-Host "Stdout: $stdout"
  Write-Host "Stderr: $stderr"
  $process = Start-Process `
    -FilePath $Rscript `
    -ArgumentList @("--vanilla", $script) `
    -WorkingDirectory $Workspace `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -WindowStyle Hidden `
    -PassThru `
    -Wait
  if ($process.ExitCode -ne 0) {
    throw "$Step failed with exit code $($process.ExitCode). See $stderr"
  }
  Write-Host "Completed $Step"
}

Invoke-RStep -Step "44_gse301741_raw_reconstruction" `
  -ScriptName "44_rebuild_gse301741_raw_cell_metadata.R"
Invoke-RStep -Step "45_gse301741_lineage_aware_validation" `
  -ScriptName "45_validate_gse301741_lineage_aware.R"
Invoke-RStep -Step "46_gse301741_raw_reconstruction_audit" `
  -ScriptName "46_audit_gse301741_raw_reconstruction.R"

Write-Host "GSE301741 RAW reconstruction, lineage-aware validation and audit complete."
