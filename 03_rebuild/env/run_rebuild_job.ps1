param(
  [ValidateSet("audit", "full", "summary", "loo", "all-light", "all-heavy")]
  [string]$Job = "summary",
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
  [string]$UseCachedObj = "true"
)

$ErrorActionPreference = "Stop"

function Convert-ToBooleanFlag {
  param([string]$Value)
  $text = $Value.Trim().ToLowerInvariant()
  switch -Regex ($text) {
    "^(true|t|1|yes|y)$" { return $true }
    "^(false|f|0|no|n)$" { return $false }
    default {
      throw "Invalid UseCachedObj value '$Value'. Use true/false, 1/0, or yes/no."
    }
  }
}

$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$UseCachedObjBool = Convert-ToBooleanFlag $UseCachedObj

$env:GSE200996_WORKSPACE = $Workspace
$env:GSE200996_BASEDIR = (Join-Path $Workspace "03_rebuild") -replace "\\", "/"
$env:GSE200996_RAWDIR = (Join-Path $Workspace "00_raw_data\GSE200996_RAW") -replace "\\", "/"
$env:GSE200996_META = (Join-Path $Workspace "00_raw_data\GSE200996_metadata\GSE200996_CD45.tumor.single.cell.meta.data.txt.gz") -replace "\\", "/"
$env:GSE200996_USE_CACHED_OBJ = if ($UseCachedObjBool) { "TRUE" } else { "FALSE" }
Write-Host "Use cached Seurat object: $UseCachedObjBool"

if (-not (Test-Path -LiteralPath $Rscript)) {
  throw "Rscript not found: $Rscript"
}

$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

function Invoke-RJob {
  param(
    [string]$Name,
    [string]$Script
  )
  if (-not (Test-Path -LiteralPath $Script)) {
    throw "Script not found: $Script"
  }
  $stamp = Get-Date -Format "yyyyMMdd_HHmmss_fff"
  $stdout = Join-Path $logDir "$Name`_$stamp.out.log"
  $stderr = Join-Path $logDir "$Name`_$stamp.err.log"
  Write-Host "Starting $Name"
  Write-Host "Script: $Script"
  Write-Host "Stdout: $stdout"
  Write-Host "Stderr: $stderr"

  $proc = Start-Process -FilePath $Rscript `
    -ArgumentList @($Script) `
    -WorkingDirectory $Workspace `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -WindowStyle Hidden `
    -Wait `
    -PassThru

  if ($proc.ExitCode -ne 0) {
    throw "$Name failed with exit code $($proc.ExitCode). Logs: $stdout ; $stderr"
  }
  Write-Host "Completed $Name"
}

$auditScript = Join-Path $Workspace "03_rebuild\analysis\01_data_audit.R"
$fullScript = Join-Path $Workspace "03_rebuild\analysis\run_pipeline_rebuild.R"
$summaryScript = Join-Path $Workspace "03_rebuild\analysis\03_summarize_rebuild_results.R"
$looScript = Join-Path $Workspace "03_rebuild\analysis\04_leave_one_patient_diagnostics.R"

switch ($Job) {
  "audit" {
    Invoke-RJob -Name "01_data_audit" -Script $auditScript
  }
  "full" {
    Invoke-RJob -Name "02_full_reanalysis" -Script $fullScript
  }
  "summary" {
    Invoke-RJob -Name "03_summary" -Script $summaryScript
  }
  "loo" {
    Invoke-RJob -Name "04_leave_one_patient" -Script $looScript
  }
  "all-light" {
    Invoke-RJob -Name "01_data_audit" -Script $auditScript
    Invoke-RJob -Name "03_summary" -Script $summaryScript
  }
  "all-heavy" {
    Invoke-RJob -Name "01_data_audit" -Script $auditScript
    Invoke-RJob -Name "02_full_reanalysis" -Script $fullScript
    Invoke-RJob -Name "04_leave_one_patient" -Script $looScript
    Invoke-RJob -Name "03_summary" -Script $summaryScript
  }
}

Write-Host "Job complete: $Job"
