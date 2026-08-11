param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
  [int]$MinimumFreeMemoryGB = 45,
  [switch]$ExtractRdsMetadata,
  [switch]$OverrideMemoryCheck
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$env:GSE200996_WORKSPACE = $Workspace

if ($ExtractRdsMetadata) {
  $extractScript = Join-Path $Workspace "03_rebuild\env\extract_gse301741_metadata_when_idle.ps1"
  $extractArgs = @(
    "-ExecutionPolicy", "Bypass",
    "-File", $extractScript,
    "-Workspace", $Workspace,
    "-Rscript", $Rscript,
    "-MinimumFreeMemoryGB", $MinimumFreeMemoryGB
  )
  if ($OverrideMemoryCheck) {
    $extractArgs += "-OverrideMemoryCheck"
  }
  & powershell @extractArgs
}

$pythonScript = Join-Path $Workspace "03_rebuild\analysis\37_recover_gse301741_response_labels.py"
& python $pythonScript --workspace $Workspace
if ($LASTEXITCODE -ne 0) {
  throw "GSE301741 response-label recovery audit failed."
}

Write-Host "GSE301741 response-label recovery audit complete."
