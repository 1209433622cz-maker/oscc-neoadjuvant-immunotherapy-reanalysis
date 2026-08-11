param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [switch]$RetryKnownBrokenRds
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
if (-not $RetryKnownBrokenRds) {
  $diagnostic = Join-Path $Workspace "03_rebuild\manifests\GSE301741_RDS_ARCHIVE_STREAM_DIAGNOSTIC.md"
  $rawRoute = Join-Path $Workspace "03_rebuild\env\run_gse301741_raw_reconstruction.ps1"
  throw "Deprecated: the deposited RDS fails before memory allocation and should not be retried by lowering a RAM threshold. See $diagnostic. The completed replacement route is $rawRoute."
}

$guarded = Join-Path $Workspace "03_rebuild\env\extract_gse301741_metadata_commit_guarded.ps1"
& powershell -ExecutionPolicy Bypass -File $guarded -RetryKnownBrokenRds
exit $LASTEXITCODE
