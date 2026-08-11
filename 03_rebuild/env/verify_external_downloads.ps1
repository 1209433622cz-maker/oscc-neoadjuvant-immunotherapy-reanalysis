param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Manifest = "",
  [string]$OutputRoot = "",
  [switch]$IncludeP2,
  [switch]$IncludeP3,
  [switch]$MetadataOnly,
  [string]$OnlyAccession = "",
  [string]$OnlyFileName = ""
)

$ErrorActionPreference = "Stop"

$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
if ($Manifest -eq "") {
  $Manifest = Join-Path $Workspace "03_rebuild\config\external_download_manifest.csv"
}
if ($OutputRoot -eq "") {
  $OutputRoot = Join-Path $Workspace "00_raw_data\external_validation"
}

if (-not (Test-Path -LiteralPath $Manifest)) {
  throw "Download manifest not found: $Manifest"
}

$OutputRoot = [System.IO.Path]::GetFullPath($OutputRoot)
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$outCsv = Join-Path $logDir "external_download_verification_$stamp.csv"
$latestCsv = Join-Path $logDir "external_download_verification_latest.csv"
$outMd = Join-Path $logDir "EXTERNAL_DOWNLOAD_VERIFICATION_LATEST.md"

$rows = Import-Csv -LiteralPath $Manifest
$selected = @($rows | Where-Object {
  $include = ($_.include_by_default -match "^(?i:true|1|yes)$")
  if ($IncludeP2 -and $_.priority -eq "P2") { $include = $true }
  if ($IncludeP3 -and $_.priority -eq "P3") { $include = $true }
  if ($MetadataOnly -and $_.download_group -notmatch "metadata") { $include = $false }
  if ($OnlyAccession -ne "" -and $_.accession -ne $OnlyAccession) { $include = $false }
  if ($OnlyFileName -ne "" -and $_.file_name -ne $OnlyFileName) { $include = $false }
  $include
})

$status = foreach ($row in $selected) {
  $targetFile = Join-Path (Join-Path $OutputRoot $row.target_subdir) $row.file_name
  $exists = Test-Path -LiteralPath $targetFile
  $size = ""
  $sha = ""
  if ($exists) {
    $info = Get-Item -LiteralPath $targetFile
    $size = $info.Length
    $sha = (Get-FileHash -Algorithm SHA256 -LiteralPath $targetFile).Hash
  }
  [pscustomobject]@{
    accession = $row.accession
    priority = $row.priority
    file_label = $row.file_label
    file_name = $row.file_name
    expected_size = $row.reported_size
    exists = $exists
    zero_byte = ($exists -and $size -eq 0)
    size_bytes = $size
    sha256 = $sha
    path = $targetFile
    role = $row.role
  }
}

$status | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding UTF8
$status | Export-Csv -LiteralPath $latestCsv -NoTypeInformation -Encoding UTF8

$missing = @($status | Where-Object { -not $_.exists })
$zeroByte = @($status | Where-Object { $_.zero_byte })
$lines = @(
  "# External Download Verification",
  "",
  "Checked: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')",
  "",
  "- Manifest: $Manifest",
  "- Output root: $OutputRoot",
  "- Selected files: $($selected.Count)",
  "- Missing files: $($missing.Count)",
  "- Zero-byte files: $($zeroByte.Count)",
  "",
  "| Accession | Priority | File | Exists | Zero byte | Size bytes | SHA256 |",
  "|---|---|---|---:|---:|---:|---|"
)
foreach ($item in $status) {
  $lines += "| $($item.accession) | $($item.priority) | $($item.file_name) | $($item.exists) | $($item.zero_byte) | $($item.size_bytes) | $($item.sha256) |"
}
$lines | Set-Content -LiteralPath $outMd -Encoding UTF8

Write-Host "Verification CSV: $outCsv"
Write-Host "Latest CSV:       $latestCsv"
Write-Host "Latest Markdown:  $outMd"
if ($missing.Count -gt 0) {
  Write-Warning "$($missing.Count) selected file(s) are missing."
  exit 1
}
if ($zeroByte.Count -gt 0) {
  Write-Warning "$($zeroByte.Count) selected file(s) are zero bytes. Review metadata endpoints before using them."
}
Write-Host "All selected external files are present."
