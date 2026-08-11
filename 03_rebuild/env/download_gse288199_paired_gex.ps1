param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [ValidateSet("Metadata", "Matrix", "All")]
  [string]$Stage = "Metadata",
  [switch]$Force
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$manifestPath = Join-Path $Workspace "03_rebuild\config\GSE288199_PAIRED_GEX_DOWNLOAD_MANIFEST.csv"
$outputRoot = Join-Path $Workspace "00_raw_data\external_validation\GSE288199\paired_gex"
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $outputRoot, $logDir | Out-Null

$rows = Import-Csv -LiteralPath $manifestPath
if ($Stage -eq "Metadata") {
  $rows = @($rows | Where-Object { $_.stage -eq "metadata" })
} elseif ($Stage -eq "Matrix") {
  $rows = @($rows | Where-Object { $_.stage -eq "matrix" })
} else {
  $rows = @($rows)
}

$curl = (Get-Command curl.exe -ErrorAction Stop).Source
$status = [System.Collections.Generic.List[object]]::new()
Write-Host "Manifest: $manifestPath"
Write-Host "Output:   $outputRoot"
Write-Host "Stage:    $Stage"
Write-Host "Files:    $($rows.Count)"

for ($index = 0; $index -lt $rows.Count; $index++) {
  $row = $rows[$index]
  $destination = Join-Path $outputRoot $row.filename
  if ((Test-Path -LiteralPath $destination -PathType Leaf) -and -not $Force) {
    $status.Add([pscustomobject]@{
      hn_id = $row.hn_id
      timepoint = $row.timepoint
      component = $row.component
      filename = $row.filename
      status = "SKIPPED_EXISTS"
      bytes = (Get-Item -LiteralPath $destination).Length
      url = $row.url
    })
    continue
  }
  if ($Force -and (Test-Path -LiteralPath $destination -PathType Leaf)) {
    Remove-Item -LiteralPath $destination -Force
  }
  Write-Host "Downloading [$($index + 1)/$($rows.Count)]: $($row.filename)"
  & $curl --fail --location --retry 5 --retry-all-errors --continue-at - `
    --output $destination $row.url
  if ($LASTEXITCODE -ne 0) {
    throw "Download failed for $($row.filename) with curl exit code $LASTEXITCODE."
  }
  $item = Get-Item -LiteralPath $destination
  if ($item.Length -le 0) {
    throw "Downloaded file is empty: $destination"
  }
  $status.Add([pscustomobject]@{
    hn_id = $row.hn_id
    timepoint = $row.timepoint
    component = $row.component
    filename = $row.filename
    status = "DOWNLOADED"
    bytes = $item.Length
    url = $row.url
  })
}

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$statusPath = Join-Path $logDir "gse288199_paired_gex_download_$stamp.csv"
$latestPath = Join-Path $logDir "gse288199_paired_gex_download_latest.csv"
$status | Export-Csv -LiteralPath $statusPath -NoTypeInformation -Encoding utf8
$status | Export-Csv -LiteralPath $latestPath -NoTypeInformation -Encoding utf8
Write-Host "Download status: $statusPath"
Write-Host "Latest status:   $latestPath"
