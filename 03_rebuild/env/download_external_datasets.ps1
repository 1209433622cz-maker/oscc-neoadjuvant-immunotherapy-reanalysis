param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Manifest = "",
  [string]$OutputRoot = "",
  [switch]$IncludeP2,
  [switch]$IncludeP3,
  [switch]$MetadataOnly,
  [switch]$Force,
  [string]$OnlyAccession = "",
  [string]$OnlyFileName = "",
  [switch]$NoCurl,
  [int]$Retries = 3
)

$ErrorActionPreference = "Stop"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

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
New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null

$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$statusCsv = Join-Path $logDir "external_download_status_$stamp.csv"
$latestCsv = Join-Path $logDir "external_download_status_latest.csv"
$curlExe = ""
if (-not $NoCurl) {
  $curlCmd = Get-Command curl.exe -ErrorAction SilentlyContinue
  if ($curlCmd) {
    $curlExe = $curlCmd.Source
  }
}

function Get-SelectedRows {
  param([array]$Rows)
  $Rows | Where-Object {
    $include = ($_.include_by_default -match "^(?i:true|1|yes)$")
    if ($IncludeP2 -and $_.priority -eq "P2") { $include = $true }
    if ($IncludeP3 -and $_.priority -eq "P3") { $include = $true }
    if ($MetadataOnly -and $_.download_group -notmatch "metadata") { $include = $false }
    if ($OnlyAccession -ne "" -and $_.accession -ne $OnlyAccession) { $include = $false }
    if ($OnlyFileName -ne "" -and $_.file_name -ne $OnlyFileName) { $include = $false }
    $include
  }
}

function Invoke-CurlDownloadWithResume {
  param(
    [string]$Url,
    [string]$OutFile,
    [string]$CurlExe,
    [int]$Retries
  )

  Write-Host "Downloading with curl resume: $Url"
  $args = @(
    "--location",
    "--fail",
    "--retry", "$Retries",
    "--retry-delay", "10",
    "--retry-all-errors",
    "--connect-timeout", "60",
    "--continue-at", "-",
    "--output", $OutFile,
    $Url
  )
  & $CurlExe @args
  if ($LASTEXITCODE -ne 0) {
    throw "curl.exe failed with exit code $LASTEXITCODE"
  }
}

function Invoke-DownloadWithRetry {
  param(
    [string]$Url,
    [string]$OutFile,
    [int]$Retries
  )

  for ($attempt = 1; $attempt -le $Retries; $attempt++) {
    try {
      Write-Host "Downloading [$attempt/$Retries]: $Url"
      Invoke-WebRequest -Uri $Url -OutFile $OutFile -UseBasicParsing
      return
    } catch {
      if ($attempt -eq $Retries) {
        throw
      }
      Write-Warning "Download failed; retrying in 10 seconds. $($_.Exception.Message)"
      Start-Sleep -Seconds 10
    }
  }
}

$rows = Import-Csv -LiteralPath $Manifest
$selected = @(Get-SelectedRows -Rows $rows)
if ($selected.Count -eq 0) {
  throw "No files selected. Remove -MetadataOnly or add -IncludeP2/-IncludeP3 as needed."
}

Write-Host "Workspace:  $Workspace"
Write-Host "Manifest:   $Manifest"
Write-Host "OutputRoot: $OutputRoot"
Write-Host "Selected files: $($selected.Count)"
Write-Host "IncludeP2: $IncludeP2; IncludeP3: $IncludeP3; MetadataOnly: $MetadataOnly; Force: $Force"
if ($OnlyAccession -ne "") { Write-Host "OnlyAccession: $OnlyAccession" }
if ($OnlyFileName -ne "") { Write-Host "OnlyFileName:  $OnlyFileName" }
if ($curlExe -ne "") {
  Write-Host "Downloader: curl.exe with resume ($curlExe)"
} else {
  Write-Host "Downloader: Invoke-WebRequest"
}

$status = New-Object System.Collections.Generic.List[object]
foreach ($row in $selected) {
  $targetDir = Join-Path $OutputRoot $row.target_subdir
  New-Item -ItemType Directory -Force -Path $targetDir | Out-Null
  $targetFile = Join-Path $targetDir $row.file_name
  $state = "downloaded"
  $errorMessage = ""

  try {
    if ((Test-Path -LiteralPath $targetFile) -and -not $Force) {
      $state = "exists"
      Write-Host "Skipping existing file: $targetFile"
    } else {
      $tmpFile = "$targetFile.partial"
      if ($Force) {
        if (Test-Path -LiteralPath $tmpFile) {
          Remove-Item -LiteralPath $tmpFile -Force
        }
        if (Test-Path -LiteralPath $targetFile) {
          Remove-Item -LiteralPath $targetFile -Force
        }
      }
      if ($curlExe -ne "") {
        Invoke-CurlDownloadWithResume -Url $row.url -OutFile $tmpFile -CurlExe $curlExe -Retries $Retries
      } else {
        if (Test-Path -LiteralPath $tmpFile) {
          Remove-Item -LiteralPath $tmpFile -Force
        }
        Invoke-DownloadWithRetry -Url $row.url -OutFile $tmpFile -Retries $Retries
      }
      Move-Item -LiteralPath $tmpFile -Destination $targetFile -Force
    }

    $fileInfo = Get-Item -LiteralPath $targetFile
    $sha = (Get-FileHash -Algorithm SHA256 -LiteralPath $targetFile).Hash
    $status.Add([pscustomobject]@{
      accession = $row.accession
      priority = $row.priority
      file_label = $row.file_label
      file_name = $row.file_name
      state = $state
      size_bytes = $fileInfo.Length
      sha256 = $sha
      path = $targetFile
      url = $row.url
      error = ""
    })
  } catch {
    $errorMessage = $_.Exception.Message
    Write-Warning "Failed: $($row.file_name): $errorMessage"
    $status.Add([pscustomobject]@{
      accession = $row.accession
      priority = $row.priority
      file_label = $row.file_label
      file_name = $row.file_name
      state = "failed"
      size_bytes = ""
      sha256 = ""
      path = $targetFile
      url = $row.url
      error = $errorMessage
    })
  }
}

$status | Export-Csv -LiteralPath $statusCsv -NoTypeInformation -Encoding UTF8
$status | Export-Csv -LiteralPath $latestCsv -NoTypeInformation -Encoding UTF8

$failed = @($status | Where-Object { $_.state -eq "failed" })
Write-Host "Download status: $statusCsv"
Write-Host "Latest status:   $latestCsv"
if ($failed.Count -gt 0) {
  throw "$($failed.Count) download(s) failed. See status CSV."
}
Write-Host "External dataset download complete."
