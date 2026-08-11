param(
  [string]$Workspace = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
  [switch]$Force,
  [switch]$ArchiveOnly
)

$ErrorActionPreference = "Stop"
$outDir = Join-Path $Workspace "00_raw_data\external_validation\GSE232240"
$archive = Join-Path $outDir "GSE232240_RAW.tar"
$partial = "$archive.partial"
$expectedBytes = 61460480L
$url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232240/suppl/GSE232240_RAW.tar"
$countName = "GSM7324294_Count_data_IMCISION.txt.gz"
$metaName = "GSM7324295_Meta_data_IMCISION.txt.gz"

New-Item -ItemType Directory -Path $outDir -Force | Out-Null

if ($Force) {
  Remove-Item -LiteralPath $partial -Force -ErrorAction SilentlyContinue
}

$archiveReady = (Test-Path -LiteralPath $archive) -and ((Get-Item -LiteralPath $archive).Length -eq $expectedBytes)
if (-not $archiveReady) {
  if ((Test-Path -LiteralPath $archive) -and ((Get-Item -LiteralPath $archive).Length -gt 0)) {
    throw "An archive with an unexpected size already exists: $archive. Move it aside or rerun with a clean target."
  }

  Write-Host "Downloading GSE232240 archive with IPv4 and resume support."
  Write-Host "Target: $partial"
  & curl.exe -4 -L --fail --retry 30 --retry-delay 10 --retry-all-errors --connect-timeout 30 --continue-at - --output $partial $url
  if ($LASTEXITCODE -ne 0) {
    throw "curl failed with exit code $LASTEXITCODE. Keep the .partial file and rerun the same command to resume."
  }
  if ((Get-Item -LiteralPath $partial).Length -ne $expectedBytes) {
    throw "Downloaded size is $((Get-Item -LiteralPath $partial).Length), expected $expectedBytes bytes. Rerun to resume."
  }
  Move-Item -LiteralPath $partial -Destination $archive -Force
}

$listing = & tar.exe -tf $archive
if ($LASTEXITCODE -ne 0) {
  throw "tar integrity listing failed for $archive"
}
foreach ($required in @($countName, $metaName)) {
  if (-not ($listing | Where-Object { $_ -eq $required -or $_.EndsWith("/$required") })) {
    throw "Required archive member not found: $required"
  }
}

if (-not $ArchiveOnly) {
  foreach ($member in @($countName, $metaName)) {
    $target = Join-Path $outDir $member
    if ($Force -or -not (Test-Path -LiteralPath $target) -or (Get-Item -LiteralPath $target).Length -eq 0) {
      & tar.exe -xf $archive -C $outDir $member
      if ($LASTEXITCODE -ne 0) {
        throw "Failed to extract $member"
      }
    }
  }
}

$python = (Get-Command python -ErrorAction Stop).Source
$gzipFiles = @($countName, $metaName) | ForEach-Object { Join-Path $outDir $_ } | Where-Object { Test-Path -LiteralPath $_ }
foreach ($file in $gzipFiles) {
  & $python -c "import gzip,sys; f=gzip.open(sys.argv[1],'rb'); [None for _ in iter(lambda:f.read(1024*1024),b'')]; f.close(); print('gzip OK:',sys.argv[1])" $file
  if ($LASTEXITCODE -ne 0) {
    throw "gzip integrity check failed: $file"
  }
}

Write-Host "GSE232240 acquisition verified."
Write-Host "Archive: $archive"
Write-Host "Next: powershell -ExecutionPolicy Bypass -File $Workspace\03_rebuild\env\run_gse232240_frozen_validation.ps1"
