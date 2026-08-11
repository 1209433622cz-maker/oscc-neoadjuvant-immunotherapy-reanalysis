param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\runtimes\R-4.3.3\bin\Rscript.exe",
  [switch]$RetryKnownBrokenRds,
  [double]$MinimumFreePhysicalGB = 20,
  [double]$MinimumFreeVirtualGB = 30,
  [double]$AbortFreePhysicalGB = 2.5,
  [double]$AbortFreeVirtualGB = 8,
  [double]$MaximumRPrivateGB = 34,
  [int]$PollSeconds = 5,
  [int]$TimeoutMinutes = 90
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$Rscript = (Resolve-Path -LiteralPath $Rscript).Path
$env:GSE200996_WORKSPACE = $Workspace
$env:R_LIBS_USER = "C:\Program Files\R\R-4.3.3\library"
$env:R_MAX_VSIZE = "38Gb"
$env:R_MAX_NUM_DLLS = "200"
$env:TMPDIR = Join-Path $Workspace "03_rebuild\tmp\gse301741_rds"
$env:TEMP = $env:TMPDIR
$env:TMP = $env:TMPDIR
New-Item -ItemType Directory -Force -Path $env:TMPDIR | Out-Null

$rds = Join-Path $Workspace "00_raw_data\external_validation\GSE301741\GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds"
$script = Join-Path $Workspace "03_rebuild\analysis\26_extract_gse301741_rds_metadata.R"
$logDir = Join-Path $Workspace "03_rebuild\logs"
$manifestDir = Join-Path $Workspace "03_rebuild\manifests"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
New-Item -ItemType Directory -Force -Path $manifestDir | Out-Null

if (-not $RetryKnownBrokenRds) {
  $diagnostic = Join-Path $manifestDir "GSE301741_RDS_ARCHIVE_STREAM_DIAGNOSTIC.md"
  throw "The deposited RDS is known to fail identically in R 4.3.3 and R 4.6.0 before memory allocation. Do not lower memory guards. See $diagnostic and use run_gse301741_raw_reconstruction.ps1."
}

if (-not (Test-Path -LiteralPath $rds)) {
  throw "GSE301741 RDS not found: $rds"
}

$expectedSha256 = "B31D782EE3BD6144065C403C2ACE32C7FAA9DA452F8F5627FF8FFBE82186DCDC"
$observedSha256 = (Get-FileHash -LiteralPath $rds -Algorithm SHA256).Hash
if ($observedSha256 -ne $expectedSha256) {
  throw "SHA256 mismatch. Expected $expectedSha256; observed $observedSha256"
}

$stream = [System.IO.File]::OpenRead($rds)
try {
  $header = New-Object byte[] 2
  [void]$stream.Read($header, 0, 2)
} finally {
  $stream.Dispose()
}
$headerHex = [BitConverter]::ToString($header)
if ($headerHex -ne "58-0A") {
  Write-Warning "RDS is not the expected uncompressed XDR stream. Header: $headerHex"
}

$os = Get-CimInstance Win32_OperatingSystem
$freePhysicalGB = [math]::Round($os.FreePhysicalMemory / 1MB, 2)
$freeVirtualGB = [math]::Round($os.FreeVirtualMemory / 1MB, 2)
$totalPhysicalGB = [math]::Round($os.TotalVisibleMemorySize / 1MB, 2)
$totalVirtualGB = [math]::Round($os.TotalVirtualMemorySize / 1MB, 2)

Write-Host "RDS SHA256: $observedSha256"
Write-Host "RDS header: $headerHex (58-0A means uncompressed XDR)"
Write-Host "Physical memory: $freePhysicalGB GB free / $totalPhysicalGB GB total"
Write-Host "Virtual memory:  $freeVirtualGB GB free / $totalVirtualGB GB total"

if ($freePhysicalGB -lt $MinimumFreePhysicalGB) {
  throw "Free physical memory is below $MinimumFreePhysicalGB GB. Close other analyses before retrying."
}
if ($freeVirtualGB -lt $MinimumFreeVirtualGB) {
  throw "Free virtual memory is below $MinimumFreeVirtualGB GB. A larger page file or a high-memory machine is required."
}

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$stdout = Join-Path $logDir "26_extract_gse301741_metadata_guarded_$stamp.out.log"
$stderr = Join-Path $logDir "26_extract_gse301741_metadata_guarded_$stamp.err.log"
$monitor = Join-Path $logDir "26_extract_gse301741_metadata_guarded_$stamp.monitor.csv"
"timestamp,free_physical_gb,free_virtual_gb,r_working_set_gb,r_private_gb,action" | Set-Content -LiteralPath $monitor -Encoding ascii

$arguments = @("--vanilla", $script)
$process = Start-Process `
  -FilePath $Rscript `
  -ArgumentList $arguments `
  -WorkingDirectory $Workspace `
  -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr `
  -WindowStyle Hidden `
  -PassThru

try {
  $process.PriorityClass = "BelowNormal"
} catch {
  Write-Warning "Could not lower R process priority: $($_.Exception.Message)"
}

$started = Get-Date
$aborted = $false
$abortReason = ""
while (-not $process.HasExited) {
  Start-Sleep -Seconds $PollSeconds
  $process.Refresh()
  $os = Get-CimInstance Win32_OperatingSystem
  $freePhysicalGB = [math]::Round($os.FreePhysicalMemory / 1MB, 3)
  $freeVirtualGB = [math]::Round($os.FreeVirtualMemory / 1MB, 3)
  $workingGB = [math]::Round($process.WorkingSet64 / 1GB, 3)
  $privateGB = [math]::Round($process.PrivateMemorySize64 / 1GB, 3)
  $elapsed = (Get-Date) - $started
  $action = "continue"

  if ($freePhysicalGB -lt $AbortFreePhysicalGB) {
    $aborted = $true
    $abortReason = "Free physical memory fell below $AbortFreePhysicalGB GB."
  } elseif ($freeVirtualGB -lt $AbortFreeVirtualGB) {
    $aborted = $true
    $abortReason = "Free virtual memory fell below $AbortFreeVirtualGB GB."
  } elseif ($privateGB -gt $MaximumRPrivateGB) {
    $aborted = $true
    $abortReason = "R private memory exceeded $MaximumRPrivateGB GB."
  } elseif ($elapsed.TotalMinutes -gt $TimeoutMinutes) {
    $aborted = $true
    $abortReason = "Runtime exceeded $TimeoutMinutes minutes."
  }

  if ($aborted) {
    $action = "terminate"
  }
  "$(Get-Date -Format o),$freePhysicalGB,$freeVirtualGB,$workingGB,$privateGB,$action" | Add-Content -LiteralPath $monitor -Encoding ascii

  if ($aborted) {
    Stop-Process -Id $process.Id -Force
    break
  }
}

$process.WaitForExit()
$process.Refresh()
Write-Host "Stdout: $stdout"
Write-Host "Stderr: $stderr"
Write-Host "Monitor: $monitor"

if ($aborted) {
  throw "Guarded extraction was stopped safely: $abortReason Use the high-memory handoff route."
}
$success = Join-Path $manifestDir "GSE301741_RDS_METADATA_EXTRACTION_SUCCESS.txt"
if (-not (Test-Path -LiteralPath $success)) {
  throw "R extraction did not produce the success marker. See stderr log: $stderr"
}

Write-Host "GSE301741 guarded metadata extraction complete."
Write-Host "Success marker: $success"
