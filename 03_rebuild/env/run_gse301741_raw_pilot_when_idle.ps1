param(
  [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
  [int]$MaxScRNA = 4,
  [int]$MinimumFreeMemoryGB = 8
)

$ErrorActionPreference = "Stop"
$Workspace = (Resolve-Path -LiteralPath $Workspace).Path
$env:GSE200996_WORKSPACE = $Workspace

$os = Get-CimInstance Win32_OperatingSystem
$freeGb = [math]::Round(($os.FreePhysicalMemory * 1024) / 1GB, 2)
Write-Host "Free physical memory: $freeGb GB"
if ($freeGb -lt $MinimumFreeMemoryGB) {
  throw "Free memory is below $MinimumFreeMemoryGB GB. Defer the GSE301741 raw-route pilot."
}

$script = Join-Path $Workspace "03_rebuild\analysis\28_pilot_gse301741_raw_route.R"
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$out = Join-Path $logDir "28_pilot_gse301741_raw_route_$stamp.out.log"
$err = Join-Path $logDir "28_pilot_gse301741_raw_route_$stamp.err.log"

$args = @($script, "--max-scrna=$MaxScRNA")
$proc = Start-Process -FilePath $Rscript `
  -ArgumentList $args `
  -WorkingDirectory $Workspace `
  -RedirectStandardOutput $out `
  -RedirectStandardError $err `
  -WindowStyle Hidden `
  -Wait `
  -PassThru

Write-Host "Stdout: $out"
Write-Host "Stderr: $err"
if ($proc.ExitCode -ne 0) {
  throw "GSE301741 raw-route pilot failed with exit code $($proc.ExitCode)."
}
Write-Host "GSE301741 raw-route pilot complete."

