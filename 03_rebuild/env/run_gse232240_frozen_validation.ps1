param(
  [string]$Workspace = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
  [switch]$MetadataOnly
)

$ErrorActionPreference = "Stop"
$python = (Get-Command python -ErrorAction Stop).Source
$script = Join-Path $Workspace "03_rebuild\analysis\180_validate_gse232240_frozen_family.py"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss_fff"
$logDir = Join-Path $Workspace "03_rebuild\logs"
New-Item -ItemType Directory -Path $logDir -Force | Out-Null
$stdout = Join-Path $logDir "180_gse232240_frozen_validation_$stamp.out.log"
$stderr = Join-Path $logDir "180_gse232240_frozen_validation_$stamp.err.log"

$arguments = @($script, "--workspace", $Workspace)
if ($MetadataOnly) { $arguments += "--metadata-only" }

Write-Host "Stdout: $stdout"
Write-Host "Stderr: $stderr"
$process = Start-Process -FilePath $python -ArgumentList $arguments -Wait -PassThru -NoNewWindow -RedirectStandardOutput $stdout -RedirectStandardError $stderr
if ($process.ExitCode -ne 0) {
  Get-Content -LiteralPath $stderr -Tail 40 -ErrorAction SilentlyContinue
  throw "GSE232240 frozen validation failed with exit code $($process.ExitCode)."
}

Get-Content -LiteralPath $stdout -Tail 30
Write-Host "GSE232240 frozen validation complete."
