param(
    [string]$Workspace = "H:\SCI2\OSCC-GSE200996-2025.12",
    [double]$MinimumFreeMemoryGB = 6,
    [int]$IndependentSeed = 20260810
)

$ErrorActionPreference = "Stop"
$workspacePath = [System.IO.Path]::GetFullPath($Workspace)
$cleanRoot = [System.IO.Path]::GetFullPath((Join-Path $workspacePath "03_rebuild\cleanroom\gse232240"))
$inputRoot = [System.IO.Path]::GetFullPath((Join-Path $cleanRoot "input"))
$resultRoot = [System.IO.Path]::GetFullPath((Join-Path $cleanRoot "results"))
$auditRoot = [System.IO.Path]::GetFullPath((Join-Path $cleanRoot "audit"))
foreach ($path in @($inputRoot, $resultRoot, $auditRoot)) {
    if (-not $path.StartsWith($cleanRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Clean-room path escapes the intended root: $path"
    }
}

$memory = Get-CimInstance Win32_OperatingSystem
$freeGB = [math]::Round($memory.FreePhysicalMemory / 1MB, 2)
Write-Host "Free physical memory: $freeGB GB"
if ($freeGB -lt $MinimumFreeMemoryGB) {
    throw "Free memory is below $MinimumFreeMemoryGB GB. Wait before starting the clean-room run."
}

$sourceFiles = [ordered]@{
    "GSM7324294_Count_data_IMCISION.txt.gz" = Join-Path $workspacePath "00_raw_data\external_validation\GSE232240\GSM7324294_Count_data_IMCISION.txt.gz"
    "GSM7324295_Meta_data_IMCISION.txt.gz" = Join-Path $workspacePath "00_raw_data\external_validation\GSE232240\GSM7324295_Meta_data_IMCISION.txt.gz"
    "gse232240_validation.json" = Join-Path $workspacePath "03_rebuild\config\gse232240_validation.json"
    "locked_family_gene_set_manifest.csv" = Join-Path $workspacePath "03_rebuild\config\locked_family_gene_set_manifest.csv"
}

New-Item -ItemType Directory -Path $inputRoot -Force | Out-Null
foreach ($entry in $sourceFiles.GetEnumerator()) {
    if (-not (Test-Path -LiteralPath $entry.Value)) {
        throw "Missing clean-room input: $($entry.Value)"
    }
    $destination = Join-Path $inputRoot $entry.Key
    Copy-Item -LiteralPath $entry.Value -Destination $destination -Force
    $sourceHash = (Get-FileHash -LiteralPath $entry.Value -Algorithm SHA256).Hash
    $destinationHash = (Get-FileHash -LiteralPath $destination -Algorithm SHA256).Hash
    if ($sourceHash -ne $destinationHash) {
        throw "Staged input hash mismatch: $($entry.Key)"
    }
}

foreach ($path in @($resultRoot, $auditRoot)) {
    if (Test-Path -LiteralPath $path) {
        Remove-Item -LiteralPath $path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $path -Force | Out-Null
}

$python = (Get-Command python -ErrorAction Stop).Source
$computeScript = Join-Path $workspacePath "03_rebuild\analysis\184_gse232240_cleanroom_recompute.py"
$compareScript = Join-Path $workspacePath "03_rebuild\analysis\185_compare_gse232240_cleanroom.py"
$baselineRoot = Join-Path $workspacePath "03_rebuild\results\external_validation\GSE232240"
$requiredBaseline = Join-Path $baselineRoot "GSE232240_pseudobulk_group_totals.csv"
if (-not (Test-Path -LiteralPath $requiredBaseline)) {
    throw "The frozen baseline group-total audit asset is missing: $requiredBaseline"
}

Write-Host "Starting isolated clean-room computation."
& $python $computeScript --input-dir $inputRoot --output-dir $resultRoot --independent-seed $IndependentSeed
if ($LASTEXITCODE -ne 0) {
    throw "Independent clean-room computation failed with exit code $LASTEXITCODE."
}

Write-Host "Computation complete. Unblinding once for baseline comparison."
& $python $compareScript --clean-dir $resultRoot --baseline-dir $baselineRoot --output-dir $auditRoot
if ($LASTEXITCODE -ne 0) {
    throw "Clean-room comparison failed with exit code $LASTEXITCODE."
}

Write-Host "GSE232240 independent clean-room workflow complete."
Write-Host "Audit: $(Join-Path $auditRoot 'GSE232240_CLEANROOM_COMPARISON_AUDIT.md')"
