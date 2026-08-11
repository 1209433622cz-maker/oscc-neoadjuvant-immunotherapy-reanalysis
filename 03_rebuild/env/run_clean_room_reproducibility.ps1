param(
  [ValidateSet("prepare", "discovery", "external", "all", "compare")]
  [string]$Stage = "prepare",
  [string]$BaselineWorkspace = "H:\SCI2\OSCC-GSE200996-2025.12",
  [string]$CleanWorkspace = "H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\clean_room\clean_workspace",
  [string]$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
  [string]$Python = "D:\bioinfor\python.exe",
  [double]$MinimumFreeMemoryGB = 18,
  [switch]$Resume,
  [switch]$DeepInputHash
)

$ErrorActionPreference = "Stop"

function Get-NormalizedPath {
  param([string]$Path)
  return [System.IO.Path]::GetFullPath($Path).TrimEnd("\")
}

function Assert-Executable {
  param([string]$Path, [string]$Name)
  if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
    throw "$Name executable not found: $Path"
  }
}

function Get-FreeMemoryGB {
  $os = Get-CimInstance Win32_OperatingSystem
  return [math]::Round(($os.FreePhysicalMemory * 1KB) / 1GB, 2)
}

function Assert-FreeMemory {
  param([string]$StepName)
  $free = Get-FreeMemoryGB
  Write-Host "Free physical memory before $StepName`: $free GB"
  if ($free -lt $MinimumFreeMemoryGB) {
    throw "Free memory is below $MinimumFreeMemoryGB GB. Resume when the machine is idle."
  }
}

function Ensure-Junction {
  param([string]$LinkPath, [string]$TargetPath)
  if (Test-Path -LiteralPath $LinkPath) {
    $item = Get-Item -LiteralPath $LinkPath -Force
    if ($item.LinkType -ne "Junction") {
      throw "Expected an NTFS junction but found a regular item: $LinkPath"
    }
    $actual = Get-NormalizedPath ([string]$item.Target)
    $expected = Get-NormalizedPath $TargetPath
    if ($actual -ne $expected) {
      throw "Junction target mismatch: $LinkPath -> $actual; expected $expected"
    }
    return
  }
  $parent = Split-Path -Parent $LinkPath
  New-Item -ItemType Directory -Force -Path $parent | Out-Null
  New-Item -ItemType Junction -Path $LinkPath -Target $TargetPath | Out-Null
}

function Copy-DirectoryContents {
  param([string]$Source, [string]$Destination)
  if (-not (Test-Path -LiteralPath $Source -PathType Container)) {
    throw "Source directory not found: $Source"
  }
  New-Item -ItemType Directory -Force -Path $Destination | Out-Null
  Get-ChildItem -LiteralPath $Source -Recurse -File -Force | Where-Object {
    $_.FullName -notmatch "[\\/]__pycache__[\\/]" -and $_.Extension -ne ".pyc"
  } | ForEach-Object {
    $relative = $_.FullName.Substring($Source.TrimEnd("\").Length).TrimStart("\")
    $target = Join-Path $Destination $relative
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $target) | Out-Null
    Copy-Item -LiteralPath $_.FullName -Destination $target -Force
  }
}

$BaselineWorkspace = Get-NormalizedPath $BaselineWorkspace
$CleanWorkspace = Get-NormalizedPath $CleanWorkspace
if (-not (Test-Path -LiteralPath $BaselineWorkspace -PathType Container)) {
  throw "Baseline workspace not found: $BaselineWorkspace"
}
if ($CleanWorkspace -eq $BaselineWorkspace) {
  throw "CleanWorkspace must differ from BaselineWorkspace."
}
$allowedRoot = Get-NormalizedPath (Join-Path $BaselineWorkspace "03_rebuild\clean_room")
if (-not $CleanWorkspace.StartsWith($allowedRoot + "\", [StringComparison]::OrdinalIgnoreCase)) {
  throw "CleanWorkspace must be inside: $allowedRoot"
}

Assert-Executable -Path $Rscript -Name "Rscript"
Assert-Executable -Path $Python -Name "Python"

New-Item -ItemType Directory -Force -Path $CleanWorkspace | Out-Null
Ensure-Junction `
  -LinkPath (Join-Path $CleanWorkspace "00_raw_data") `
  -TargetPath (Join-Path $BaselineWorkspace "00_raw_data")
Ensure-Junction `
  -LinkPath (Join-Path $CleanWorkspace "02_references") `
  -TargetPath (Join-Path $BaselineWorkspace "02_references")

$cleanRebuild = Join-Path $CleanWorkspace "03_rebuild"
$analysisDir = Join-Path $cleanRebuild "analysis"
$configDir = Join-Path $cleanRebuild "config"
$logDir = Join-Path $cleanRebuild "logs"
$markerDir = Join-Path $logDir "clean_room_markers"
New-Item -ItemType Directory -Force -Path $logDir, $markerDir | Out-Null
Copy-DirectoryContents `
  -Source (Join-Path $BaselineWorkspace "03_rebuild\analysis") `
  -Destination $analysisDir
Copy-DirectoryContents `
  -Source (Join-Path $BaselineWorkspace "03_rebuild\config") `
  -Destination $configDir

$env:GSE200996_WORKSPACE = $CleanWorkspace
$env:GSE200996_BASEDIR = $cleanRebuild -replace "\\", "/"
$env:GSE200996_RAWDIR = (
  Join-Path $CleanWorkspace "00_raw_data\GSE200996_RAW"
) -replace "\\", "/"
$env:GSE200996_META = (
  Join-Path $CleanWorkspace "00_raw_data\GSE200996_metadata\GSE200996_CD45.tumor.single.cell.meta.data.txt.gz"
) -replace "\\", "/"
$env:GSE200996_USE_CACHED_OBJ = "FALSE"
$env:OMP_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$bundledPdftoppm = Join-Path $env:USERPROFILE `
  ".cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe"
if (Test-Path -LiteralPath $bundledPdftoppm -PathType Leaf) {
  $env:POPPLER_PDFTOPPM = $bundledPdftoppm
}

$provenanceScript = Join-Path $analysisDir "70_clean_room_provenance_and_compare.py"
$provenanceOutput = Join-Path $logDir "CLEAN_ROOM_PROVENANCE_MANIFEST.csv"
$manifestArgs = @(
  $provenanceScript,
  "manifest",
  "--workspace", $CleanWorkspace,
  "--output", $provenanceOutput
)
if ($DeepInputHash) {
  $manifestArgs += "--deep-input-hash"
}
& $Python @manifestArgs
if ($LASTEXITCODE -ne 0) {
  throw "Clean-room provenance manifest generation failed."
}

Write-Host "Baseline workspace: $BaselineWorkspace"
Write-Host "Clean workspace:    $CleanWorkspace"
Write-Host "Stage:              $Stage"
Write-Host "Resume:             $($Resume.IsPresent)"
Write-Host "Deep input hash:    $($DeepInputHash.IsPresent)"

if ($Stage -eq "prepare") {
  Write-Host "Clean-room preparation complete."
  exit 0
}

$runManifest = Join-Path $logDir "CLEAN_ROOM_RUN_MANIFEST.csv"

function Write-RunRecord {
  param([pscustomobject]$Record)
  if (Test-Path -LiteralPath $runManifest) {
    $Record | Export-Csv -LiteralPath $runManifest -NoTypeInformation -Append
  } else {
    $Record | Export-Csv -LiteralPath $runManifest -NoTypeInformation
  }
}

function Invoke-CleanStep {
  param(
    [string]$Name,
    [ValidateSet("R", "Python")]
    [string]$Engine,
    [string]$ScriptName,
    [string[]]$Arguments = @(),
    [switch]$MemoryGuard
  )

  $script = Join-Path $analysisDir $ScriptName
  if (-not (Test-Path -LiteralPath $script -PathType Leaf)) {
    throw "Analysis script not found: $script"
  }
  $marker = Join-Path $markerDir "$Name.success.json"
  if ($Resume -and (Test-Path -LiteralPath $marker -PathType Leaf)) {
    Write-Host "Skipping completed step: $Name"
    return
  }
  if ((-not $Resume) -and (Test-Path -LiteralPath $marker -PathType Leaf)) {
    throw "Existing success marker for $Name. Use -Resume or choose a new CleanWorkspace."
  }
  if ($MemoryGuard) {
    Assert-FreeMemory -StepName $Name
  }

  $stamp = Get-Date -Format "yyyyMMdd_HHmmss_fff"
  $stdout = Join-Path $logDir "$Name`_$stamp.out.log"
  $stderr = Join-Path $logDir "$Name`_$stamp.err.log"
  $executable = if ($Engine -eq "R") { $Rscript } else { $Python }
  $processArgs = @($script) + $Arguments
  $started = Get-Date
  Write-Host "Starting $Name"
  Write-Host "Script: $script"

  $process = Start-Process `
    -FilePath $executable `
    -ArgumentList $processArgs `
    -WorkingDirectory $CleanWorkspace `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -WindowStyle Hidden `
    -Wait `
    -PassThru

  $finished = Get-Date
  $record = [pscustomobject]@{
    step = $Name
    engine = $Engine
    script = $ScriptName
    script_sha256 = (Get-FileHash -LiteralPath $script -Algorithm SHA256).Hash
    started = $started.ToString("o")
    finished = $finished.ToString("o")
    elapsed_seconds = [math]::Round(($finished - $started).TotalSeconds, 3)
    exit_code = $process.ExitCode
    stdout = $stdout
    stderr = $stderr
  }
  Write-RunRecord -Record $record

  if ($process.ExitCode -ne 0) {
    throw "$Name failed with exit code $($process.ExitCode). Logs: $stdout ; $stderr"
  }
  $record | ConvertTo-Json | Set-Content -LiteralPath $marker -Encoding UTF8
  Write-Host "Completed $Name in $($record.elapsed_seconds) seconds"
}

function Invoke-Discovery {
  Invoke-CleanStep -Name "01_data_audit" -Engine R -ScriptName "01_data_audit.R"
  Invoke-CleanStep -Name "02_full_reanalysis" -Engine R `
    -ScriptName "run_pipeline_rebuild.R" -MemoryGuard
  Invoke-CleanStep -Name "03_leave_one_patient" -Engine R `
    -ScriptName "04_leave_one_patient_diagnostics.R" -MemoryGuard
  Invoke-CleanStep -Name "04_summary" -Engine R `
    -ScriptName "03_summarize_rebuild_results.R"
  Invoke-CleanStep -Name "05_exact_abundance" -Engine Python `
    -ScriptName "42_exact_permutation_abundance_sensitivity.py"
  Invoke-CleanStep -Name "06_cohort_adjusted_pseudobulk" -Engine R `
    -ScriptName "43_discovery_cohort_adjusted_pseudobulk_sensitivity.R" `
    -MemoryGuard
  Invoke-CleanStep -Name "07_core_figures" -Engine R `
    -ScriptName "05_make_core_figures.R"
  Invoke-CleanStep -Name "08_extended_data_figures" -Engine R `
    -ScriptName "06_make_extended_data_figures.R"
}

function Invoke-External {
  $dynamicInput = Join-Path $cleanRebuild `
    "results\dynamic_paired\Fig4B_T_cell_interaction_DE_trend.csv"
  if (-not (Test-Path -LiteralPath $dynamicInput -PathType Leaf)) {
    throw "External validation requires clean-room discovery outputs. Run -Stage discovery or -Stage all first."
  }

  Invoke-CleanStep -Name "09_metadata_manifests" -Engine Python `
    -ScriptName "24_build_metadata_manifests.py"
  Invoke-CleanStep -Name "10_response_freeze" -Engine Python `
    -ScriptName "25_freeze_response_metadata.py"
  Invoke-CleanStep -Name "11_gse123813_modules" -Engine Python `
    -ScriptName "17_external_validation_gse123813.py"
  Invoke-CleanStep -Name "12_gse123813_tcr" -Engine Python `
    -ScriptName "20_external_tcr_orthogonal_validation_gse123813.py"
  Invoke-CleanStep -Name "13_gse281729_modules" -Engine Python `
    -ScriptName "30_validate_modules_gse281729_bulk.py"
  Invoke-CleanStep -Name "14_response_freeze" -Engine Python `
    -ScriptName "31_freeze_response_metadata.py"
  Invoke-CleanStep -Name "15_gse281729_robustness" -Engine Python `
    -ScriptName "41_gse281729_external_robustness.py"
  Invoke-CleanStep -Name "16_gse281729_adaptive_timing" -Engine Python `
    -ScriptName "54_gse281729_response_adaptive_timing_sensitivity.py"
  Invoke-CleanStep -Name "17_extended_data7" -Engine Python `
    -ScriptName "56_make_extended_data7_response_adaptive_robustness.py"
  Invoke-CleanStep -Name "18_gse179730_direction" -Engine Python `
    -ScriptName "33_score_gse179730_treatment_direction.py"
  Invoke-CleanStep -Name "19_gse179730_response" -Engine Python `
    -ScriptName "55_recover_and_validate_gse179730_response.py"
  Invoke-CleanStep -Name "20_gse301741_raw_manifest" -Engine R `
    -ScriptName "27_prepare_gse301741_raw_route.R"
  Invoke-CleanStep -Name "21_gse301741_strata" -Engine Python `
    -ScriptName "29_build_gse301741_validation_strata.py"
  Invoke-CleanStep -Name "22_gse301741_figure_outcomes" -Engine Python `
    -ScriptName "38_extract_gse301741_supplement_figure_outcomes.py"
  Invoke-CleanStep -Name "23_gse301741_response_labels" -Engine Python `
    -ScriptName "37_recover_gse301741_response_labels.py"
  Invoke-CleanStep -Name "24_gse301741_tierA_pairs" -Engine Python `
    -ScriptName "69_select_gse301741_tierA_pairs.py"
  Invoke-CleanStep -Name "25_gse301741_raw_metadata" -Engine R `
    -ScriptName "44_rebuild_gse301741_raw_cell_metadata.R" -MemoryGuard
  Invoke-CleanStep -Name "26_gse301741_lineage_validation" -Engine R `
    -ScriptName "45_validate_gse301741_lineage_aware.R" -MemoryGuard
  Invoke-CleanStep -Name "27_gse301741_reconstruction_audit" -Engine R `
    -ScriptName "46_audit_gse301741_raw_reconstruction.R"
  Invoke-CleanStep -Name "28_gse301741_boundary_figure" -Engine R `
    -ScriptName "47_make_gse301741_boundary_extended_data.R"
  Invoke-CleanStep -Name "29_gse281729_validation_figure" -Engine Python `
    -ScriptName "32_make_gse281729_validation_figure.py"
  Invoke-CleanStep -Name "30_integrated_validation_figure" -Engine Python `
    -ScriptName "34_make_integrated_validation_figure.py"
  Invoke-CleanStep -Name "31_locked_family_tests" -Engine Python `
    -ScriptName "61_locked_family_robustness_external_cohorts.py"
  Invoke-CleanStep -Name "32_locked_family_figure" -Engine Python `
    -ScriptName "62_make_locked_family_robustness_figure.py"
}

function Invoke-Comparison {
  $comparisonCsv = Join-Path $logDir "CLEAN_ROOM_RESULT_COMPARISON.csv"
  $comparisonMd = Join-Path $logDir "CLEAN_ROOM_RESULT_COMPARISON.md"
  & $Python $provenanceScript compare `
    --baseline $BaselineWorkspace `
    --clean $CleanWorkspace `
    --output-csv $comparisonCsv `
    --output-md $comparisonMd
  if ($LASTEXITCODE -ne 0) {
    throw "Clean-room comparison found non-PASS tables. Review: $comparisonMd"
  }
}

switch ($Stage) {
  "discovery" { Invoke-Discovery }
  "external" { Invoke-External }
  "all" {
    Invoke-Discovery
    Invoke-External
    Invoke-Comparison
  }
  "compare" { Invoke-Comparison }
}

Write-Host "Clean-room stage complete: $Stage"
