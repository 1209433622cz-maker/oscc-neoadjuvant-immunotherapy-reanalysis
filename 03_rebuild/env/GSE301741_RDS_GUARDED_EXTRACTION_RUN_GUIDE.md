# GSE301741 RDS Guarded Metadata Extraction

Date: 2026-07-26

## Current Verdict

The verified 12.354 GB RDS begins with `58-0A` (`X\n`), identifying an uncompressed XDR R serialization stream written by R 4.3.1. Both R 4.6.0 and project-local R 4.3.3 fail with `ReadItem: unknown type 66` at the same byte offset before meaningful memory allocation. A remote GEO byte range spanning the failure offset is identical to the local file.

This is therefore not a 31.8 GB workstation-memory failure. The deposited stream is malformed or corrupted. Do not lower the memory guard and do not move the unchanged file to a larger machine expecting RAM alone to fix it.

The guarded route runs R in an isolated below-normal-priority process and monitors:

- Free physical memory.
- Free virtual/commit memory.
- R working-set memory.
- R private memory.
- Elapsed time.

The process is terminated before system commit exhaustion. Partial output files use a `.partial` suffix and are not promoted to final outputs.

## Replacement Local Command

Use the sample-wise RAW H5 reconstruction route:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_gse301741_raw_reconstruction.ps1
```

This route requires 8 GB free physical memory by default and never holds all 137,020 cells in one R object.

## Expected Outputs

- `03_rebuild/manifests/GSE301741_RDS_CELL_METADATA.csv.gz`
- `03_rebuild/manifests/GSE301741_RDS_CELL_METADATA.rds`
- `03_rebuild/manifests/GSE301741_RDS_METADATA_FIELD_SUMMARY.csv`
- `03_rebuild/manifests/GSE301741_RDS_METADATA_VALUE_COUNTS.csv`
- `03_rebuild/manifests/GSE301741_RDS_COMPACT_CLINICAL_CELLTYPE_FIELDS.csv`
- `03_rebuild/manifests/GSE301741_RDS_OBJECT_STRUCTURE.csv`
- `03_rebuild/manifests/GSE301741_RDS_METADATA_EXTRACTION_REPORT.md`
- `03_rebuild/manifests/GSE301741_RDS_METADATA_EXTRACTION_SUCCESS.txt`

## Corrected-File Handoff

If the authors or GEO provide a corrected RDS or a cell-level metadata table, preserve it under a new file name, record its checksum, and then use the guarded extractor. The guarded extractor now refuses the known-broken object unless `-RetryKnownBrokenRds` is passed explicitly.

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\extract_gse301741_metadata_commit_guarded.ps1 -RetryKnownBrokenRds
```

Only use that switch for a deliberate diagnostic retry or after updating the configured RDS path to a corrected source.
