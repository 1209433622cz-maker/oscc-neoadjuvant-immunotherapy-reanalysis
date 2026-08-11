# External Dataset Download Handoff

Date: 2026-07-22

This handoff prepares data acquisition for the manuscript rebuild. It does not run analysis.

## Default download

The default command downloads:

- GSE200996 metadata refresh files.
- GSE301741 raw archive, large Seurat object, GEO metadata, and SRA RunInfo.
- GSE281729 processed bulk RNA-seq matrix, GEO metadata, and SRA RunInfo.
- GSE179730 processed bulk RNA-seq matrix, GEO metadata, and SRA RunInfo.

Run:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1
```

## Metadata-only check

Use this first if you want to test the network path without downloading the 12.4 GB GSE301741 RDS.

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -MetadataOnly
```

## Optional spatial and secondary cohorts

Download P2 datasets after the P1 signal review or when storage is available:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -IncludeP2
```

Download P2 plus the P3 targeted-expression baseline sensitivity cohort:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -IncludeP2 -IncludeP3
```

## Verification

After download, run:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\verify_external_downloads.ps1
```

## Resume the large GSE301741 RDS only

If the 12.4 GB RDS fails with an unexpected EOF, rerun only that file. The downloader now uses `curl.exe` with resume support and will continue from the `.partial` file when present.

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -OnlyAccession GSE301741 -OnlyFileName GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds -Retries 20
```

Then verify only that file:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\verify_external_downloads.ps1 -OnlyAccession GSE301741 -OnlyFileName GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds
```

For optional cohorts:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\verify_external_downloads.ps1 -IncludeP2 -IncludeP3
```

## Output locations

- Downloaded files: `H:\SCI2\OSCC-GSE200996-2025.12\00_raw_data\external_validation`
- Download manifest: `H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\config\external_download_manifest.csv`
- Download status: `H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\logs\external_download_status_latest.csv`
- Verification status: `H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\logs\external_download_verification_latest.csv`

## Notes

- The script writes SHA256 checksums for every completed file.
- Existing files are skipped unless `-Force` is provided.
- SRA entries are RunInfo metadata only. The script does not download FASTQ or SRA files.
- `GSE200996_PRJNA827834_SRA_RunInfo.csv` is retained in the manifest but is not a default acquisition gate because the endpoint currently returns zero bytes.
- GSE301741 has a 12.4 GB RDS file and may require roughly 40-60 GB RAM after loading in R.
