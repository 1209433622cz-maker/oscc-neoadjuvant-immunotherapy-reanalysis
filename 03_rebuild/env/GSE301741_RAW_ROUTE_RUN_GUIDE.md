# GSE301741 Raw-Route Run Guide

Date: 2026-07-22

Use this route when the full 12.4 GB Seurat RDS cannot be loaded safely.

## Why

The local machine has about 34 GB RAM. The GSE301741 RDS is 12.4 GB on disk and may need more memory after loading. The raw archive is about 586.5 MB and contains 58 per-sample files, so it is safer to process sample-by-sample.

## First check only

This writes a manifest report without extracting the tar archive.

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\prepare_gse301741_raw_route_when_idle.ps1
```

## Extract and index

Run when current local analyses are finished. This extracts `GSE301741_RAW.tar` and creates per-sample H5/TCR QC summaries.

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\prepare_gse301741_raw_route_when_idle.ps1 -Extract
```

## Outputs

- `H:\SCI2\OSCC-GSE200996-2025.12\00_raw_data\external_validation\GSE301741\RAW_extracted`
- `H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\manifests\GSE301741_RAW_ROUTE_SAMPLE_QC.csv`
- `H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\manifests\GSE301741_RAW_ROUTE_SUMMARY.csv`
- `H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\manifests\GSE301741_RAW_ROUTE_PREP_REPORT.md`

## Decision

Use the raw route for validation unless a high-memory machine is available for extracting response metadata from the full RDS.
