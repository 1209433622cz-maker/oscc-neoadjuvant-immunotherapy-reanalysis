# GSE301741 Response Labels Recovery Runbook

Created: 2026-07-22

## Principle

GSE301741 response labels must be recovered at patient level with provenance. Do not infer labels from plot colours or narrative snippets unless they can be cross-checked against a supplementary table or deposited metadata.

## Route 1: Publication Supplement

Open these pages in a normal browser and download the supplementary material, especially the PDF/table containing pTR or outcome annotations:

- GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741
- Article page: https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(26)00132-1
- Candidate supplement PDF: https://www.cell.com/cms/10.1016/j.xcrm.2026.102715/attachment/391a66fd-a2a3-422b-ab6f-bd6cd895187f/mmc5.pdf

Save downloaded files here:

```powershell
H:\SCI2\OSCC-GSE200996-2025.12\02_references\external_supplements\GSE301741
```

Check that the downloaded file is a real PDF rather than a blocked HTML page:

```powershell
$pdf = "H:\SCI2\OSCC-GSE200996-2025.12\02_references\external_supplements\GSE301741\mmc5.pdf"
[System.Text.Encoding]::ASCII.GetString([System.IO.File]::ReadAllBytes($pdf)[0..3])
```

Expected output:

```text
%PDF
```

If the output starts with `< !DO` or `<htm`, the request was blocked and the file is not usable.

After the patient-level labels are identified, curate them into:

```powershell
H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\config\GSE301741_RESPONSE_LABELS_MANUAL.csv
```

Allowed label fields:

- `response_label`: `responder` or `non_responder`
- `pTR_class`: `pTR-0`, `pTR-1` or `pTR-2` if available
- `pTR_percent`: numeric tumor regression percentage if available
- `source_file` and `source_detail`: exact supplement file and page/table/figure location

Then audit:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\recover_gse301741_response_labels_when_ready.ps1
```

## Route 2: High-Memory RDS Metadata

Use this only when the workstation is idle and has enough free memory. A 64 GB machine is preferred.

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\recover_gse301741_response_labels_when_ready.ps1 -ExtractRdsMetadata -MinimumFreeMemoryGB 45
```

The wrapper first extracts metadata from the deposited RDS, then runs the response-label audit.

## Current Gate

Current local status:

- Patients audited: 16
- Response labels recovered: 0
- Candidate evidence rows: 0

Until this changes, GSE301741 remains a response-pending single-cell validation candidate and must not be used for response-stratified claims.
