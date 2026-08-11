# User-side download commands for the OSCC/HNSCC manuscript rebuild.
# Run selected blocks in PowerShell from any directory.

# 1. Resume or download only the failed 12.4 GB GSE301741 Seurat RDS.
# This uses curl.exe with resume support and continues from the existing .partial file.
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -OnlyAccession GSE301741 -OnlyFileName GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds -Retries 20

# 2. Verify only the large RDS.
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\verify_external_downloads.ps1 -OnlyAccession GSE301741 -OnlyFileName GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds

# 3. Verify all default P0/P1 acquisition-gate files.
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\verify_external_downloads.ps1

# 4. Optional P2 datasets for spatial localization and secondary paired sensitivity.
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -IncludeP2 -Retries 20

# 5. Optional P2 plus P3 datasets, including baseline-only targeted expression sensitivity.
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\download_external_datasets.ps1 -IncludeP2 -IncludeP3 -Retries 20

# 6. Verify optional P2/P3 files after downloading them.
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\verify_external_downloads.ps1 -IncludeP2 -IncludeP3
