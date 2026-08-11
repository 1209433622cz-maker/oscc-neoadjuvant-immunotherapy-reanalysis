# Data availability

No primary patient data are distributed in this repository. Reproduction requires downloading the public source files from their authoritative deposits.

| Accession | Role | Primary source |
|---|---|---|
| GSE200996 | OSCC single-cell discovery cohort | NCBI GEO |
| GSE281729 | Paired HNSCC bulk response-association boundary | NCBI GEO |
| GSE195832 | Independently eligible paired HNSCC bulk boundary | NCBI GEO and the author-linked Mendeley deposit |
| GSE179730 | Paired OCSCC exact-test boundary | NCBI GEO |
| GSE301741 | Disease-matched single-cell boundary | NCBI GEO |
| GSE123813 | Cross-disease anti-PD-1 expression and TCR context | NCBI GEO |
| GSE232240 | Prospective frozen independent longitudinal immune validation | NCBI GEO |

The original registered file URLs and priorities are provided in `03_rebuild/config/external_download_manifest.csv`. The GSE232240 archive, count matrix and cell metadata are registered separately in `03_rebuild/config/external_download_manifest.csv`. Download helpers write files below `00_raw_data/external_validation/`, which is excluded from Git.

Response endpoints were retained as defined by each source. They were not silently harmonized across cohorts. GSE301741 response labels recovered from publication-supplement graphics remain provisional and are used only for boundary analysis.

The deposited 12.35-GB GSE301741 Seurat RDS is not required for the memory-safe route. The registered workflow instead reconstructs eligible broad lineages from Cell Ranger-filtered H5 matrices one sample at a time.

GSE232240 is acquired through its 61,460,480-byte GEO raw archive. The registered script uses IPv4, retries and resume support, verifies the tar members and checks the extracted gzip streams. The deposited cell metadata have 32,399 unique cells from 18 patients; no expression-association result is stored in this repository at registration time.
