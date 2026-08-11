from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import tarfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
RAW_DIR = WORKSPACE / "00_raw_data" / "GSE200996_RAW"
EXT_DIR = WORKSPACE / "00_raw_data" / "external_validation"
OUT_DIR = WORKSPACE / "03_rebuild" / "manifests"
LOG_DIR = WORKSPACE / "03_rebuild" / "logs"
TABLE_DIR = WORKSPACE / "03_rebuild" / "tables" / "submission" / "csv"

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def parse_series_matrix(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    series: dict[str, str] = {}
    sample_rows: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith("!"):
                continue
            parts = line.split("\t")
            key = parts[0].lstrip("!")
            values = [strip_quotes(x) for x in parts[1:]]
            if key.startswith("Series_"):
                series.setdefault(key, " | ".join(values))
            elif key == "Sample_characteristics_ch1":
                char_rows.append(values)
            elif key.startswith("Sample_"):
                sample_rows[key] = values

    accessions = sample_rows.get("Sample_geo_accession", [])
    records: list[dict[str, str]] = []
    for idx, gsm in enumerate(accessions):
        rec: dict[str, str] = {
            "geo_accession": gsm,
            "source_matrix": str(path),
            "series_accession": series.get("Series_geo_accession", ""),
            "series_title": series.get("Series_title", ""),
            "series_pubmed_id": series.get("Series_pubmed_id", ""),
        }
        for key, values in sample_rows.items():
            if idx < len(values):
                rec[key.replace("Sample_", "").lower()] = values[idx]
        for values in char_rows:
            if idx >= len(values):
                continue
            raw = values[idx]
            if ":" in raw:
                field, value = raw.split(":", 1)
                field = re.sub(r"[^A-Za-z0-9]+", "_", field.strip().lower()).strip("_")
                rec[f"characteristic_{field}"] = value.strip()
            else:
                rec.setdefault("characteristic_unparsed", "")
                rec["characteristic_unparsed"] += (" | " if rec["characteristic_unparsed"] else "") + raw
        records.append(rec)
    return series, records


def parse_patient_tokens(text: str) -> list[str]:
    tokens = re.findall(r"P\d{2}", text)
    return sorted(set(tokens))


def harmonize_timepoint(value: str) -> str:
    v = value.lower()
    if re.search(r"\bpre\b|pre-tx|baseline|first of", v):
        return "pre"
    if re.search(r"\bpost\b|post-tx|surgery|second of", v):
        return "post"
    if re.search(r"\bb1\b", v):
        return "B1"
    if re.search(r"\bb2\b", v):
        return "B2"
    if re.search(r"\bb3\b", v):
        return "B3"
    if "recur" in v:
        return "recurrence"
    return ""


def parse_gse200996_raw_file(path: Path) -> list[dict[str, str]]:
    name = path.name
    gsm = re.match(r"^(GSM\d+)", name)
    gsm_id = gsm.group(1) if gsm else ""
    patients = parse_patient_tokens(name)
    if not patients:
        patients = [""]

    lower = name.lower()
    compartment = "PBMC" if "pbmc" in lower else "tumor" if "tumor" in lower else ""
    modality = ""
    if "gex_sc" in lower or "feature_bc_matrix" in lower:
        modality = "scRNA"
    if "tcr_sc" in lower or "contig_annotations" in lower:
        modality = "scTCR"
    if "tcr_bulk" in lower:
        modality = "bulk_TCR"

    timepoint_raw = ""
    if "pre-tx" in lower:
        timepoint_raw = "pre-Tx"
    elif "post-tx" in lower:
        timepoint_raw = "post-Tx"
    else:
        b = re.search(r"_(B[123])_", name)
        if b:
            timepoint_raw = b.group(1)

    sort_type = ""
    for token in ["CD4", "CD8", "CD45", "CD45pos", "CD45neg", "CD45ratio", "CD3"]:
        if token.lower() in lower:
            sort_type = token
            break

    rows = []
    for patient in patients:
        rows.append(
            {
                "accession": "GSE200996",
                "geo_accession": gsm_id,
                "patient_id": patient,
                "specimen_id": f"{patient}_{timepoint_raw}_{compartment}_{modality}".strip("_"),
                "compartment": compartment,
                "timepoint_raw": timepoint_raw,
                "timepoint_harmonized": harmonize_timepoint(timepoint_raw),
                "modality": modality,
                "assay_or_library": "",
                "sort_or_cell_fraction": sort_type,
                "source_file": str(path),
                "file_name": name,
                "size_bytes": path.stat().st_size,
                "inclusion_flag": "include_pending_qc",
                "exclusion_reason": "",
                "notes": "expanded from local GSE200996 raw file name",
            }
        )
    return rows


def parse_gse301741_tar_listing(tar_path: Path) -> list[dict[str, str]]:
    rows = []
    if not tar_path.exists():
        return rows
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            name = Path(member.name).name
            m = re.match(
                r"^(GSM\d+)_(?P<patient>[A-Z]+SCC\d+)-(?P<timepoint>pre|post)(?:-(?P<sort>[^_]+))?-(?P<library>3|5|TCR)",
                name,
                flags=re.IGNORECASE,
            )
            if not m:
                continue
            lib = m.group("library")
            modality = "scTCR" if lib.upper() == "TCR" else "scRNA"
            patient = m.group("patient")
            rows.append(
                {
                    "accession": "GSE301741",
                    "geo_accession": m.group(1),
                    "patient_id": patient,
                    "specimen_id": f"{patient}_{m.group('timepoint')}_{m.group('sort') or 'unsorted'}_{lib}",
                    "compartment": "tumor",
                    "timepoint_raw": m.group("timepoint"),
                    "timepoint_harmonized": harmonize_timepoint(m.group("timepoint")),
                    "modality": modality,
                    "assay_or_library": f"{lib}prime" if lib in {"3", "5"} else "TCR",
                    "sort_or_cell_fraction": m.group("sort") or "Unsorted",
                    "source_file": str(tar_path),
                    "file_name": name,
                    "size_bytes": member.size,
                    "inclusion_flag": "include_pending_qc",
                    "exclusion_reason": "",
                    "notes": "expanded from GSE301741_RAW.tar listing",
                }
            )
    return rows


def parse_external_sample_title(accession: str, rec: dict[str, str]) -> dict[str, str]:
    title = rec.get("title", "")
    patient = ""
    specimen = title
    timepoint = rec.get("characteristic_timepoint", "")
    treatment = rec.get("characteristic_treatment", "")
    response = ""
    site = rec.get("characteristic_tissue", "") or rec.get("source_name_ch1", "")
    modality = "bulk_RNA" if accession in {"GSE281729", "GSE179730"} else ""

    if accession == "GSE301741":
        m = re.match(r"([A-Z]+SCC\d+),\s*(Pre|Post),\s*([^,]+),\s*([^,]+)", title, flags=re.I)
        if m:
            patient, timepoint, sort_type, library = m.groups()
            modality = "scTCR" if library.upper() == "TCR" else "scRNA"
            return {
                "patient_id": patient,
                "specimen_id": f"{patient}_{timepoint}_{sort_type}_{library}",
                "timepoint_raw": timepoint,
                "timepoint_harmonized": harmonize_timepoint(timepoint),
                "modality": modality,
                "assay_or_library": library,
                "sort_or_cell_fraction": sort_type,
                "site": site,
                "treatment": "neoadjuvant pembrolizumab",
                "response_raw": "",
            }

    if accession == "GSE281729":
        m = re.match(r"(NI_TJ3_\d+)_(\d+)", title)
        if m:
            patient = m.group(1)
        if "Baseline" in title:
            timepoint = "Baseline"
        elif "Surgery" in title:
            timepoint = "Surgery"
        r = re.search(r"Primary Path Response_([^\]]+)", title)
        if r:
            response = r.group(1).replace("Respoder", "Responder")
        return {
            "patient_id": patient,
            "specimen_id": specimen,
            "timepoint_raw": timepoint,
            "timepoint_harmonized": harmonize_timepoint(timepoint),
            "modality": modality,
            "assay_or_library": "bulk RNA-seq",
            "sort_or_cell_fraction": "bulk tumor",
            "site": site,
            "treatment": treatment,
            "response_raw": response,
        }

    if accession == "GSE179730":
        m = re.match(r"(MUSC-HN\d+)-(Pre|Post|Recur)-Tumor", title, flags=re.I)
        if m:
            patient = m.group(1)
            timepoint = m.group(2)
        return {
            "patient_id": patient,
            "specimen_id": specimen,
            "timepoint_raw": timepoint,
            "timepoint_harmonized": harmonize_timepoint(timepoint),
            "modality": modality,
            "assay_or_library": "bulk RNA-seq",
            "sort_or_cell_fraction": "bulk tumor",
            "site": site,
            "treatment": "neoadjuvant anti-PD-1",
            "response_raw": "",
        }

    return {
        "patient_id": patient,
        "specimen_id": specimen,
        "timepoint_raw": timepoint,
        "timepoint_harmonized": harmonize_timepoint(timepoint),
        "modality": modality,
        "assay_or_library": "",
        "sort_or_cell_fraction": "",
        "site": site,
        "treatment": treatment,
        "response_raw": response,
    }


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        keys = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_acquisition_inventory() -> list[dict[str, object]]:
    rows = []
    manifest = read_csv(WORKSPACE / "03_rebuild" / "config" / "external_download_manifest.csv")
    for item in manifest:
        target = EXT_DIR / item["target_subdir"] / item["file_name"]
        exists = target.exists()
        rows.append(
            {
                **item,
                "exists": exists,
                "zero_byte": exists and target.stat().st_size == 0,
                "size_bytes": target.stat().st_size if exists else "",
                "sha256": sha256_file(target) if exists and target.stat().st_size > 0 else "",
                "local_path": str(target),
            }
        )
    return rows


def build_gse200996_manifest() -> list[dict[str, object]]:
    response_path = TABLE_DIR / "S1_patient_response.csv"
    if not response_path.exists():
        response_path = (
            WORKSPACE
            / "03_rebuild"
            / "results"
            / "data_audit"
            / "cd45_tumor_patient_response_table.csv"
        )
    response_rows = read_csv(response_path)
    response_by_patient = {r["patient_id"]: r for r in response_rows}
    rows: list[dict[str, object]] = []
    for path in sorted(RAW_DIR.glob("*")):
        if path.is_file():
            for row in parse_gse200996_raw_file(path):
                resp = response_by_patient.get(row["patient_id"], {})
                row.update(
                    {
                        "treatment_arm": resp.get("cohort", ""),
                        "path_response_category": resp.get("path_response", ""),
                        "path_response_binary": resp.get("response_bin", ""),
                        "path_response_ord_num": resp.get("response_ord_num", ""),
                        "path_response_percent": "",
                        "response_source": "S1_patient_response.csv; continuous percent not yet recovered",
                    }
                )
                rows.append(row)
    return rows


def build_external_manifest() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sample_meta_rows: list[dict[str, object]] = []
    harmonized_rows: list[dict[str, object]] = []
    matrix_files = [
        EXT_DIR / "GSE301741" / "GSE301741_series_matrix.txt.gz",
        EXT_DIR / "GSE281729" / "GSE281729_series_matrix.txt.gz",
        EXT_DIR / "GSE179730" / "GSE179730_series_matrix.txt.gz",
    ]
    for matrix in matrix_files:
        if not matrix.exists():
            continue
        series, records = parse_series_matrix(matrix)
        accession = series.get("Series_geo_accession", "")
        for rec in records:
            sample_meta_rows.append(rec)
            parsed = parse_external_sample_title(accession, rec)
            harmonized_rows.append(
                {
                    "accession": accession,
                    "geo_accession": rec.get("geo_accession", ""),
                    "patient_id": parsed["patient_id"],
                    "specimen_id": parsed["specimen_id"],
                    "compartment": "tumor",
                    "timepoint_raw": parsed["timepoint_raw"],
                    "timepoint_harmonized": parsed["timepoint_harmonized"],
                    "modality": parsed["modality"],
                    "assay_or_library": parsed["assay_or_library"],
                    "sort_or_cell_fraction": parsed["sort_or_cell_fraction"],
                    "site": parsed["site"],
                    "treatment": parsed["treatment"],
                    "response_raw": parsed["response_raw"],
                    "response_harmonized": "",
                    "response_source": "GEO series matrix; continuous response not yet recovered",
                    "source_file": str(matrix),
                    "inclusion_flag": "include_pending_qc",
                    "exclusion_reason": "",
                    "notes": "",
                }
            )

    tar_rows = parse_gse301741_tar_listing(EXT_DIR / "GSE301741" / "GSE301741_RAW.tar")
    tar_by_gsm = {row["geo_accession"]: row for row in tar_rows}
    matched = set()
    for row in harmonized_rows:
        if row.get("accession") != "GSE301741":
            continue
        tar_row = tar_by_gsm.get(str(row.get("geo_accession", "")))
        if not tar_row:
            continue
        matched.add(str(row.get("geo_accession", "")))
        row["raw_archive_file_name"] = tar_row["file_name"]
        row["raw_archive_member_size_bytes"] = tar_row["size_bytes"]
        row["source_file"] = tar_row["source_file"]
        row["notes"] = "GEO sample metadata merged with GSE301741_RAW.tar listing"

    for gsm, tar_row in sorted(tar_by_gsm.items()):
        if gsm in matched:
            continue
        tar_row = dict(tar_row)
        tar_row["raw_archive_file_name"] = tar_row["file_name"]
        tar_row["raw_archive_member_size_bytes"] = tar_row["size_bytes"]
        tar_row["site"] = ""
        tar_row["treatment"] = "neoadjuvant pembrolizumab"
        tar_row["response_raw"] = ""
        tar_row["response_harmonized"] = ""
        tar_row["response_source"] = "GSE301741_RAW.tar listing only; GEO metadata row not matched"
        harmonized_rows.append(tar_row)
    return sample_meta_rows, harmonized_rows


def summarize(rows: list[dict[str, object]], fields: list[str]) -> list[dict[str, object]]:
    out = []
    for field in fields:
        counts = Counter(str(row.get(field, "")) for row in rows)
        for value, count in sorted(counts.items()):
            out.append({"field": field, "value": value, "count": count})
    return out


def main() -> int:
    acquisition = build_acquisition_inventory()
    gse200996 = build_gse200996_manifest()
    external_meta, external = build_external_manifest()

    write_csv(OUT_DIR / "ACQUISITION_FILE_INVENTORY.csv", acquisition)
    write_csv(OUT_DIR / "GSE200996_LOCAL_MULTIMODAL_MANIFEST.csv", gse200996)
    write_csv(OUT_DIR / "EXTERNAL_GEO_SAMPLE_METADATA_LONG.csv", external_meta)
    write_csv(OUT_DIR / "EXTERNAL_COHORT_SAMPLE_MANIFEST.csv", external)

    summary_rows = []
    summary_rows.extend({"manifest": "GSE200996", **x} for x in summarize(gse200996, ["modality", "compartment", "timepoint_harmonized", "path_response_category", "treatment_arm"]))
    summary_rows.extend({"manifest": "external", **x} for x in summarize(external, ["accession", "modality", "timepoint_harmonized", "response_raw", "site"]))
    write_csv(OUT_DIR / "METADATA_MANIFEST_SUMMARY.csv", summary_rows)

    patient_summary: list[dict[str, object]] = []
    for accession, rows in [
        ("GSE200996", gse200996),
        ("external_all", external),
    ]:
        by_patient: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            patient = str(row.get("patient_id", ""))
            if patient:
                by_patient[patient].append(row)
        for patient, patient_rows in sorted(by_patient.items()):
            patient_summary.append(
                {
                    "manifest": accession,
                    "patient_id": patient,
                    "n_records": len(patient_rows),
                    "accessions": ";".join(sorted(set(str(x.get("accession", "")) for x in patient_rows))),
                    "modalities": ";".join(sorted(set(str(x.get("modality", "")) for x in patient_rows if x.get("modality")))),
                    "timepoints": ";".join(sorted(set(str(x.get("timepoint_harmonized", "")) for x in patient_rows if x.get("timepoint_harmonized")))),
                }
            )
    write_csv(OUT_DIR / "PATIENT_MODALITY_COVERAGE.csv", patient_summary)

    report = {
        "created": datetime.now().astimezone().isoformat(timespec="seconds"),
        "workspace": str(WORKSPACE),
        "acquisition_files": len(acquisition),
        "acquisition_present": sum(1 for r in acquisition if r["exists"]),
        "acquisition_zero_byte": sum(1 for r in acquisition if r["zero_byte"]),
        "gse200996_manifest_rows": len(gse200996),
        "gse200996_patients": len(set(r["patient_id"] for r in gse200996 if r["patient_id"])),
        "external_geo_metadata_rows": len(external_meta),
        "external_manifest_rows": len(external),
        "external_patients": {
            acc: len(set(r["patient_id"] for r in external if r.get("accession") == acc and r.get("patient_id")))
            for acc in sorted(set(str(r.get("accession", "")) for r in external))
        },
        "limitations": [
            "Continuous pathological response values are not yet recovered.",
            "GSE301741 RDS was not loaded because local RAM is below the expected 40-60 GB working memory range.",
            "GSE301741 raw tar listing and GEO matrix metadata are sufficient for sample-level acquisition manifesting.",
            "External response harmonization remains locked as pending until authoritative clinical fields are parsed.",
        ],
    }
    (OUT_DIR / "METADATA_MANIFEST_BUILD_REPORT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    md = [
        "# Metadata Manifest Build Report",
        "",
        f"Created: {report['created']}",
        "",
        f"- Acquisition manifest entries: {report['acquisition_files']}",
        f"- Acquisition files present: {report['acquisition_present']}",
        f"- Acquisition zero-byte files: {report['acquisition_zero_byte']}",
        f"- GSE200996 multimodal rows: {report['gse200996_manifest_rows']}",
        f"- GSE200996 patients represented: {report['gse200996_patients']}",
        f"- External GEO metadata rows: {report['external_geo_metadata_rows']}",
        f"- External harmonized rows: {report['external_manifest_rows']}",
        "",
        "## External patient counts",
        "",
        "| Accession | Patients |",
        "|---|---:|",
    ]
    for acc, count in report["external_patients"].items():
        md.append(f"| {acc} | {count} |")
    md.extend(
        [
            "",
            "## Limitations",
            "",
        ]
    )
    for item in report["limitations"]:
        md.append(f"- {item}")
    md.extend(
        [
            "",
            "## Files written",
            "",
            "- `03_rebuild/manifests/ACQUISITION_FILE_INVENTORY.csv`",
            "- `03_rebuild/manifests/GSE200996_LOCAL_MULTIMODAL_MANIFEST.csv`",
            "- `03_rebuild/manifests/EXTERNAL_GEO_SAMPLE_METADATA_LONG.csv`",
            "- `03_rebuild/manifests/EXTERNAL_COHORT_SAMPLE_MANIFEST.csv`",
            "- `03_rebuild/manifests/PATIENT_MODALITY_COVERAGE.csv`",
            "- `03_rebuild/manifests/METADATA_MANIFEST_SUMMARY.csv`",
            "- `03_rebuild/manifests/METADATA_MANIFEST_BUILD_REPORT.json`",
        ]
    )
    (OUT_DIR / "METADATA_MANIFEST_BUILD_REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
