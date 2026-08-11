from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
MANIFEST_DIR = WORKSPACE / "03_rebuild" / "manifests"
OUT_DIR = WORKSPACE / "03_rebuild" / "validation" / "GSE301741_stratification"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def join_unique(values: list[str]) -> str:
    return ";".join(sorted({v for v in values if v and v != "NA"}))


def combos(rows: list[dict[str, str]]) -> set[tuple[str, str]]:
    return {
        (row.get("sort_or_cell_fraction", ""), row.get("assay_or_library", ""))
        for row in rows
        if row.get("sort_or_cell_fraction") and row.get("assay_or_library")
    }


def fractions(rows: list[dict[str, str]]) -> set[str]:
    return {row.get("sort_or_cell_fraction", "") for row in rows if row.get("sort_or_cell_fraction")}


def validation_tier(pre_rows: list[dict[str, str]], post_rows: list[dict[str, str]]) -> tuple[str, str]:
    if not pre_rows or not post_rows:
        return "D_post_or_pre_only", "Not usable for paired response validation."
    shared_combo = combos(pre_rows) & combos(post_rows)
    if shared_combo:
        detail = ", ".join(f"{fraction}/{library}" for fraction, library in sorted(shared_combo))
        return "A_same_fraction_same_library", f"Best paired scRNA validation subset available: {detail}."
    shared_fraction = fractions(pre_rows) & fractions(post_rows)
    if shared_fraction:
        detail = ", ".join(sorted(shared_fraction))
        return "B_same_fraction_mixed_library", f"Paired fraction available but library chemistry differs: {detail}."
    return "C_prepost_mixed_fraction", "Pre/post available, but sorting fractions differ; avoid abundance claims."


def main() -> int:
    qc_path = MANIFEST_DIR / "GSE301741_RAW_ROUTE_SAMPLE_QC.csv"
    response_path = MANIFEST_DIR / "RESPONSE_METADATA_FREEZE.csv"
    if not qc_path.exists():
        raise FileNotFoundError(qc_path)

    rows = read_csv(qc_path)
    response_rows = read_csv(response_path) if response_path.exists() else []
    response_by_patient = {
        row.get("patient_id", ""): row
        for row in response_rows
        if row.get("accession") == "GSE301741" and row.get("patient_id")
    }

    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        patient = row.get("patient_id", "")
        if patient:
            by_patient[patient].append(row)

    out: list[dict[str, object]] = []
    for patient, patient_rows in sorted(by_patient.items()):
        scrna = [
            row for row in patient_rows
            if row.get("modality") == "scRNA" and row.get("h5_status") == "ok"
        ]
        tcr = [
            row for row in patient_rows
            if row.get("modality") == "scTCR" and row.get("tcr_status") == "ok"
        ]
        pre_scrna = [row for row in scrna if row.get("timepoint_harmonized") == "pre"]
        post_scrna = [row for row in scrna if row.get("timepoint_harmonized") == "post"]
        pre_tcr = [row for row in tcr if row.get("timepoint_harmonized") == "pre"]
        post_tcr = [row for row in tcr if row.get("timepoint_harmonized") == "post"]

        tier, recommendation = validation_tier(pre_scrna, post_scrna)
        shared_scrna_combos = combos(pre_scrna) & combos(post_scrna)
        shared_tcr_fractions = fractions(pre_tcr) & fractions(post_tcr)
        response = response_by_patient.get(patient, {})

        out.append(
            {
                "patient_id": patient,
                "site": join_unique([row.get("site", "") for row in patient_rows]),
                "n_scrna_h5_ok": len(scrna),
                "n_tcr_ok": len(tcr),
                "has_pre_scrna": bool(pre_scrna),
                "has_post_scrna": bool(post_scrna),
                "has_prepost_scrna": bool(pre_scrna and post_scrna),
                "pre_scrna_fractions": join_unique([row.get("sort_or_cell_fraction", "") for row in pre_scrna]),
                "post_scrna_fractions": join_unique([row.get("sort_or_cell_fraction", "") for row in post_scrna]),
                "pre_scrna_libraries": join_unique([row.get("assay_or_library", "") for row in pre_scrna]),
                "post_scrna_libraries": join_unique([row.get("assay_or_library", "") for row in post_scrna]),
                "shared_scrna_fraction_library": ";".join(
                    f"{fraction}/{library}" for fraction, library in sorted(shared_scrna_combos)
                ),
                "has_pre_tcr": bool(pre_tcr),
                "has_post_tcr": bool(post_tcr),
                "has_prepost_tcr": bool(pre_tcr and post_tcr),
                "pre_tcr_fractions": join_unique([row.get("sort_or_cell_fraction", "") for row in pre_tcr]),
                "post_tcr_fractions": join_unique([row.get("sort_or_cell_fraction", "") for row in post_tcr]),
                "shared_tcr_fraction": ";".join(sorted(shared_tcr_fractions)),
                "validation_tier": tier,
                "recommended_use": recommendation,
                "response_status": response.get("response_status", "pending_response_source_required"),
                "response_harmonized_ordinal": response.get("response_harmonized_ordinal", ""),
                "response_source": response.get("response_source", "pending: GSE301741 RDS metadata or publication supplementary tables"),
            }
        )

    out_csv = OUT_DIR / "GSE301741_VALIDATION_STRATA.csv"
    write_csv(out_csv, out)

    tier_counts = Counter(row["validation_tier"] for row in out)
    response_counts = Counter(row["response_status"] for row in out)
    candidate_a = [row for row in out if row["validation_tier"] == "A_same_fraction_same_library"]

    md = [
        "# GSE301741 Validation Strata",
        "",
        f"Created: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        f"- Patients represented: {len(out)}",
        f"- Tier A same-fraction/same-library paired scRNA patients: {len(candidate_a)}",
        f"- Patients with any pre/post scRNA: {sum(1 for row in out if row['has_prepost_scrna'])}",
        f"- Patients with any pre/post TCR: {sum(1 for row in out if row['has_prepost_tcr'])}",
        "",
        "## Validation Tier Counts",
        "",
        "| Tier | Patients |",
        "|---|---:|",
    ]
    for tier, count in sorted(tier_counts.items()):
        md.append(f"| {tier} | {count} |")
    md.extend(
        [
            "",
            "## Response Status Counts",
            "",
            "| Response status | Patients |",
            "|---|---:|",
        ]
    )
    for status, count in sorted(response_counts.items()):
        md.append(f"| {status} | {count} |")
    md.extend(
        [
            "",
            "## Tier A Candidates",
            "",
            "| Patient | Site | Shared scRNA fraction/library | Pre TCR | Post TCR | Shared TCR fraction |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in candidate_a:
        md.append(
            "| {patient_id} | {site} | {shared_scrna_fraction_library} | {has_pre_tcr} | {has_post_tcr} | {shared_tcr_fraction} |".format(
                **row
            )
        )
    md.extend(
        [
            "",
            "## Decision",
            "",
            "Use Tier A patients as the cleanest GSE301741 paired single-cell validation denominator.",
            "Tier B can support sensitivity analyses with explicit chemistry control.",
            "Tier C should not be used for abundance-shift claims because sorting fractions differ across timepoints.",
            "Tier D is not paired and should be used only for cross-sectional context.",
            "",
            "Response labels remain pending and must be recovered before response-stratified GSE301741 validation.",
            "",
            "## Files Written",
            "",
            "- `03_rebuild/validation/GSE301741_stratification/GSE301741_VALIDATION_STRATA.csv`",
            "- `03_rebuild/validation/GSE301741_stratification/GSE301741_VALIDATION_STRATA_REPORT.md`",
        ]
    )
    (OUT_DIR / "GSE301741_VALIDATION_STRATA_REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

