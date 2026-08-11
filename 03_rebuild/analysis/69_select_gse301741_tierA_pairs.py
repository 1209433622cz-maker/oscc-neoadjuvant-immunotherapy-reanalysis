from __future__ import annotations

import csv
import os
from collections import Counter
from pathlib import Path


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
MANIFEST_PATH = (
    WORKSPACE
    / "03_rebuild"
    / "manifests"
    / "GSE301741_RAW_ROUTE_SAMPLE_QC.csv"
)
STRATA_PATH = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE301741_response_recovery"
    / "GSE301741_VALIDATION_STRATA_WITH_RESPONSE.csv"
)
DEFAULT_OUTPUT = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE301741_tierA_module_validation"
    / "GSE301741_TIERA_SELECTED_SAMPLE_PAIRS.csv"
)
OUTPUT_PATH = Path(
    os.environ.get("GSE301741_PAIR_SELECTION_OUTPUT", DEFAULT_OUTPUT)
).resolve()

FIELDNAMES = [
    "pair_id",
    "patient_id",
    "response_label",
    "response_status",
    "response_source",
    "provenance_status",
    "timepoint",
    "sort_or_cell_fraction",
    "assay_or_library",
    "geo_accession",
    "extracted_path",
    "n_barcodes_manifest",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def is_true(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def main() -> int:
    manifest = read_csv(MANIFEST_PATH)
    strata = read_csv(STRATA_PATH)

    tier_a = [
        row
        for row in strata
        if row.get("validation_tier") == "A_same_fraction_same_library"
        and row.get("response_label") in {"responder", "non_responder"}
    ]
    if not tier_a:
        raise RuntimeError("No response-labelled Tier-A patients are available.")

    scrna = [
        row
        for row in manifest
        if row.get("modality") == "scRNA"
        and row.get("h5_status") == "ok"
        and is_true(row.get("extracted_exists", ""))
    ]

    selected: list[dict[str, str]] = []
    for patient_row in tier_a:
        patient_id = patient_row.get("patient_id", "")
        shared = [
            value
            for value in patient_row.get(
                "shared_scrna_fraction_library", ""
            ).split(";")
            if value
        ]
        for combination in shared:
            parts = combination.split("/")
            if len(parts) != 2:
                continue
            fraction, library = parts
            candidates = [
                row
                for row in scrna
                if row.get("patient_id") == patient_id
                and row.get("sort_or_cell_fraction") == fraction
                and row.get("assay_or_library") == library
                and row.get("timepoint_harmonized") in {"pre", "post"}
            ]
            if {row.get("timepoint_harmonized") for row in candidates} != {
                "pre",
                "post",
            }:
                continue

            for timepoint in ("pre", "post"):
                sample = next(
                    row
                    for row in candidates
                    if row.get("timepoint_harmonized") == timepoint
                )
                selected.append(
                    {
                        "pair_id": f"{patient_id}__{fraction}__{library}",
                        "patient_id": patient_id,
                        "response_label": patient_row.get("response_label", ""),
                        "response_status": patient_row.get("response_status", ""),
                        "response_source": patient_row.get("response_source", ""),
                        "provenance_status": patient_row.get(
                            "provenance_status", ""
                        ),
                        "timepoint": timepoint,
                        "sort_or_cell_fraction": fraction,
                        "assay_or_library": library,
                        "geo_accession": sample.get("geo_accession", ""),
                        "extracted_path": sample.get("extracted_path", ""),
                        "n_barcodes_manifest": sample.get("n_barcodes", ""),
                    }
                )

    pair_counts = Counter(row["pair_id"] for row in selected)
    selected = [
        row for row in selected if pair_counts.get(row["pair_id"], 0) == 2
    ]
    selected.sort(
        key=lambda row: (
            row["patient_id"],
            row["pair_id"],
            row["timepoint"],
        )
    )
    if not selected:
        raise RuntimeError("No matched Tier-A scRNA H5 pairs were selected.")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(selected)

    print(
        f"Selected {len(selected)} H5 samples across "
        f"{len(pair_counts)} matched pair IDs."
    )
    print(f"Output: {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
