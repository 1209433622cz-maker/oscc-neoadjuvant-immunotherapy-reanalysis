#!/usr/bin/env python
"""Run the locked-family external workflow with the refined parallel manifest."""

from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve()
REBUILD = SCRIPT.parents[1]
BASE_SCRIPT = REBUILD / "analysis" / "61_locked_family_robustness_external_cohorts.py"
MODULE_PATH = (
    REBUILD
    / "results"
    / "refined_module_manifest_sensitivity"
    / "REFINED_MODULE_MANIFEST.csv"
)
OUT_DIR = REBUILD / "validation" / "refined_locked_family_sensitivity"
SOURCE_DIR = OUT_DIR / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def import_base():
    spec = importlib.util.spec_from_file_location("locked_family_refined", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    base = import_base()
    base.MODULE_PATH = MODULE_PATH
    base.OUT_DIR = OUT_DIR
    base.SOURCE_DIR = SOURCE_DIR
    base.SAMPLE_OUT = OUT_DIR / "REFINED_FAMILY_SAMPLE_SCORES.csv"
    base.DELTA_OUT = OUT_DIR / "REFINED_FAMILY_PAIRED_DELTAS.csv"
    base.TEST_OUT = OUT_DIR / "REFINED_FAMILY_TESTS.csv"
    base.RANDOM_OUT = OUT_DIR / "REFINED_FAMILY_MATCHED_RANDOM_EFFECTS.csv"
    base.EXACT_NULL_OUT = OUT_DIR / "GSE179730_REFINED_FAMILY_EXACT_NULL.csv"
    base.COVERAGE_OUT = OUT_DIR / "REFINED_FAMILY_MODULE_COVERAGE.csv"
    base.REPORT_OUT = OUT_DIR / "REFINED_FAMILY_EXTERNAL_SENSITIVITY_REPORT.md"
    base.main()
    print(base.REPORT_OUT)


if __name__ == "__main__":
    main()
