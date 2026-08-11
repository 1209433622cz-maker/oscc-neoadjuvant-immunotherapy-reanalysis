#!/usr/bin/env python
"""Apply the existing response-blind global-PC stress test to the refined family."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
REBUILD = SCRIPT.parents[1]
GLOBAL_SCRIPT = REBUILD / "analysis" / "94_gse281729_global_pc_composition_stress_test.py"
STRESS_SCRIPT = REBUILD / "analysis" / "86_family_composite_global_shift_stress_test.py"
CONFIG = REBUILD / "config" / "gse281729_global_pc_composition.json"
REFINED_DIR = REBUILD / "validation" / "refined_locked_family_sensitivity"
OLD_DIR = REBUILD / "validation" / "locked_family_robustness"
PC_DIR = REBUILD / "validation" / "gse281729_global_pc_composition"
OUT_DIR = REBUILD / "validation" / "refined_family_global_pc_sensitivity"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def import_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bh_adjust(values: pd.Series) -> pd.Series:
    raw = values.to_numpy(float)
    order = np.argsort(raw)
    adjusted = np.empty_like(raw)
    running = 1.0
    for rank_index in range(len(raw) - 1, -1, -1):
        original_index = order[rank_index]
        running = min(running, raw[original_index] * len(raw) / (rank_index + 1))
        adjusted[original_index] = running
    return pd.Series(np.minimum(adjusted, 1.0), index=values.index)


def main() -> None:
    global_analysis = import_script(GLOBAL_SCRIPT, "global_pc_refined")
    stress = import_script(STRESS_SCRIPT, "global_stress_refined")
    config = json.loads(CONFIG.read_text(encoding="utf-8"))

    refined = pd.read_csv(REFINED_DIR / "REFINED_FAMILY_PAIRED_DELTAS.csv")
    refined = refined[refined["cohort"].eq("GSE281729")].copy()
    old = pd.read_csv(OLD_DIR / "LOCKED_FAMILY_PAIRED_DELTAS.csv")
    old = old[old["cohort"].eq("GSE281729")][
        ["patient_id", "scoring_method", "delta_post_minus_pre"]
    ].rename(columns={"delta_post_minus_pre": "old_family_delta"})
    pcs = pd.read_csv(PC_DIR / "GSE281729_GLOBAL_PC_PATIENT_SCORES.csv")
    pcs = pcs[pcs["gene_universe"].eq("primary_nonlocked_mean_ge_1")].drop(
        columns="gene_universe"
    )
    composition = pd.read_csv(PC_DIR / "GSE281729_COMPOSITION_PC_SCORES.csv")
    data = (
        refined.merge(old, on=["patient_id", "scoring_method"], validate="one_to_one")
        .merge(pcs, on="patient_id", validate="many_to_one")
        .merge(composition, on="patient_id", validate="many_to_one")
    )

    rows: list[dict[str, object]] = []
    bootstrap_rows: list[pd.DataFrame] = []
    iterations = int(config["wild_bootstrap_iterations"])
    for method, current in data.groupby("scoring_method", sort=True):
        current = current.sort_values("patient_id")
        for model_spec, extras in config["model_specs"].items():
            seed = stress.stable_seed(
                int(config["wild_bootstrap_seed"]),
                "refined_family",
                method,
                model_spec,
            )
            result, bootstrap = global_analysis.fit_response_model(
                current,
                "delta_post_minus_pre",
                list(extras),
                iterations,
                seed,
                stress,
            )
            rows.append(
                {
                    "scoring_method": method,
                    "model_spec": model_spec,
                    "extra_covariates": ";".join(extras),
                    **result,
                }
            )
            bootstrap.insert(0, "model_spec", model_spec)
            bootstrap.insert(0, "scoring_method", method)
            bootstrap_rows.append(bootstrap)
    models = pd.DataFrame(rows)
    models["model_spec_bh_fdr"] = models.groupby("scoring_method", group_keys=False)[
        "p_value"
    ].apply(bh_adjust)
    models["wild_model_spec_bh_fdr"] = models.groupby(
        "scoring_method", group_keys=False
    )["wild_bootstrap_p"].apply(bh_adjust)
    models.to_csv(OUT_DIR / "REFINED_FAMILY_GLOBAL_PC_MODELS.csv", index=False)
    pd.concat(bootstrap_rows, ignore_index=True).to_csv(
        OUT_DIR / "REFINED_FAMILY_GLOBAL_PC_WILD_BOOTSTRAP.csv", index=False
    )

    correlation_rows: list[dict[str, object]] = []
    for method, current in data.groupby("scoring_method", sort=True):
        for feature in [
            "old_family_delta",
            "global_pc1",
            "global_pc2",
            "composition_pc1",
        ]:
            correlation_rows.append(
                {
                    "scoring_method": method,
                    "feature": feature,
                    "pearson_r": current["delta_post_minus_pre"].corr(current[feature]),
                }
            )
    correlations = pd.DataFrame(correlation_rows)
    correlations.to_csv(
        OUT_DIR / "REFINED_FAMILY_GLOBAL_PC_CORRELATIONS.csv", index=False
    )

    display = models[
        [
            "scoring_method",
            "model_spec",
            "effect",
            "ci95_low",
            "ci95_high",
            "p_value",
            "wild_bootstrap_p",
            "model_spec_bh_fdr",
            "wild_model_spec_bh_fdr",
        ]
    ].copy()
    report = [
        "# Refined-Family Global-PC Sensitivity",
        "",
        "- The refined 16-module manifest is a post-freeze sensitivity family.",
        "- Global and composition PCs are reused from the response-blind GSE281729 analysis.",
        "- Models retain HPV, second-drug and response-adaptive dose/timing covariates.",
        f"- Wild bootstrap uses {iterations:,} HC3-studentized Rademacher draws per model.",
        "",
        "## Models",
        "",
        display.to_string(index=False),
        "",
        "## Correlations",
        "",
        correlations.to_string(index=False),
        "",
        "The refined family may be described as annotation-definition robust only if it remains "
        "highly correlated with the original frozen family and is similarly attenuated by the "
        "response-blind global PC. This sensitivity does not replace the original freeze.",
    ]
    (OUT_DIR / "REFINED_FAMILY_GLOBAL_PC_SENSITIVITY_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(OUT_DIR / "REFINED_FAMILY_GLOBAL_PC_SENSITIVITY_REPORT.md")


if __name__ == "__main__":
    main()
