from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
REBUILD = ROOT / "03_rebuild"
BASE_SCRIPT = REBUILD / "analysis" / "61_locked_family_robustness_external_cohorts.py"
CONFIG_PATH = REBUILD / "config" / "family_composite_stress_test.json"
OUT_DIR = REBUILD / "validation" / "family_composite_stress_test"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SCORES_OUT = OUT_DIR / "FAMILY_SCORE_VARIANTS.csv"
DELTAS_OUT = OUT_DIR / "FAMILY_VARIANT_DELTAS.csv"
TESTS_OUT = OUT_DIR / "FAMILY_VARIANT_TESTS.csv"
BOOTSTRAP_OUT = OUT_DIR / "GSE281729_FAMILY_VARIANT_WILD_BOOTSTRAP.csv"
BACKGROUND_OUT = OUT_DIR / "MATCHED_BACKGROUND_CONTROL_DELTAS.csv"
WEIGHTS_OUT = OUT_DIR / "FAMILY_VARIANT_WEIGHT_AUDIT.csv"
REPORT_OUT = OUT_DIR / "FAMILY_COMPOSITE_STRESS_TEST_REPORT.md"


def load_base_module():
    spec = importlib.util.spec_from_file_location("locked_family", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_config() -> dict:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    if config["variant_reselection_allowed"]:
        raise RuntimeError("Stress-test variants must be fixed without reselection")
    return config


def score_components(
    expression: pd.DataFrame,
    modules: pd.DataFrame,
    cohort: str,
    base,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    finite = expression.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    finite = finite.loc[finite.std(axis=1, ddof=0) > 0]
    arrays = {
        "z_score": finite.sub(finite.mean(axis=1), axis=0)
        .div(finite.std(axis=1, ddof=0), axis=0)
        .to_numpy(float),
        "rank_mean": finite.rank(axis=0, method="average", pct=True).to_numpy(float),
    }
    gene_index = finite.index
    present_modules: dict[str, np.ndarray] = {}
    lineage: dict[str, str] = {}
    for _, row in modules.iterrows():
        genes = [gene for gene in row["genes"] if gene in gene_index]
        if not genes:
            raise RuntimeError(f"No measured genes for {cohort}/{row['signature']}")
        present_modules[str(row["signature"])] = gene_index.get_indexer(genes)
        lineage[str(row["signature"])] = str(row["target_lineage"])

    memberships = Counter(
        int(position)
        for positions in present_modules.values()
        for position in positions
    )
    unique_positions = np.asarray(sorted(memberships), dtype=int)
    no_union = {
        name: positions
        for name, positions in present_modules.items()
        if not name.endswith("_union_core")
    }
    hallmark = {
        name: positions
        for name, positions in present_modules.items()
        if "_LE_" in name and not name.endswith("_union_core")
    }
    dynamic = {
        name: positions
        for name, positions in present_modules.items()
        if "_LE_" not in name
    }

    def module_mean(matrix: np.ndarray, selected: dict[str, np.ndarray]) -> np.ndarray:
        return np.vstack(
            [matrix[positions, :].mean(axis=0) for positions in selected.values()]
        ).mean(axis=0)

    def inverse_membership_mean(matrix: np.ndarray) -> np.ndarray:
        module_values = []
        for positions in present_modules.values():
            weights = np.asarray([1.0 / memberships[int(pos)] for pos in positions])
            module_values.append(np.average(matrix[positions, :], axis=0, weights=weights))
        return np.vstack(module_values).mean(axis=0)

    bins = base.expression_bins(finite)
    locked = set(unique_positions.tolist())
    locked_bin_counts = Counter(int(bins.iloc[position]) for position in unique_positions)
    background_by_bin = {
        bin_id: np.asarray(
            [
                position
                for position in range(len(gene_index))
                if position not in locked and int(bins.iloc[position]) == bin_id
            ],
            dtype=int,
        )
        for bin_id in sorted(locked_bin_counts)
    }
    if any(len(values) == 0 for values in background_by_bin.values()):
        raise RuntimeError(f"{cohort} has an empty matched-background expression bin")

    sample_rows: list[dict[str, object]] = []
    background_rows: list[dict[str, object]] = []
    weight_rows: list[dict[str, object]] = []
    variants: dict[str, dict[str, np.ndarray]] = {}
    for method, matrix in arrays.items():
        background = sum(
            locked_bin_counts[bin_id] * matrix[background_by_bin[bin_id], :].mean(axis=0)
            for bin_id in locked_bin_counts
        ) / sum(locked_bin_counts.values())
        variants[method] = {
            "module_mean_16": module_mean(matrix, present_modules),
            "unique_gene_equal": matrix[unique_positions, :].mean(axis=0),
            "inverse_membership_module_mean": inverse_membership_mean(matrix),
            "no_union_module_mean": module_mean(matrix, no_union),
            "hallmark_only_module_mean": module_mean(matrix, hallmark),
            "dynamic_only_module_mean": module_mean(matrix, dynamic),
        }
        for sample_id, value in zip(finite.columns, background):
            background_rows.append(
                {
                    "cohort": cohort,
                    "sample_id": sample_id,
                    "scoring_method": method,
                    "matched_background_score": float(value),
                    "n_background_genes": sum(
                        len(values) for values in background_by_bin.values()
                    ),
                }
            )
        for variant, values in variants[method].items():
            for sample_id, value in zip(finite.columns, values):
                sample_rows.append(
                    {
                        "cohort": cohort,
                        "sample_id": sample_id,
                        "scoring_method": method,
                        "score_variant": variant,
                        "family_score": float(value),
                    }
                )

    current_weights: Counter[int] = Counter()
    for positions in present_modules.values():
        for position in positions:
            current_weights[int(position)] += 1.0 / (
                len(present_modules) * len(positions)
            )
    inverse_weights: Counter[int] = Counter()
    for positions in present_modules.values():
        raw = np.asarray([1.0 / memberships[int(pos)] for pos in positions])
        normalized = raw / raw.sum() / len(present_modules)
        for position, weight in zip(positions, normalized):
            inverse_weights[int(position)] += float(weight)
    for variant, weights in {
        "module_mean_16": current_weights,
        "unique_gene_equal": Counter(
            {int(position): 1.0 / len(unique_positions) for position in unique_positions}
        ),
        "inverse_membership_module_mean": inverse_weights,
    }.items():
        values = np.asarray(list(weights.values()), dtype=float)
        weight_rows.append(
            {
                "cohort": cohort,
                "score_variant": variant,
                "n_unique_genes": len(weights),
                "weight_sum": float(values.sum()),
                "minimum_gene_weight": float(values.min()),
                "median_gene_weight": float(np.median(values)),
                "maximum_gene_weight": float(values.max()),
                "max_to_min_weight_ratio": float(values.max() / values.min()),
                "effective_gene_number": float(1.0 / np.sum(values**2)),
                "genes_with_multiple_module_memberships": sum(
                    memberships[position] > 1 for position in weights
                ),
                "maximum_membership": max(memberships[position] for position in weights),
            }
        )
    return (
        pd.DataFrame(sample_rows),
        pd.DataFrame(background_rows),
        pd.DataFrame(weight_rows),
    )


def paired_deltas(
    scores: pd.DataFrame,
    background: pd.DataFrame,
    annotation: pd.DataFrame,
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_meta = scores.merge(annotation, on="sample_id", how="left", validate="many_to_one")
    index = ["patient_id", "scoring_method", "score_variant"]
    wide = score_meta.pivot_table(
        index=index,
        columns="timepoint",
        values="family_score",
        aggfunc="mean",
    ).reset_index()
    wide["delta_post_minus_pre"] = wide["post"] - wide["pre"]

    background_meta = background.merge(
        annotation, on="sample_id", how="left", validate="many_to_one"
    )
    background_wide = background_meta.pivot_table(
        index=["patient_id", "scoring_method"],
        columns="timepoint",
        values="matched_background_score",
        aggfunc="mean",
    ).reset_index()
    background_wide["background_delta_post_minus_pre"] = (
        background_wide["post"] - background_wide["pre"]
    )
    metadata_columns = [
        column
        for column in [
            "patient_id",
            "response_binary",
            "response_harmonized_ordinal",
            "response_ord_num",
            "hpv",
            "second_drug",
            "doses",
            "primary_response_percent",
        ]
        if column in annotation.columns
    ]
    metadata = annotation[metadata_columns].drop_duplicates("patient_id")
    wide = wide.merge(
        background_wide[
            ["patient_id", "scoring_method", "background_delta_post_minus_pre"]
        ],
        on=["patient_id", "scoring_method"],
        how="left",
        validate="many_to_one",
    ).merge(metadata, on="patient_id", how="left", validate="many_to_one")
    wide.insert(0, "cohort", cohort)
    background_wide = background_wide.merge(
        metadata, on="patient_id", how="left", validate="many_to_one"
    )
    background_wide.insert(0, "cohort", cohort)
    return wide, background_wide


def exact_effect(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    responder = labels == "responder"
    n_responder = int(responder.sum())
    observed = float(values[responder].mean() - values[~responder].mean())
    null = []
    for indices in itertools.combinations(range(len(values)), n_responder):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(indices)] = True
        null.append(float(values[mask].mean() - values[~mask].mean()))
    null_values = np.asarray(null)
    p_value = float(np.mean(np.abs(null_values) >= abs(observed) - 1e-12))
    return observed, p_value, len(null_values)


def gse179_tests(deltas: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (method, variant), subset in deltas.groupby(
        ["scoring_method", "score_variant"]
    ):
        current = subset.sort_values("patient_id").reset_index(drop=True)
        labels = current["response_binary"].to_numpy()
        y = current["delta_post_minus_pre"].to_numpy(float)
        background = current["background_delta_post_minus_pre"].to_numpy(float)
        for adjusted in [False, True]:
            tested = y
            background_coefficient = math.nan
            if adjusted:
                design = sm.add_constant(background, has_constant="add")
                fit = sm.OLS(y, design).fit()
                tested = fit.resid
                background_coefficient = float(fit.params[1])
            effect, p_value, assignments = exact_effect(tested, labels)
            rows.append(
                {
                    "cohort": "GSE179730",
                    "scoring_method": method,
                    "score_variant": variant,
                    "background_adjusted": adjusted,
                    "analysis": (
                        "reduced_model_background_residual_exact"
                        if adjusted
                        else "unadjusted_response_exact"
                    ),
                    "n_patients": len(current),
                    "effect": effect,
                    "std_error": math.nan,
                    "ci95_low": math.nan,
                    "ci95_high": math.nan,
                    "p_value": p_value,
                    "background_delta_coefficient": background_coefficient,
                    "family_background_delta_correlation": float(
                        np.corrcoef(y, background)[0, 1]
                    ),
                    "design_condition_number": math.nan,
                    "exact_assignments": assignments,
                }
            )
    return rows


def gse281_design(
    current: pd.DataFrame,
    include_background: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "delta_post_minus_pre",
        "response_ord_num",
        "hpv",
        "second_drug",
        "doses",
    ]
    if include_background:
        columns.append("background_delta_post_minus_pre")
    model_data = current[columns].replace([np.inf, -np.inf], np.nan).dropna()
    parts = [model_data[["response_ord_num"]].astype(float)]
    for column in ["hpv", "second_drug"]:
        if model_data[column].nunique() > 1:
            parts.append(
                pd.get_dummies(
                    model_data[column].astype(str),
                    prefix=column,
                    drop_first=True,
                    dtype=float,
                )
            )
    parts.append(model_data[["doses"]].astype(float))
    if include_background:
        parts.append(model_data[["background_delta_post_minus_pre"]].astype(float))
    return model_data, sm.add_constant(pd.concat(parts, axis=1), has_constant="add")


def stable_seed(base_seed: int, *parts: str) -> int:
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).digest()
    return (base_seed + int.from_bytes(digest[:4], "little")) % (2**32 - 1)


def wild_bootstrap_p(
    y: np.ndarray,
    design: pd.DataFrame,
    iterations: int,
    seed: int,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    response_index = list(design.columns).index("response_ord_num")
    x = design.to_numpy(float)
    full_pinv = np.linalg.pinv(x)
    full_leverage = np.diag(x @ full_pinv)
    observed_beta = full_pinv @ y
    observed_residual = y - x @ observed_beta
    response_projection = full_pinv[response_index, :]
    observed_se = float(
        np.sqrt(
            np.sum(
                response_projection**2
                * (observed_residual / np.maximum(1.0 - full_leverage, 1e-8)) ** 2
            )
        )
    )
    observed_t = float(observed_beta[response_index] / observed_se)

    reduced = np.delete(x, response_index, axis=1)
    reduced_fit = sm.OLS(y, reduced).fit()
    leverage = reduced_fit.get_influence().hat_matrix_diag
    scaled_residual = reduced_fit.resid / np.sqrt(np.maximum(1.0 - leverage, 1e-8))
    rng = np.random.default_rng(seed)
    multipliers = rng.choice(
        np.asarray([-1.0, 1.0]),
        size=(len(y), iterations),
        replace=True,
    )
    y_star = reduced_fit.fittedvalues[:, None] + scaled_residual[:, None] * multipliers
    coefficients = full_pinv @ y_star
    bootstrap = coefficients[response_index, :]
    bootstrap_residual = y_star - x @ coefficients
    bootstrap_se = np.sqrt(
        np.sum(
            response_projection[:, None] ** 2
            * (
                bootstrap_residual
                / np.maximum(1.0 - full_leverage[:, None], 1e-8)
            )
            ** 2,
            axis=0,
        )
    )
    bootstrap_t = bootstrap / bootstrap_se
    p_value = float(
        (1 + np.sum(np.abs(bootstrap_t) >= abs(observed_t) - 1e-12))
        / (iterations + 1)
    )
    return p_value, observed_t, bootstrap, bootstrap_t


def gse281_tests(
    deltas: pd.DataFrame,
    config: dict,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    iterations = int(config["wild_bootstrap_iterations"])
    for (method, variant), subset in deltas.groupby(
        ["scoring_method", "score_variant"]
    ):
        current = subset[
            subset["response_harmonized_ordinal"].isin(["Low", "Medium", "High"])
        ].sort_values("patient_id")
        for adjusted in [False, True]:
            model_data, design = gse281_design(current, adjusted)
            y = model_data["delta_post_minus_pre"].to_numpy(float)
            fit = sm.OLS(y, design.astype(float)).fit(cov_type="HC3")
            term = "response_ord_num"
            ci = fit.conf_int(alpha=0.05).loc[term]
            seed = stable_seed(
                int(config["wild_bootstrap_seed"]),
                method,
                variant,
                str(adjusted),
            )
            bootstrap_p, bootstrap_observed_t, bootstrap, bootstrap_t = wild_bootstrap_p(
                y,
                design.astype(float),
                iterations,
                seed,
            )
            background_coefficient = (
                float(fit.params["background_delta_post_minus_pre"])
                if adjusted
                else math.nan
            )
            rows.append(
                {
                    "cohort": "GSE281729",
                    "scoring_method": method,
                    "score_variant": variant,
                    "background_adjusted": adjusted,
                    "analysis": (
                        "hpv_second_drug_doses_background_HC3"
                        if adjusted
                        else "hpv_second_drug_doses_HC3"
                    ),
                    "n_patients": len(model_data),
                    "effect": float(fit.params[term]),
                    "std_error": float(fit.bse[term]),
                    "ci95_low": float(ci.iloc[0]),
                    "ci95_high": float(ci.iloc[1]),
                    "p_value": float(fit.pvalues[term]),
                    "wild_bootstrap_p": bootstrap_p,
                    "wild_bootstrap_observed_t": bootstrap_observed_t,
                    "background_delta_coefficient": background_coefficient,
                    "family_background_delta_correlation": float(
                        model_data["delta_post_minus_pre"].corr(
                            model_data.get(
                                "background_delta_post_minus_pre",
                                current.loc[
                                    model_data.index,
                                    "background_delta_post_minus_pre",
                                ],
                            )
                        )
                    ),
                    "design_condition_number": float(np.linalg.cond(design)),
                    "exact_assignments": math.nan,
                }
            )
            for iteration, (value, t_value) in enumerate(
                zip(bootstrap, bootstrap_t),
                start=1,
            ):
                bootstrap_rows.append(
                    {
                        "scoring_method": method,
                        "score_variant": variant,
                        "background_adjusted": adjusted,
                        "iteration": iteration,
                        "random_seed": seed,
                        "response_coefficient_under_null": float(value),
                        "response_hc3_t_under_null": float(t_value),
                    }
                )
    return rows, pd.DataFrame(bootstrap_rows)


def bh_adjust(values: pd.Series) -> pd.Series:
    raw = values.to_numpy(float)
    order = np.argsort(raw)
    adjusted = np.empty_like(raw)
    running = 1.0
    for rank_index in range(len(raw) - 1, -1, -1):
        original_index = order[rank_index]
        candidate = raw[original_index] * len(raw) / (rank_index + 1)
        running = min(running, candidate)
        adjusted[original_index] = running
    return pd.Series(np.minimum(adjusted, 1.0), index=values.index)


def write_report(tests: pd.DataFrame, weights: pd.DataFrame) -> None:
    gse281 = tests[
        (tests["cohort"] == "GSE281729")
        & (tests["scoring_method"] == "z_score")
    ][
        [
            "score_variant",
            "background_adjusted",
            "effect",
            "p_value",
            "wild_bootstrap_p",
            "variant_bh_fdr",
            "family_background_delta_correlation",
        ]
    ]
    gse179 = tests[
        (tests["cohort"] == "GSE179730")
        & (tests["scoring_method"] == "z_score")
    ][
        [
            "score_variant",
            "background_adjusted",
            "effect",
            "p_value",
            "variant_bh_fdr",
            "family_background_delta_correlation",
        ]
    ]

    def markdown(frame: pd.DataFrame) -> str:
        def render(value: object) -> str:
            if pd.isna(value):
                return ""
            if isinstance(value, (float, np.floating)):
                return f"{float(value):.6g}"
            return str(value).replace("|", r"\|")

        columns = [str(column) for column in frame.columns]
        lines = [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join(["---"] * len(columns)) + " |",
        ]
        lines.extend(
            "| " + " | ".join(render(value) for value in row) + " |"
            for row in frame.itertuples(index=False, name=None)
        )
        return "\n".join(lines)

    report = [
        "# Family-Composite and Global-Shift Stress Test",
        "",
        "## Purpose",
        "",
        (
            "This post hoc stress test evaluates whether the external locked-family "
            "association depends on repeated gene membership, union modules or a "
            "broad expression-decile-matched transcriptomic shift."
        ),
        "",
        "## Weight audit",
        "",
        markdown(weights),
        "",
        "## GSE281729 primary z-score sensitivities",
        "",
        markdown(gse281),
        "",
        "## GSE179730 primary z-score sensitivities",
        "",
        markdown(gse179),
        "",
        "## Interpretation rules",
        "",
        "- All six variants are reported without outcome-driven selection.",
        "- Rank scoring is a scale sensitivity, not independent replication.",
        (
            "- Background adjustment uses a deterministic non-locked transcriptome "
            "control matched to the locked unique-gene expression-decile distribution."
        ),
        (
            "- Wild-bootstrap P values use 9,999 Rademacher draws under the reduced "
            "covariate model with HC2 leverage-scaled residuals and an HC3-studentized "
            "response statistic."
        ),
        (
            "- These post hoc analyses can constrain the family claim but cannot "
            "remove the source study's response-adaptive exposure and timing limitation."
        ),
    ]
    REPORT_OUT.write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> int:
    config = load_config()
    base = load_base_module()
    modules = base.load_modules()

    all_scores = []
    all_deltas = []
    all_background = []
    all_weights = []
    cohort_payload = [
        ("GSE179730", *base.load_gse179730()),
        ("GSE281729", *base.load_gse281729()),
    ]
    for cohort, expression, annotation in cohort_payload:
        scores, background, weights = score_components(
            expression,
            modules,
            cohort,
            base,
        )
        deltas, background_deltas = paired_deltas(
            scores,
            background,
            annotation,
            cohort,
        )
        all_scores.append(scores)
        all_deltas.append(deltas)
        all_background.append(background_deltas)
        all_weights.append(weights)

    scores = pd.concat(all_scores, ignore_index=True)
    deltas = pd.concat(all_deltas, ignore_index=True)
    background = pd.concat(all_background, ignore_index=True)
    weights = pd.concat(all_weights, ignore_index=True)

    test_rows = gse179_tests(deltas[deltas["cohort"] == "GSE179730"])
    gse281_rows, bootstrap = gse281_tests(
        deltas[deltas["cohort"] == "GSE281729"],
        config,
    )
    tests = pd.DataFrame(test_rows + gse281_rows)
    tests["variant_bh_fdr"] = tests.groupby(
        ["cohort", "scoring_method", "background_adjusted"],
        group_keys=False,
    )["p_value"].apply(bh_adjust)
    tests["wild_bootstrap_variant_bh_fdr"] = math.nan
    mask = tests["wild_bootstrap_p"].notna()
    tests.loc[mask, "wild_bootstrap_variant_bh_fdr"] = (
        tests.loc[mask]
        .groupby(
            ["cohort", "scoring_method", "background_adjusted"],
            group_keys=False,
        )["wild_bootstrap_p"]
        .apply(bh_adjust)
    )

    expected_variants = set(config["score_variants"])
    if set(tests["score_variant"]) != expected_variants:
        raise RuntimeError("Score-variant set differs from the frozen configuration")
    if len(tests) != 48:
        raise RuntimeError(f"Expected 48 test rows, observed {len(tests)}")
    if len(bootstrap) != 24 * int(config["wild_bootstrap_iterations"]):
        raise RuntimeError("Wild-bootstrap row count mismatch")

    scores.to_csv(SCORES_OUT, index=False)
    deltas.to_csv(DELTAS_OUT, index=False)
    tests.to_csv(TESTS_OUT, index=False)
    bootstrap.to_csv(BOOTSTRAP_OUT, index=False)
    background.to_csv(BACKGROUND_OUT, index=False)
    weights.to_csv(WEIGHTS_OUT, index=False)
    write_report(tests, weights)
    print(TESTS_OUT)
    print(REPORT_OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
