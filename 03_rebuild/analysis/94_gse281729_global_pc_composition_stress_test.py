from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
REBUILD = ROOT / "03_rebuild"
BASE_SCRIPT = REBUILD / "analysis" / "61_locked_family_robustness_external_cohorts.py"
STRESS_SCRIPT = REBUILD / "analysis" / "86_family_composite_global_shift_stress_test.py"
CONFIG_PATH = REBUILD / "config" / "gse281729_global_pc_composition.json"
OUT_DIR = REBUILD / "validation" / "gse281729_global_pc_composition"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_QC_OUT = OUT_DIR / "GSE281729_PROCESSED_EXPRESSION_SAMPLE_QC.csv"
PAIRED_QC_OUT = OUT_DIR / "GSE281729_PROCESSED_EXPRESSION_PAIRED_QC_DELTAS.csv"
PC_QC_CORRELATIONS_OUT = OUT_DIR / "GSE281729_GLOBAL_PC_QC_CORRELATIONS.csv"
PC_SCORES_OUT = OUT_DIR / "GSE281729_GLOBAL_PC_PATIENT_SCORES.csv"
PC_LOADINGS_OUT = OUT_DIR / "GSE281729_GLOBAL_PC_GENE_LOADINGS.csv"
PC_VARIANCE_OUT = OUT_DIR / "GSE281729_GLOBAL_PC_VARIANCE.csv"
PC_STABILITY_OUT = OUT_DIR / "GSE281729_GLOBAL_PC_STABILITY.csv"
MCP_COVERAGE_OUT = OUT_DIR / "GSE281729_MCP_COUNTER_MARKER_COVERAGE.csv"
MCP_SAMPLE_OUT = OUT_DIR / "GSE281729_MCP_COUNTER_SAMPLE_SCORES.csv"
MCP_DELTA_OUT = OUT_DIR / "GSE281729_MCP_COUNTER_PAIRED_DELTAS.csv"
COMPOSITION_PC_OUT = OUT_DIR / "GSE281729_COMPOSITION_PC_SCORES.csv"
COMPOSITION_LOADINGS_OUT = OUT_DIR / "GSE281729_COMPOSITION_PC_LOADINGS.csv"
FAMILY_DELTA_OUT = OUT_DIR / "GSE281729_FAMILY_GLOBAL_PC_COMPOSITION_DELTAS.csv"
FAMILY_MODELS_OUT = OUT_DIR / "GSE281729_FAMILY_GLOBAL_PC_COMPOSITION_MODELS.csv"
FAMILY_BOOTSTRAP_OUT = OUT_DIR / "GSE281729_FAMILY_GLOBAL_PC_COMPOSITION_WILD_BOOTSTRAP.csv"
MCP_MODELS_OUT = OUT_DIR / "GSE281729_MCP_COUNTER_RESPONSE_MODELS.csv"
PC_MODELS_OUT = OUT_DIR / "GSE281729_PC_RESPONSE_MODELS.csv"
CORRELATIONS_OUT = OUT_DIR / "GSE281729_FAMILY_PC_COMPOSITION_CORRELATIONS.csv"
REPORT_OUT = OUT_DIR / "GSE281729_GLOBAL_PC_COMPOSITION_REPORT.md"


def import_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def load_config() -> dict:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    if config["response_dependent_feature_selection_allowed"]:
        raise RuntimeError("Response-dependent feature selection must remain disabled")
    signature = ROOT / config["composition_signature_file"]
    if sha256(signature).lower() != config["composition_signature_sha256"].lower():
        raise RuntimeError("MCP-counter signature SHA256 does not match frozen configuration")
    return config


def paired_inputs(
    expression: pd.DataFrame,
    annotation: pd.DataFrame,
    expected_patients: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = (
        annotation[annotation["timepoint"].isin(["pre", "post"])]
        .groupby(["patient_id", "timepoint"])
        .size()
        .unstack(fill_value=0)
    )
    eligible = counts[(counts.get("pre", 0) == 1) & (counts.get("post", 0) == 1)].index
    paired_annotation = annotation[
        annotation["patient_id"].isin(eligible)
        & annotation["timepoint"].isin(["pre", "post"])
        & annotation["response_harmonized_ordinal"].isin(["Low", "Medium", "High"])
    ].copy()
    retained = paired_annotation["patient_id"].nunique()
    if retained != expected_patients:
        raise RuntimeError(f"Expected {expected_patients} paired patients; found {retained}")
    paired_annotation = paired_annotation.sort_values(
        ["patient_id", "timepoint"]
    ).reset_index(drop=True)
    paired_expression = expression[paired_annotation["sample_id"].tolist()].copy()
    sample_map = paired_annotation.pivot(
        index="patient_id", columns="timepoint", values="sample_id"
    )
    deltas = pd.DataFrame(
        {
            patient: paired_expression[sample_map.loc[patient, "post"]]
            - paired_expression[sample_map.loc[patient, "pre"]]
            for patient in sample_map.index
        }
    )
    return paired_expression, paired_annotation, deltas


def sample_qc(expression: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample_id in expression.columns:
        values = expression[sample_id].to_numpy(float)
        meta = annotation[annotation["sample_id"].eq(sample_id)].iloc[0]
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": meta["patient_id"],
                "timepoint": meta["timepoint"],
                "response_harmonized_ordinal": meta["response_harmonized_ordinal"],
                "response_ord_num": meta["response_ord_num"],
                "hpv": meta["hpv"],
                "second_drug": meta["second_drug"],
                "doses": meta["doses"],
                "n_genes_total": len(values),
                "n_genes_positive": int(np.sum(values > 0)),
                "mean_logtpm": float(np.mean(values)),
                "median_logtpm": float(np.median(values)),
                "p10_logtpm": float(np.quantile(values, 0.10)),
                "p90_logtpm": float(np.quantile(values, 0.90)),
            }
        )
    return pd.DataFrame(rows)


def paired_qc_deltas(qc: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "n_genes_positive",
        "mean_logtpm",
        "median_logtpm",
        "p10_logtpm",
        "p90_logtpm",
    ]
    result = pd.DataFrame({"patient_id": sorted(qc["patient_id"].unique())})
    for metric in metrics:
        wide = qc.pivot(index="patient_id", columns="timepoint", values=metric)
        delta = (wide["post"] - wide["pre"]).rename(f"delta_{metric}")
        result = result.merge(
            delta.reset_index(), on="patient_id", how="left", validate="one_to_one"
        )
    return result


def deterministic_pca(
    delta_matrix: pd.DataFrame,
    genes: list[str],
    universe: str,
    n_components: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = delta_matrix.loc[genes].copy()
    means = selected.mean(axis=1)
    sds = selected.std(axis=1, ddof=0)
    standardized = selected.sub(means, axis=0).div(sds, axis=0)
    x = standardized.T.to_numpy(float)
    u, singular, vt = np.linalg.svd(x, full_matrices=False)
    total = float(np.sum(singular**2))
    score_rows = {"patient_id": standardized.columns.tolist(), "gene_universe": universe}
    loading_frames = []
    variance_rows = []
    for index in range(n_components):
        scores = u[:, index] * singular[index]
        loadings = vt[index, :].copy()
        anchor = int(np.argmax(np.abs(loadings)))
        if loadings[anchor] < 0:
            scores *= -1
            loadings *= -1
        scaled_scores = (scores - scores.mean()) / scores.std(ddof=0)
        pc = f"global_pc{index + 1}"
        score_rows[pc] = scaled_scores
        loading_frames.append(
            pd.DataFrame(
                {
                    "gene_universe": universe,
                    "gene": standardized.index,
                    "pc": pc,
                    "loading": loadings,
                    "absolute_loading": np.abs(loadings),
                }
            )
        )
        variance_rows.append(
            {
                "gene_universe": universe,
                "pc": pc,
                "n_genes": len(genes),
                "explained_variance_ratio": float(singular[index] ** 2 / total),
                "singular_value": float(singular[index]),
                "orientation_anchor_gene": standardized.index[anchor],
            }
        )
    return (
        pd.DataFrame(score_rows),
        pd.concat(loading_frames, ignore_index=True),
        pd.DataFrame(variance_rows),
    )


def mcp_counter_scores(
    expression: pd.DataFrame,
    annotation: pd.DataFrame,
    signature_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    signatures = pd.read_csv(signature_path, sep="\t")
    signatures["HUGO symbols"] = signatures["HUGO symbols"].astype(str)
    coverage_rows = []
    score_rows = []
    present_by_population: dict[str, list[str]] = {}
    for population, current in signatures.groupby("Cell population", sort=False):
        defined = list(dict.fromkeys(current["HUGO symbols"].tolist()))
        present = [gene for gene in defined if gene in expression.index]
        if not present:
            raise RuntimeError(f"No measured MCP-counter markers for {population}")
        present_by_population[population] = present
        coverage_rows.append(
            {
                "population": population,
                "markers_defined": len(defined),
                "markers_present": len(present),
                "coverage_fraction": len(present) / len(defined),
                "present_markers": ";".join(present),
                "missing_markers": ";".join(sorted(set(defined) - set(present))),
            }
        )
        values = expression.loc[present].mean(axis=0)
        for sample_id, value in values.items():
            score_rows.append(
                {
                    "sample_id": sample_id,
                    "population": population,
                    "mcp_counter_score": float(value),
                }
            )
    scores = pd.DataFrame(score_rows).merge(
        annotation[
            [
                "sample_id",
                "patient_id",
                "timepoint",
                "response_harmonized_ordinal",
                "response_ord_num",
                "hpv",
                "second_drug",
                "doses",
            ]
        ],
        on="sample_id",
        how="left",
        validate="many_to_one",
    )
    wide = scores.pivot_table(
        index=["patient_id", "population"],
        columns="timepoint",
        values="mcp_counter_score",
        aggfunc="mean",
    ).reset_index()
    wide["delta_post_minus_pre"] = wide["post"] - wide["pre"]
    metadata = annotation[
        [
            "patient_id",
            "response_harmonized_ordinal",
            "response_ord_num",
            "hpv",
            "second_drug",
            "doses",
        ]
    ].drop_duplicates("patient_id")
    deltas = wide.merge(metadata, on="patient_id", how="left", validate="many_to_one")
    return pd.DataFrame(coverage_rows), scores, deltas


def composition_pca(
    mcp_deltas: pd.DataFrame,
    n_components: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = mcp_deltas.pivot(
        index="patient_id", columns="population", values="delta_post_minus_pre"
    )
    standardized = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
    u, singular, vt = np.linalg.svd(standardized.to_numpy(float), full_matrices=False)
    total = float(np.sum(singular**2))
    scores = pd.DataFrame({"patient_id": standardized.index})
    loadings = []
    variance = []
    for index in range(n_components):
        current_scores = u[:, index] * singular[index]
        current_loadings = vt[index, :].copy()
        anchor = int(np.argmax(np.abs(current_loadings)))
        if current_loadings[anchor] < 0:
            current_scores *= -1
            current_loadings *= -1
        scores[f"composition_pc{index + 1}"] = (
            current_scores - current_scores.mean()
        ) / current_scores.std(ddof=0)
        for population, loading in zip(standardized.columns, current_loadings):
            loadings.append(
                {
                    "pc": f"composition_pc{index + 1}",
                    "population": population,
                    "loading": float(loading),
                    "absolute_loading": float(abs(loading)),
                }
            )
        variance.append(
            {
                "pc": f"composition_pc{index + 1}",
                "explained_variance_ratio": float(singular[index] ** 2 / total),
                "singular_value": float(singular[index]),
                "orientation_anchor_population": standardized.columns[anchor],
            }
        )
    return scores, pd.DataFrame(loadings), pd.DataFrame(variance)


def model_design(
    data: pd.DataFrame,
    outcome: str,
    extra_covariates: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        outcome,
        "response_ord_num",
        "hpv",
        "second_drug",
        "doses",
        *extra_covariates,
    ]
    model_data = data[columns].replace([np.inf, -np.inf], np.nan).dropna()
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
    if extra_covariates:
        parts.append(model_data[extra_covariates].astype(float))
    design = sm.add_constant(pd.concat(parts, axis=1), has_constant="add")
    return model_data, design


def fit_response_model(
    data: pd.DataFrame,
    outcome: str,
    extra_covariates: list[str],
    iterations: int,
    seed: int,
    stress,
) -> tuple[dict[str, object], pd.DataFrame]:
    model_data, design = model_design(data, outcome, extra_covariates)
    y = model_data[outcome].to_numpy(float)
    fit = sm.OLS(y, design.astype(float)).fit(cov_type="HC3")
    term = "response_ord_num"
    ci = fit.conf_int(alpha=0.05).loc[term]
    wild_p, observed_t, bootstrap, bootstrap_t = stress.wild_bootstrap_p(
        y, design.astype(float), iterations, seed
    )
    response_index = list(design.columns).index(term)
    reduced = np.delete(design.to_numpy(float), response_index, axis=1)
    reduced_fit = sm.OLS(y, reduced).fit()
    full_fit = sm.OLS(y, design.astype(float)).fit()
    partial_r2 = max(
        0.0,
        float((np.sum(reduced_fit.resid**2) - np.sum(full_fit.resid**2)) / np.sum(reduced_fit.resid**2)),
    )
    extras = {
        column: float(fit.params[column])
        for column in extra_covariates
        if column in fit.params.index
    }
    result = {
        "n_patients": len(model_data),
        "effect": float(fit.params[term]),
        "std_error": float(fit.bse[term]),
        "ci95_low": float(ci.iloc[0]),
        "ci95_high": float(ci.iloc[1]),
        "p_value": float(fit.pvalues[term]),
        "wild_bootstrap_p": wild_p,
        "wild_bootstrap_observed_t": observed_t,
        "response_partial_r2": partial_r2,
        "model_r2": float(full_fit.rsquared),
        "design_condition_number": float(np.linalg.cond(design.to_numpy(float))),
        "extra_covariate_coefficients": json.dumps(extras, sort_keys=True),
    }
    bootstrap_frame = pd.DataFrame(
        {
            "iteration": np.arange(1, iterations + 1),
            "random_seed": seed,
            "response_coefficient_under_null": bootstrap,
            "response_hc3_t_under_null": bootstrap_t,
        }
    )
    return result, bootstrap_frame


def family_models(
    family_deltas: pd.DataFrame,
    config: dict,
    stress,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    bootstrap_frames = []
    iterations = int(config["wild_bootstrap_iterations"])
    for (method, variant), current in family_deltas.groupby(
        ["scoring_method", "score_variant"], sort=True
    ):
        current = current.sort_values("patient_id")
        for model_spec, extras in config["model_specs"].items():
            seed = stress.stable_seed(
                int(config["wild_bootstrap_seed"]), method, variant, model_spec
            )
            result, bootstrap = fit_response_model(
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
                    "score_variant": variant,
                    "model_spec": model_spec,
                    "extra_covariates": ";".join(extras),
                    **result,
                }
            )
            bootstrap.insert(0, "model_spec", model_spec)
            bootstrap.insert(0, "score_variant", variant)
            bootstrap.insert(0, "scoring_method", method)
            bootstrap_frames.append(bootstrap)
    models = pd.DataFrame(rows)
    models["variant_bh_fdr"] = models.groupby(
        ["scoring_method", "model_spec"], group_keys=False
    )["p_value"].apply(bh_adjust)
    models["wild_variant_bh_fdr"] = models.groupby(
        ["scoring_method", "model_spec"], group_keys=False
    )["wild_bootstrap_p"].apply(bh_adjust)
    models["model_spec_bh_fdr"] = models.groupby(
        ["scoring_method", "score_variant"], group_keys=False
    )["p_value"].apply(bh_adjust)
    models["wild_model_spec_bh_fdr"] = models.groupby(
        ["scoring_method", "score_variant"], group_keys=False
    )["wild_bootstrap_p"].apply(bh_adjust)
    return models, pd.concat(bootstrap_frames, ignore_index=True)


def feature_response_models(
    data: pd.DataFrame,
    feature_column: str,
    outcome_column: str,
    family: str,
    config: dict,
    stress,
) -> pd.DataFrame:
    rows = []
    for feature, current in data.groupby(feature_column, sort=True):
        seed = stress.stable_seed(
            int(config["wild_bootstrap_seed"]), family, str(feature)
        )
        result, _ = fit_response_model(
            current.sort_values("patient_id"),
            outcome_column,
            [],
            int(config["wild_bootstrap_iterations"]),
            seed,
            stress,
        )
        rows.append({feature_column: feature, "model_family": family, **result})
    models = pd.DataFrame(rows)
    models["bh_fdr"] = bh_adjust(models["p_value"])
    models["wild_bh_fdr"] = bh_adjust(models["wild_bootstrap_p"])
    return models


def write_report(
    config: dict,
    pc_variance: pd.DataFrame,
    composition_variance: pd.DataFrame,
    family_models_frame: pd.DataFrame,
    mcp_models: pd.DataFrame,
    pc_models: pd.DataFrame,
    pc_stability: pd.DataFrame,
    pc_qc_correlations: pd.DataFrame,
) -> None:
    def markdown_table(frame: pd.DataFrame) -> str:
        display = frame.copy()
        for column in display.select_dtypes(include=[np.number]).columns:
            display[column] = display[column].map(
                lambda value: "" if pd.isna(value) else f"{float(value):.6g}"
            )
        headers = [str(column) for column in display.columns]
        lines = [
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join(["---"] * len(headers)) + " |",
        ]
        for row in display.itertuples(index=False, name=None):
            lines.append("| " + " | ".join(str(value) for value in row) + " |")
        return "\n".join(lines)

    primary = family_models_frame[
        family_models_frame["score_variant"].eq("module_mean_16")
    ][
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
            "design_condition_number",
        ]
    ]
    lines = [
        "# GSE281729 Global-PC and MCP-counter Composition Stress Test",
        "",
        "## Frozen design",
        "",
        f"- Paired patients: {config['expected_paired_patients']}.",
        "- Global PCs were derived without response labels from paired post-minus-pre expression of non-locked genes.",
        "- MCP-counter v1.1 marker means were used as same-platform relative abundance scores, not absolute cell fractions.",
        f"- Wild bootstrap: {config['wild_bootstrap_iterations']:,} Rademacher draws with HC2 reduced-model residual scaling and an HC3-studentized statistic.",
        "",
        "## Global expression PCs",
        "",
        markdown_table(pc_variance),
        "",
        "## PC-universe stability",
        "",
        markdown_table(pc_stability),
        "",
        "## Global-PC relation to processed-expression distribution",
        "",
        markdown_table(pc_qc_correlations),
        "",
        "## Composition PCs",
        "",
        markdown_table(composition_variance),
        "",
        "## Frozen 16-module family across adjustment models",
        "",
        markdown_table(primary),
        "",
        "## MCP-counter response associations",
        "",
        markdown_table(mcp_models[
            [
                "population",
                "effect",
                "p_value",
                "wild_bootstrap_p",
                "bh_fdr",
                "wild_bh_fdr",
            ]
        ]),
        "",
        "## Global and composition PC response associations",
        "",
        markdown_table(pc_models[
            [
                "pc_source",
                "pc",
                "effect",
                "p_value",
                "wild_bootstrap_p",
                "bh_fdr",
                "wild_bh_fdr",
            ]
        ]),
        "",
        "## Interpretation rule",
        "",
        "The family association may be promoted above a background-sensitive association only if it remains coherent after both global-PC and composition adjustment without unstable design conditioning. All models are reported regardless of direction or significance.",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    config = load_config()
    base = import_script(BASE_SCRIPT, "locked_family")
    stress = import_script(STRESS_SCRIPT, "family_stress")
    expression, annotation = base.load_gse281729()
    paired_expression, paired_annotation, expression_deltas = paired_inputs(
        expression, annotation, int(config["expected_paired_patients"])
    )
    modules = base.load_modules()
    locked_genes = {
        gene
        for genes in modules["genes"]
        for gene in genes
        if gene in paired_expression.index
    }

    finite = paired_expression.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    nonlocked = [gene for gene in finite.index if gene not in locked_genes]
    delta_sds = expression_deltas.loc[nonlocked].std(axis=1, ddof=0)
    means = finite.loc[nonlocked].mean(axis=1)
    primary_genes = [
        gene
        for gene in nonlocked
        if means.loc[gene] >= 1.0 and delta_sds.loc[gene] > 0
    ]
    variable_genes = delta_sds[delta_sds > 0].sort_values(ascending=False).head(5000).index.tolist()
    if len(primary_genes) < 5000 or len(variable_genes) != 5000:
        raise RuntimeError("Global-PC gene universes are unexpectedly small")

    primary_scores, primary_loadings, primary_variance = deterministic_pca(
        expression_deltas, primary_genes, "primary_nonlocked_mean_ge_1", int(config["global_pc_count"])
    )
    variable_scores, variable_loadings, variable_variance = deterministic_pca(
        expression_deltas, variable_genes, "sensitivity_top5000_delta_variance", int(config["global_pc_count"])
    )
    stability_rows = []
    for left_pc in [f"global_pc{i}" for i in range(1, 4)]:
        for right_pc in [f"global_pc{i}" for i in range(1, 4)]:
            stability_rows.append(
                {
                    "primary_pc": left_pc,
                    "sensitivity_pc": right_pc,
                    "pearson_correlation": float(
                        primary_scores[left_pc].corr(variable_scores[right_pc])
                    ),
                    "absolute_correlation": float(
                        abs(primary_scores[left_pc].corr(variable_scores[right_pc]))
                    ),
                }
            )
    pc_stability = pd.DataFrame(stability_rows)
    qc = sample_qc(paired_expression, paired_annotation)
    paired_qc = paired_qc_deltas(qc)
    pc_qc_rows = []
    qc_with_pc = paired_qc.merge(
        primary_scores.drop(columns="gene_universe"),
        on="patient_id",
        how="left",
        validate="one_to_one",
    )
    for pc in [f"global_pc{i}" for i in range(1, 4)]:
        for metric in [column for column in paired_qc.columns if column != "patient_id"]:
            pc_qc_rows.append(
                {
                    "pc": pc,
                    "processed_expression_qc_delta": metric,
                    "pearson_correlation": float(qc_with_pc[pc].corr(qc_with_pc[metric])),
                    "spearman_correlation": float(
                        qc_with_pc[pc].corr(qc_with_pc[metric], method="spearman")
                    ),
                }
            )
    pc_qc_correlations = pd.DataFrame(pc_qc_rows)

    signature_path = ROOT / config["composition_signature_file"]
    mcp_coverage, mcp_samples, mcp_deltas = mcp_counter_scores(
        paired_expression, paired_annotation, signature_path
    )
    composition_scores, composition_loadings, composition_variance = composition_pca(
        mcp_deltas, int(config["composition_pc_count"])
    )

    scores, background, _ = stress.score_components(
        paired_expression, modules, "GSE281729", base
    )
    family_deltas, _ = stress.paired_deltas(
        scores, background, paired_annotation, "GSE281729"
    )
    family_deltas = (
        family_deltas.merge(
            primary_scores.drop(columns="gene_universe"),
            on="patient_id",
            how="left",
            validate="many_to_one",
        )
        .merge(
            composition_scores,
            on="patient_id",
            how="left",
            validate="many_to_one",
        )
        .sort_values(["scoring_method", "score_variant", "patient_id"])
    )
    models, bootstrap = family_models(family_deltas, config, stress)

    mcp_models = feature_response_models(
        mcp_deltas,
        "population",
        "delta_post_minus_pre",
        "mcp_counter_population",
        config,
        stress,
    )
    metadata = paired_annotation[
        [
            "patient_id",
            "response_harmonized_ordinal",
            "response_ord_num",
            "hpv",
            "second_drug",
            "doses",
        ]
    ].drop_duplicates("patient_id")
    pc_long_frames = []
    for source, frame, prefix in [
        ("global", primary_scores.drop(columns="gene_universe"), "global_pc"),
        ("composition", composition_scores, "composition_pc"),
    ]:
        current = frame.merge(metadata, on="patient_id", how="left", validate="one_to_one")
        long = current.melt(
            id_vars=list(metadata.columns),
            value_vars=[f"{prefix}{i}" for i in range(1, 4)],
            var_name="pc",
            value_name="pc_score",
        )
        long["pc_source"] = source
        pc_long_frames.append(long)
    pc_long = pd.concat(pc_long_frames, ignore_index=True)
    pc_model_frames = []
    for source, current in pc_long.groupby("pc_source"):
        fitted = feature_response_models(
            current,
            "pc",
            "pc_score",
            f"{source}_pc",
            config,
            stress,
        )
        fitted.insert(0, "pc_source", source)
        pc_model_frames.append(fitted)
    pc_models = pd.concat(pc_model_frames, ignore_index=True)

    correlation_rows = []
    primary_family = family_deltas[
        family_deltas["score_variant"].eq("module_mean_16")
    ]
    mcp_wide = mcp_deltas.pivot(
        index="patient_id", columns="population", values="delta_post_minus_pre"
    )
    correlates = (
        primary_scores.drop(columns="gene_universe")
        .merge(composition_scores, on="patient_id", validate="one_to_one")
        .merge(mcp_wide.reset_index(), on="patient_id", validate="one_to_one")
        .set_index("patient_id")
    )
    for method, current in primary_family.groupby("scoring_method"):
        family_series = current.set_index("patient_id")["delta_post_minus_pre"]
        for variable in correlates.columns:
            correlation_rows.append(
                {
                    "scoring_method": method,
                    "family_score_variant": "module_mean_16",
                    "correlate": variable,
                    "pearson_correlation": float(
                        family_series.corr(correlates[variable], method="pearson")
                    ),
                    "spearman_correlation": float(
                        family_series.corr(correlates[variable], method="spearman")
                    ),
                }
            )
    correlations = pd.DataFrame(correlation_rows)

    qc.to_csv(SAMPLE_QC_OUT, index=False)
    paired_qc.to_csv(PAIRED_QC_OUT, index=False)
    pc_qc_correlations.to_csv(PC_QC_CORRELATIONS_OUT, index=False)
    pd.concat([primary_scores, variable_scores], ignore_index=True).to_csv(
        PC_SCORES_OUT, index=False
    )
    pd.concat([primary_loadings, variable_loadings], ignore_index=True).to_csv(
        PC_LOADINGS_OUT, index=False
    )
    pd.concat([primary_variance, variable_variance], ignore_index=True).to_csv(
        PC_VARIANCE_OUT, index=False
    )
    pc_stability.to_csv(PC_STABILITY_OUT, index=False)
    mcp_coverage.to_csv(MCP_COVERAGE_OUT, index=False)
    mcp_samples.to_csv(MCP_SAMPLE_OUT, index=False)
    mcp_deltas.to_csv(MCP_DELTA_OUT, index=False)
    composition_scores.to_csv(COMPOSITION_PC_OUT, index=False)
    composition_loadings.to_csv(COMPOSITION_LOADINGS_OUT, index=False)
    family_deltas.to_csv(FAMILY_DELTA_OUT, index=False)
    models.to_csv(FAMILY_MODELS_OUT, index=False)
    bootstrap.to_csv(FAMILY_BOOTSTRAP_OUT, index=False)
    mcp_models.to_csv(MCP_MODELS_OUT, index=False)
    pc_models.to_csv(PC_MODELS_OUT, index=False)
    correlations.to_csv(CORRELATIONS_OUT, index=False)
    write_report(
        config,
        pd.concat([primary_variance, variable_variance], ignore_index=True),
        composition_variance,
        models,
        mcp_models,
        pc_models,
        pc_stability,
        pc_qc_correlations,
    )

    if len(models) != 60 or len(bootstrap) != 60 * int(config["wild_bootstrap_iterations"]):
        raise RuntimeError("Family model or bootstrap output shape is incorrect")
    if mcp_coverage["coverage_fraction"].min() < 0.80:
        raise RuntimeError("MCP-counter marker coverage fell below the frozen 80% gate")
    print(REPORT_OUT)


if __name__ == "__main__":
    main()
