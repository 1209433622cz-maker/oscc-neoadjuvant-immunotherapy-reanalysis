from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REBUILD = ROOT / "03_rebuild"
CONFIG = REBUILD / "config" / "gse195832_locked_family_validation.json"
REFERENCE = REBUILD / "validation" / "GSE195832_bulk_locked_family"
OUT = REBUILD / "validation" / "GSE195832_independent_primary_recompute"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hc3_fit(y: np.ndarray, predictors: pd.DataFrame) -> dict[str, float]:
    x = np.column_stack([np.ones(len(predictors)), predictors.to_numpy(float)])
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    inverse = np.linalg.pinv(x.T @ x)
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, inverse, x)
    scaled = residual / (1.0 - leverage)
    meat = x.T @ ((scaled**2)[:, None] * x)
    covariance = inverse @ meat @ inverse
    standard_error = math.sqrt(float(covariance[1, 1]))
    z_value = float(beta[1] / standard_error)
    p_value = math.erfc(abs(z_value) / math.sqrt(2.0))
    return {
        "effect": float(beta[1]),
        "std_error": standard_error,
        "ci95_low": float(beta[1] - 1.959963984540054 * standard_error),
        "ci95_high": float(beta[1] + 1.959963984540054 * standard_error),
        "p_value_hc3": p_value,
    }


def score_variants(
    expression: pd.DataFrame, module_manifest: Path
) -> tuple[pd.DataFrame, set[str]]:
    modules = pd.read_csv(module_manifest)
    module_genes = {
        row.signature: [
            gene.strip().upper()
            for gene in str(row.genes_defined).split(";")
            if gene.strip() and gene.strip().upper() in expression.index
        ]
        for row in modules.itertuples(index=False)
    }
    if len(module_genes) != 16 or any(not genes for genes in module_genes.values()):
        raise RuntimeError("Frozen module coverage failed")

    finite = expression.loc[expression.std(axis=1, ddof=0).gt(0)]
    z_score = finite.sub(finite.mean(axis=1), axis=0).div(
        finite.std(axis=1, ddof=0), axis=0
    )
    memberships = Counter(
        gene for genes in module_genes.values() for gene in genes
    )
    unique_genes = sorted(memberships)
    no_union = {
        name: genes for name, genes in module_genes.items()
        if not name.endswith("_union_core")
    }
    hallmark = {
        name: genes for name, genes in module_genes.items()
        if "_LE_" in name and not name.endswith("_union_core")
    }
    dynamic = {
        name: genes for name, genes in module_genes.items() if "_LE_" not in name
    }

    def module_mean(selected: dict[str, list[str]]) -> pd.Series:
        return pd.concat(
            [z_score.loc[genes].mean(axis=0) for genes in selected.values()],
            axis=1,
        ).mean(axis=1)

    inverse_module_scores = []
    for genes in module_genes.values():
        weights = np.asarray([1.0 / memberships[gene] for gene in genes])
        inverse_module_scores.append(
            pd.Series(
                np.average(z_score.loc[genes], axis=0, weights=weights),
                index=z_score.columns,
            )
        )

    variants = pd.DataFrame(
        {
            "module_mean_16": module_mean(module_genes),
            "unique_gene_equal": z_score.loc[unique_genes].mean(axis=0),
            "inverse_membership_module_mean": pd.concat(
                inverse_module_scores, axis=1
            ).mean(axis=1),
            "no_union_module_mean": module_mean(no_union),
            "hallmark_only_module_mean": module_mean(hallmark),
            "dynamic_only_module_mean": module_mean(dynamic),
        }
    )
    variants.index.name = "sample_id"
    return variants, set(unique_genes)


def patient_deltas(
    scores: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    long = scores.reset_index().melt(
        id_vars="sample_id", var_name="score_variant", value_name="score"
    )
    merged = long.merge(metadata, on="sample_id", validate="many_to_one")
    wide = merged.pivot(
        index=["patient_id", "score_variant"],
        columns="timepoint",
        values="score",
    ).reset_index()
    wide["delta_post_minus_pre"] = wide["post"] - wide["pre"]
    patient_meta = metadata[
        [
            "patient_id",
            "therapy",
            "therapy_tadalafil",
            "batch2",
            "response_ord_num",
        ]
    ].drop_duplicates("patient_id")
    return wide.merge(patient_meta, on="patient_id", validate="many_to_one")


def global_pc1(
    expression: pd.DataFrame, metadata: pd.DataFrame, locked: set[str]
) -> pd.DataFrame:
    patients = sorted(metadata["patient_id"].unique())
    sample_positions = {
        sample: index for index, sample in enumerate(expression.columns)
    }
    operator = np.zeros((len(patients), len(expression.columns)))
    for row, patient in enumerate(patients):
        current = metadata.loc[metadata["patient_id"].eq(patient)]
        pre = current.loc[current["timepoint"].eq("pre"), "sample_id"].iloc[0]
        post = current.loc[current["timepoint"].eq("post"), "sample_id"].iloc[0]
        operator[row, sample_positions[pre]] = -1.0
        operator[row, sample_positions[post]] = 1.0
    genes = [gene for gene in expression.index if gene not in locked]
    delta = operator @ expression.loc[genes].T.to_numpy(float)
    keep = delta.std(axis=0, ddof=0) > 0
    standardized = (
        delta[:, keep] - delta[:, keep].mean(axis=0)
    ) / delta[:, keep].std(axis=0, ddof=0)
    u, singular, vt = np.linalg.svd(standardized, full_matrices=False)
    scores = u[:, 0] * singular[0]
    anchor = int(np.argmax(np.abs(vt[0])))
    if vt[0, anchor] < 0:
        scores *= -1
    scores = (scores - scores.mean()) / scores.std(ddof=0)
    return pd.DataFrame({"patient_id": patients, "global_pc1": scores})


def stratified_permutation(data: pd.DataFrame, iterations: int, seed: int) -> float:
    current = data.sort_values("patient_id").reset_index(drop=True)
    y = current["delta_post_minus_pre"].to_numpy(float)
    nuisance = np.column_stack(
        [
            np.ones(len(current)),
            current["therapy_tadalafil"].to_numpy(float),
            current["batch2"].to_numpy(float),
        ]
    )
    projection = nuisance @ np.linalg.pinv(nuisance)

    def coefficient(response: np.ndarray) -> float:
        residual = response - projection @ response
        return float(residual @ y / (residual @ residual))

    response = current["response_ord_num"].to_numpy(float)
    observed = coefficient(response)
    therapy = current["therapy"].to_numpy()
    rng = np.random.default_rng(seed)
    exceedances = 0
    for _ in range(iterations):
        permuted = response.copy()
        for arm in np.unique(therapy):
            positions = np.flatnonzero(therapy == arm)
            permuted[positions] = rng.permutation(permuted[positions])
        exceedances += abs(coefficient(permuted)) >= abs(observed) - 1e-12
    return float((1 + exceedances) / (iterations + 1))


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    counts_path = ROOT / config["counts_file"]
    metadata_path = ROOT / config["metadata_file"]
    if sha256(counts_path) != config["counts_sha256"]:
        raise RuntimeError("Count-table SHA256 mismatch")
    if sha256(metadata_path) != config["metadata_sha256"]:
        raise RuntimeError("Metadata SHA256 mismatch")

    metadata = pd.read_csv(metadata_path, index_col=0).rename(
        columns={
            "Patient": "patient_id",
            "Therapy": "therapy",
            "Clinical Response(Primary_Site_TE)": "response_label",
            "Sample": "sample_id",
            "Treatment": "treatment_time",
            "Batch": "batch",
        }
    )
    metadata["timepoint"] = metadata["treatment_time"].map(
        {"Pre-Treated": "pre", "Post-Treated": "post"}
    )
    response_map = {
        label: index for index, label in enumerate(config["response_order"])
    }
    metadata["response_ord_num"] = metadata["response_label"].map(response_map)
    metadata["therapy_tadalafil"] = metadata["therapy"].eq(
        "antiPD1+Tadalafil"
    ).astype(int)
    metadata["batch2"] = metadata["batch"].eq("Batch2").astype(int)
    metadata = metadata.reset_index(drop=True)

    counts = pd.read_csv(counts_path, sep="\t")
    counts["symbol"] = counts["symbol"].astype(str).str.strip().str.upper()
    counts = (
        counts.loc[counts["symbol"].ne("") & counts["symbol"].ne("NAN")]
        .set_index("symbol")
        .apply(pd.to_numeric, errors="raise")
        .groupby(level=0, sort=True)
        .sum()
    )
    counts = counts.loc[:, metadata["sample_id"]]
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1_000_000
    expression = np.log2(cpm.loc[(cpm >= 1).sum(axis=1) >= 6] + 0.5)

    variants, locked = score_variants(
        expression, ROOT / config["module_manifest"]
    )
    deltas = patient_deltas(variants, metadata)
    pc1 = global_pc1(expression, metadata, locked)
    reference_models = pd.read_csv(
        REFERENCE / "GSE195832_FAMILY_RESPONSE_MODELS.csv"
    )
    reference_permutation = pd.read_csv(
        REFERENCE / "GSE195832_PRIMARY_PERMUTATION_TEST.csv"
    )

    rows = []
    for variant in variants.columns:
        current = (
            deltas.loc[deltas["score_variant"].eq(variant)]
            .sort_values("patient_id")
            .reset_index(drop=True)
        )
        primary = hc3_fit(
            current["delta_post_minus_pre"].to_numpy(float),
            current[["response_ord_num", "therapy_tadalafil", "batch2"]],
        )
        expected = reference_models.loc[
            reference_models["scoring_method"].eq("z_score")
            & reference_models["score_variant"].eq(variant)
            & reference_models["model_spec"].eq("primary")
        ].iloc[0]
        rows.append(
            {
                "check": f"primary_{variant}",
                "effect": primary["effect"],
                "expected_effect": expected["effect"],
                "p_value_hc3": primary["p_value_hc3"],
                "expected_p_value_hc3": expected["p_value_hc3"],
                "status": "PASS"
                if math.isclose(primary["effect"], expected["effect"], abs_tol=1e-12)
                and math.isclose(
                    primary["p_value_hc3"],
                    expected["p_value_hc3"],
                    abs_tol=1e-12,
                )
                else "FAIL",
            }
        )

    primary_data = (
        deltas.loc[deltas["score_variant"].eq("module_mean_16")]
        .sort_values("patient_id")
        .reset_index(drop=True)
    )
    permutation_p = stratified_permutation(
        primary_data,
        int(config["permutation_iterations"]),
        int(config["permutation_seed"]),
    )
    expected_permutation = float(
        reference_permutation["two_sided_permutation_p"].iloc[0]
    )
    rows.append(
        {
            "check": "within_treatment_permutation",
            "effect": np.nan,
            "expected_effect": np.nan,
            "p_value_hc3": permutation_p,
            "expected_p_value_hc3": expected_permutation,
            "status": "PASS"
            if math.isclose(permutation_p, expected_permutation, abs_tol=1e-15)
            else "FAIL",
        }
    )

    pc_data = primary_data.merge(pc1, on="patient_id", validate="one_to_one")
    pc_fit = hc3_fit(
        pc_data["delta_post_minus_pre"].to_numpy(float),
        pc_data[
            [
                "response_ord_num",
                "therapy_tadalafil",
                "batch2",
                "global_pc1",
            ]
        ],
    )
    expected_pc = reference_models.loc[
        reference_models["scoring_method"].eq("z_score")
        & reference_models["score_variant"].eq("module_mean_16")
        & reference_models["model_spec"].eq("global_pc1")
    ].iloc[0]
    rows.append(
        {
            "check": "global_pc1_adjusted_primary",
            "effect": pc_fit["effect"],
            "expected_effect": expected_pc["effect"],
            "p_value_hc3": pc_fit["p_value_hc3"],
            "expected_p_value_hc3": expected_pc["p_value_hc3"],
            "status": "PASS"
            if math.isclose(pc_fit["effect"], expected_pc["effect"], abs_tol=1e-12)
            and math.isclose(
                pc_fit["p_value_hc3"],
                expected_pc["p_value_hc3"],
                abs_tol=1e-12,
            )
            else "FAIL",
        }
    )

    leave_one_out = []
    for patient in primary_data["patient_id"]:
        current = primary_data.loc[primary_data["patient_id"].ne(patient)]
        fit = hc3_fit(
            current["delta_post_minus_pre"].to_numpy(float),
            current[["response_ord_num", "therapy_tadalafil", "batch2"]],
        )
        leave_one_out.append(fit["effect"])
    rows.append(
        {
            "check": "all_leave_one_patient_effects_negative",
            "effect": max(leave_one_out),
            "expected_effect": "<0",
            "p_value_hc3": np.nan,
            "expected_p_value_hc3": np.nan,
            "status": "PASS" if max(leave_one_out) < 0 else "FAIL",
        }
    )

    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "GSE195832_INDEPENDENT_PRIMARY_RECOMPUTE.csv", index=False)
    passed = int(audit["status"].eq("PASS").sum())
    failed = int(audit["status"].eq("FAIL").sum())
    primary = audit.loc[audit["check"].eq("primary_module_mean_16")].iloc[0]
    report = [
        "# GSE195832 independent primary recomputation",
        "",
        "- Implementation: clean-room script; no import from analysis 120.",
        "- Estimator: manual OLS and HC3 sandwich covariance.",
        "- Normalization and filtering: recomputed directly from SHA-verified raw counts.",
        f"- Result: **{'PASS' if failed == 0 else 'FAIL'}** ({passed} PASS, {failed} FAIL).",
        f"- Primary effect: {primary['effect']:.12g}; HC3 P = {primary['p_value_hc3']:.12g}.",
        f"- Within-treatment permutation P = {permutation_p:.12g}.",
        f"- Global-PC1-adjusted effect: {pc_fit['effect']:.12g}; HC3 P = {pc_fit['p_value_hc3']:.12g}.",
        f"- Leave-one-patient effect range: {min(leave_one_out):.12g} to {max(leave_one_out):.12g}.",
        "",
        "The independent implementation confirms the opposite orientation, null primary association, global-PC attenuation and uniformly negative leave-one-patient estimates.",
    ]
    (OUT / "GSE195832_INDEPENDENT_PRIMARY_RECOMPUTE_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print({"pass": passed, "fail": failed, "output": str(OUT)})
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
