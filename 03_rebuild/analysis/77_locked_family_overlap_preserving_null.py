from __future__ import annotations

import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REBUILD = ROOT / "03_rebuild"
BASE_SCRIPT = REBUILD / "analysis" / "61_locked_family_robustness_external_cohorts.py"
CONFIG_PATH = REBUILD / "config" / "locked_family_overlap_null.json"
OUT_DIR = REBUILD / "validation" / "locked_family_overlap_null"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_OUT = OUT_DIR / "OVERLAP_PRESERVING_MATCHED_RANDOM_EFFECTS.csv"
TEST_OUT = OUT_DIR / "OVERLAP_PRESERVING_FAMILY_TESTS.csv"
AUDIT_OUT = OUT_DIR / "OVERLAP_PRESERVING_NULL_AUDIT.csv"
REPORT_OUT = OUT_DIR / "OVERLAP_PRESERVING_NULL_REPORT.md"


def load_base_module():
    spec = importlib.util.spec_from_file_location("locked_family", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_config() -> dict:
    if not CONFIG_PATH.exists():
        raise FileNotFoundError(CONFIG_PATH)
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def membership_template(
    expression: pd.DataFrame,
    modules: pd.DataFrame,
    base,
) -> tuple[
    dict[str, tuple[str, ...]],
    dict[int, np.ndarray],
    dict[int, list[str]],
    pd.Series,
]:
    gene_index = expression.index
    bins = base.expression_bins(expression)
    membership: dict[str, list[str]] = defaultdict(list)
    for _, module in modules.iterrows():
        signature = str(module["signature"])
        for gene in module["genes"]:
            if gene in gene_index:
                membership[gene].append(signature)

    locked_present = set(membership)
    background_by_bin: dict[int, np.ndarray] = {}
    for bin_id in sorted(bins.unique()):
        positions = np.asarray(
            [
                position
                for position, gene in enumerate(gene_index)
                if gene not in locked_present
                and int(bins.iloc[position]) == int(bin_id)
            ],
            dtype=int,
        )
        background_by_bin[int(bin_id)] = positions

    template_by_bin: dict[int, list[str]] = defaultdict(list)
    for gene in sorted(membership):
        template_by_bin[int(bins.loc[gene])].append(gene)
    for bin_id, genes in template_by_bin.items():
        if len(background_by_bin[bin_id]) < len(genes):
            raise RuntimeError(
                f"Insufficient background genes in expression bin {bin_id}: "
                f"need {len(genes)}, have {len(background_by_bin[bin_id])}."
            )

    frozen_membership = {
        gene: tuple(sorted(signatures))
        for gene, signatures in membership.items()
    }
    return frozen_membership, background_by_bin, dict(template_by_bin), bins


def draw_overlap_preserving_modules(
    rng: np.random.Generator,
    membership: dict[str, tuple[str, ...]],
    background_by_bin: dict[int, np.ndarray],
    template_by_bin: dict[int, list[str]],
) -> dict[str, np.ndarray]:
    selected_by_signature: dict[str, list[int]] = defaultdict(list)
    for bin_id, original_genes in template_by_bin.items():
        replacements = rng.choice(
            background_by_bin[bin_id],
            size=len(original_genes),
            replace=False,
        )
        for original_gene, replacement in zip(original_genes, replacements):
            for signature in membership[original_gene]:
                selected_by_signature[signature].append(int(replacement))
    return {
        signature: np.asarray(positions, dtype=int)
        for signature, positions in selected_by_signature.items()
    }


def effect_weights(
    cohort: str,
    methods: list[str],
    expression: pd.DataFrame,
    annotation: pd.DataFrame,
    observed_deltas: pd.DataFrame,
    gse281_weights: dict[str, np.ndarray] | None,
    base,
) -> dict[str, np.ndarray]:
    patient_order = sorted(observed_deltas["patient_id"].unique())
    patient_to_sample_delta = base.delta_operator(
        expression.columns.tolist(),
        annotation,
        patient_order,
    )
    if cohort == "GSE281729":
        if gse281_weights is None:
            raise RuntimeError("Missing GSE281729 effect weights.")
        return {
            method: gse281_weights[method] @ patient_to_sample_delta
            for method in methods
        }

    patient_meta = (
        observed_deltas.drop_duplicates("patient_id")
        .set_index("patient_id")
        .loc[patient_order]
    )
    labels = patient_meta["response_binary"].to_numpy() == "responder"
    patient_effect = np.where(
        labels,
        1.0 / labels.sum(),
        -1.0 / (~labels).sum(),
    )
    return {
        method: patient_effect @ patient_to_sample_delta
        for method in methods
    }


def membership_audit(
    cohort: str,
    membership: dict[str, tuple[str, ...]],
    indices: dict[str, np.ndarray],
) -> dict[str, object]:
    observed_sizes = Counter(
        signature
        for signatures in membership.values()
        for signature in signatures
    )
    random_sizes = {signature: len(values) for signature, values in indices.items()}
    multiplicities = Counter(len(signatures) for signatures in membership.values())
    return {
        "cohort": cohort,
        "unique_locked_genes_present": len(membership),
        "total_module_memberships": sum(observed_sizes.values()),
        "genes_reused_across_modules": sum(
            len(signatures) > 1 for signatures in membership.values()
        ),
        "maximum_module_membership_per_gene": max(map(len, membership.values())),
        "membership_multiplicity_profile": ";".join(
            f"{multiplicity}:{count}"
            for multiplicity, count in sorted(multiplicities.items())
        ),
        "module_sizes_exactly_preserved": all(
            random_sizes.get(signature, 0) == size
            for signature, size in observed_sizes.items()
        ),
        "unique_replacements_preserved": (
            len(set(np.concatenate(list(indices.values()))))
            == len(membership)
        ),
        "null_definition": (
            "one expression-decile-matched unique replacement per measured "
            "locked gene, inherited by every module containing that gene"
        ),
    }


def overlap_random_controls(
    config: dict,
    cohort: str,
    expression: pd.DataFrame,
    annotation: pd.DataFrame,
    modules: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    observed_deltas: pd.DataFrame,
    gse281_weights: dict[str, np.ndarray] | None,
    base,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(int(config["random_seed"]))
    (
        membership,
        background_by_bin,
        template_by_bin,
        _,
    ) = membership_template(expression, modules, base)
    weights = effect_weights(
        cohort,
        list(arrays),
        expression,
        annotation,
        observed_deltas,
        gse281_weights,
        base,
    )

    rows: list[dict[str, object]] = []
    first_indices: dict[str, np.ndarray] | None = None
    for iteration in range(1, int(config["random_control_iterations"]) + 1):
        indices = draw_overlap_preserving_modules(
            rng,
            membership,
            background_by_bin,
            template_by_bin,
        )
        if first_indices is None:
            first_indices = indices
        sample_scores = base.random_family_sample_score(arrays, indices)
        for method, values in sample_scores.items():
            rows.append(
                {
                    "cohort": cohort,
                    "scoring_method": method,
                    "iteration": iteration,
                    "random_seed": int(config["random_seed"]),
                    "null_type": "overlap_preserving_expression_decile",
                    "effect": float(weights[method] @ values),
                }
            )
    if first_indices is None:
        raise RuntimeError("No random families were generated.")
    return pd.DataFrame(rows), membership_audit(
        cohort,
        membership,
        first_indices,
    )


def attach_specificity(
    tests: pd.DataFrame,
    random_effects: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, test in tests.iterrows():
        null = random_effects[
            (random_effects["cohort"] == test["cohort"])
            & (random_effects["scoring_method"] == test["scoring_method"])
        ]["effect"].to_numpy(float)
        observed = float(test["effect"])
        empirical_p = float(
            (1 + np.sum(np.abs(null) >= abs(observed) - 1e-12))
            / (len(null) + 1)
        )
        row = test.to_dict()
        row.update(
            {
                "null_type": "overlap_preserving_expression_decile",
                "matched_random_iterations": len(null),
                "matched_random_mean": float(null.mean()),
                "matched_random_sd": float(null.std(ddof=1)),
                "matched_random_q025": float(np.quantile(null, 0.025)),
                "matched_random_q975": float(np.quantile(null, 0.975)),
                "empirical_specificity_p": empirical_p,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def markdown_table(frame: pd.DataFrame) -> str:
    def render(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.6g}"
        return str(value).replace("|", r"\|").replace("\n", " ")

    columns = [str(column) for column in frame.columns]
    rows = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    rows.extend(
        "| " + " | ".join(render(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(rows)


def write_report(
    config: dict,
    tests: pd.DataFrame,
    audit: pd.DataFrame,
) -> None:
    lines = [
        "# Locked-Family Overlap-Preserving Null",
        "",
        "## Null design",
        "",
        (
            "Each measured locked gene is replaced by one unique background gene "
            "from the same mean-expression decile. The replacement inherits every "
            "module membership of the original gene, exactly preserving module "
            "sizes, pairwise and higher-order overlaps, and the gene-reuse profile."
        ),
        "",
        f"- Iterations per cohort/method: {config['random_control_iterations']}",
        f"- Random seed: {config['random_seed']}",
        "",
        "## Membership audit",
        "",
        markdown_table(audit),
        "",
        "## Family tests",
        "",
        markdown_table(
            tests[
                [
                    "cohort",
                    "scoring_method",
                    "effect",
                    "p_value",
                    "matched_random_q025",
                    "matched_random_q975",
                    "empirical_specificity_p",
                ]
            ]
        ),
        "",
        "## Interpretation",
        "",
        (
            "This is a stricter specificity sensitivity analysis than independent "
            "module-wise randomization. It tests whether the observed family effect "
            "is extreme after preserving the complete measured gene-to-module "
            "incidence structure as well as expression decile."
        ),
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    config = load_config()
    base = load_base_module()
    modules = base.load_modules()

    expr179, annotation179 = base.load_gse179730()
    scores179, _, _, arrays179, eligible179 = base.score_family(
        expr179, modules, "GSE179730"
    )
    deltas179 = base.paired_deltas(scores179, annotation179, "GSE179730")
    tests179, _ = base.exact_gse179730(deltas179)

    expr281, annotation281 = base.load_gse281729()
    scores281, _, _, arrays281, eligible281 = base.score_family(
        expr281, modules, "GSE281729"
    )
    deltas281 = base.paired_deltas(scores281, annotation281, "GSE281729")
    deltas281 = deltas281[
        deltas281["response_harmonized_ordinal"].isin(["Low", "Medium", "High"])
    ].copy()
    tests281, model_weights281 = base.model_gse281729(deltas281)

    random179, audit179 = overlap_random_controls(
        config,
        "GSE179730",
        eligible179,
        annotation179,
        modules,
        arrays179,
        deltas179,
        None,
        base,
    )
    random281, audit281 = overlap_random_controls(
        config,
        "GSE281729",
        eligible281,
        annotation281,
        modules,
        arrays281,
        deltas281,
        model_weights281,
        base,
    )

    random_effects = pd.concat([random179, random281], ignore_index=True)
    tests = attach_specificity(
        pd.DataFrame(tests179 + tests281),
        random_effects,
    )
    audit = pd.DataFrame([audit179, audit281])
    random_effects.to_csv(RANDOM_OUT, index=False)
    tests.to_csv(TEST_OUT, index=False)
    audit.to_csv(AUDIT_OUT, index=False)
    write_report(config, tests, audit)
    print(TEST_OUT)
    print(REPORT_OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
