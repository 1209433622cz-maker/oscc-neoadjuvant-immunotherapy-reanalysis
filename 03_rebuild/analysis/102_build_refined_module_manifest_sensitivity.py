#!/usr/bin/env python
"""Build a parallel locked-module manifest from refined discovery lineages."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve()
REBUILD = SCRIPT.parents[1]
REFINED = REBUILD / "results" / "discovery_lineage_composition_sensitivity"
OLD_MANIFEST = REBUILD / "results" / "external_validation" / "GSE123813_gene_set_manifest.csv"
OUT_DIR = REBUILD / "results" / "refined_module_manifest_sensitivity"
OUT_DIR.mkdir(parents=True, exist_ok=True)


T_PATHWAYS = [
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]
M_PATHWAYS = [
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
]


def is_interpretable_symbol(gene: str) -> bool:
    if not isinstance(gene, str) or not gene:
        return False
    if gene.startswith(("AC", "AL", "AP", "LINC", "MIR", "MT-", "RPL", "RPS")):
        return False
    if "." in gene:
        return False
    return bool(re.match(r"^[A-Za-z0-9][A-Za-z0-9.-]*$", gene))


def parse_leading_edge(value: object) -> list[str]:
    if pd.isna(value):
        return []
    genes = re.split(r"[;,]\s*", str(value).strip())
    return [gene for gene in genes if is_interpretable_symbol(gene)]


def top_positive_de(
    frame: pd.DataFrame,
    padj_cutoff: float | None,
    n: int,
    pvalue_cutoff: float | None = None,
) -> list[str]:
    work = frame.copy()
    for column in ["log2FoldChange", "padj", "pvalue"]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work[work["log2FoldChange"] > 0]
    if padj_cutoff is not None:
        work = work[work["padj"] <= padj_cutoff]
    if pvalue_cutoff is not None:
        work = work[work["pvalue"] <= pvalue_cutoff]
    work = work[work["gene"].map(is_interpretable_symbol)]
    work = work.sort_values(["padj", "pvalue", "log2FoldChange"], ascending=[True, True, False])
    return work["gene"].drop_duplicates().head(n).tolist()


def leading_edge_signature(frame: pd.DataFrame, pathway: str) -> list[str]:
    row = frame[frame["pathway"] == pathway]
    if row.empty:
        raise RuntimeError(f"Missing pathway: {pathway}")
    return list(dict.fromkeys(parse_leading_edge(row.iloc[0]["leadingEdge"])))


def build_refined_signatures() -> list[dict[str, object]]:
    t_de = pd.read_csv(REFINED / "T_cell_primary_GENE_RESULTS.csv")
    m_de = pd.read_csv(REFINED / "Myeloid_primary_GENE_RESULTS.csv")
    t_gsea = pd.read_csv(REFINED / "T_cell_primary_HALLMARK_GSEA.csv")
    m_gsea = pd.read_csv(REFINED / "Myeloid_primary_HALLMARK_GSEA.csv")
    signatures: list[dict[str, object]] = [
        {
            "signature": "T_DE_FDR05_positive",
            "target_lineage": "T_cell",
            "source": "Refined OSCC T-cell dynamic DE, positive log2FC, FDR < 0.05",
            "genes": top_positive_de(t_de, 0.05, 50),
        },
        {
            "signature": "M_DE_FDR10_positive",
            "target_lineage": "Myeloid",
            "source": "Refined OSCC myeloid dynamic DE, positive log2FC, FDR < 0.10",
            "genes": top_positive_de(m_de, 0.10, 50),
        },
        {
            "signature": "M_DE_nominal_top30_positive",
            "target_lineage": "Myeloid",
            "source": "Refined OSCC myeloid dynamic DE, positive log2FC, nominal P < 0.01, top 30",
            "genes": top_positive_de(m_de, None, 30, 0.01),
        },
    ]
    t_union: list[str] = []
    for pathway in T_PATHWAYS:
        genes = leading_edge_signature(t_gsea, pathway)
        t_union.extend(genes)
        signatures.append(
            {
                "signature": pathway.replace("HALLMARK_", "T_LE_"),
                "target_lineage": "T_cell",
                "source": f"Refined OSCC T-cell Hallmark leading edge: {pathway}",
                "genes": genes,
            }
        )
    signatures.append(
        {
            "signature": "T_LE_union_core",
            "target_lineage": "T_cell",
            "source": "Union of selected refined OSCC T-cell Hallmark leading-edge genes",
            "genes": list(dict.fromkeys(t_union)),
        }
    )
    m_union: list[str] = []
    for pathway in M_PATHWAYS:
        genes = leading_edge_signature(m_gsea, pathway)
        m_union.extend(genes)
        signatures.append(
            {
                "signature": pathway.replace("HALLMARK_", "M_LE_"),
                "target_lineage": "Myeloid",
                "source": f"Refined OSCC myeloid Hallmark leading edge: {pathway}",
                "genes": genes,
            }
        )
    signatures.append(
        {
            "signature": "M_LE_union_core",
            "target_lineage": "Myeloid",
            "source": "Union of selected refined OSCC myeloid Hallmark leading-edge genes",
            "genes": list(dict.fromkeys(m_union)),
        }
    )
    return signatures


def parse_genes(value: object) -> set[str]:
    if pd.isna(value) or not str(value).strip():
        return set()
    return {gene for gene in str(value).split(";") if gene}


def main() -> None:
    signatures = build_refined_signatures()
    if len(signatures) != 16:
        raise RuntimeError(f"Expected 16 refined modules, found {len(signatures)}")
    refined = pd.DataFrame(
        [
            {
                "signature": item["signature"],
                "target_lineage": item["target_lineage"],
                "source": item["source"],
                "n_genes_defined": len(item["genes"]),
                "genes_defined": ";".join(item["genes"]),
            }
            for item in signatures
        ]
    )
    refined.to_csv(OUT_DIR / "REFINED_MODULE_MANIFEST.csv", index=False)

    old = pd.read_csv(OLD_MANIFEST)
    merged = old[["signature", "target_lineage", "n_genes_defined", "genes_defined"]].merge(
        refined[["signature", "n_genes_defined", "genes_defined"]],
        on="signature",
        suffixes=("", ""),
        validate="one_to_one",
    )
    rows: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        old_genes = parse_genes(row["genes_defined"])
        new_genes = parse_genes(row["genes_defined"])
        intersection = old_genes & new_genes
        union = old_genes | new_genes
        rows.append(
            {
                "signature": row["signature"],
                "target_lineage": row["target_lineage"],
                "n_genes": len(old_genes),
                "n_genes": len(new_genes),
                "n_intersection": len(intersection),
                "n_union": len(union),
                "jaccard": len(intersection) / len(union) if union else 1.0,
                "v1_covered_by": len(intersection) / len(old_genes) if old_genes else 1.0,
                "v2_covered_by": len(intersection) / len(new_genes) if new_genes else 1.0,
                "genes_only": ";".join(sorted(old_genes - new_genes)),
                "genes_only": ";".join(sorted(new_genes - old_genes)),
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(OUT_DIR / "MODULE_MANIFEST_COMPARISON.csv", index=False)

    old_family = set().union(*(parse_genes(value) for value in old["genes_defined"]))
    new_family = set().union(*(parse_genes(value) for value in refined["genes_defined"]))
    family_intersection = old_family & new_family
    family_union = old_family | new_family
    family = pd.DataFrame(
        [
            {
                "n_unique_genes": len(old_family),
                "n_unique_genes": len(new_family),
                "n_intersection": len(family_intersection),
                "n_union": len(family_union),
                "jaccard": len(family_intersection) / len(family_union),
                "v1_covered_by": len(family_intersection) / len(old_family),
                "v2_covered_by": len(family_intersection) / len(new_family),
            }
        ]
    )
    family.to_csv(OUT_DIR / "FAMILY_GENE_OVERLAP_SUMMARY.csv", index=False)

    report = [
        "# Refined Module Manifest Sensitivity",
        "",
        f"- Old manifest: `{OLD_MANIFEST}`",
        f"- Refined manifest: `{OUT_DIR / 'REFINED_MODULE_MANIFEST.csv'}`",
        "- The module recipe and pathway panel are unchanged; only discovery DE/GSEA inputs differ.",
        "- This is a post-freeze sensitivity manifest and does not overwrite the original locked family.",
        "",
        "## Family overlap",
        "",
        family.to_string(index=False),
        "",
        "## Module overlap",
        "",
        comparison[
            [
                "signature",
                "target_lineage",
                "n_genes",
                "n_genes",
                "n_intersection",
                "jaccard",
            ]
        ].to_string(index=False),
        "",
        "A full external rerun is required if the refined family is used as the manuscript's locked validation family. "
        "The original family remains the only genuinely pre-sensitivity frozen family.",
    ]
    (OUT_DIR / "REFINED_MODULE_MANIFEST_SENSITIVITY_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(OUT_DIR / "REFINED_MODULE_MANIFEST_SENSITIVITY_REPORT.md")


if __name__ == "__main__":
    main()
