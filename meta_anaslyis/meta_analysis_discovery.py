#!/usr/bin/env python3
"""
Fixed-effect inverse-variance-weighted meta-analysis of UK Biobank (UKB)
and All of Us (AoU) rare-variant burden results.

The script expects one row per gene in each cohort file with these columns:
    Gene
    Cases_with_variant
    Cases_without_variant
    Controls_with_variant
    Controls_without_variant
    Odds_Ratio
    95%_CI_Lower
    95%_CI_Upper
    P_value

For each cohort, beta is reconstructed as log(OR) and the standard error is
reconstructed from the reported 95% confidence interval:

    SE = [log(CI_upper) - log(CI_lower)] / (2 * 1.96)

Only genes with valid results in both UKB and AoU are included. The script
performs a two-cohort fixed-effect IVW meta-analysis and reports Cochran's Q
and the corresponding heterogeneity P value.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import mpmath as mp
import numpy as np
import pandas as pd
from scipy.stats import chi2


Z975 = 1.96
mp.prec = 200

REQUIRED_COLUMNS = {
    "Gene",
    "Cases_with_variant",
    "Cases_without_variant",
    "Controls_with_variant",
    "Controls_without_variant",
    "Odds_Ratio",
    "95%_CI_Lower",
    "95%_CI_Upper",
    "P_value",
}

NUMERIC_COLUMNS = [
    "Cases_with_variant",
    "Cases_without_variant",
    "Controls_with_variant",
    "Controls_without_variant",
    "Odds_Ratio",
    "95%_CI_Lower",
    "95%_CI_Upper",
    "P_value",
]


def compute_p_value_mpmath(z: float) -> mp.mpf:
    """Return a two-sided normal-theory P value with high precision."""
    z_mp = mp.mpf(z)
    return mp.erfc(mp.fabs(z_mp) / mp.sqrt(2))


def format_p_value(value: object) -> str:
    """Format P values in scientific notation while retaining very small values."""
    try:
        if isinstance(value, mp.mpf):
            return mp.nstr(value, n=3, min_fixed=0, max_fixed=0)
        return f"{float(value):.3e}"
    except (TypeError, ValueError):
        return str(value)


def _validate_unique_genes(df: pd.DataFrame, label: str) -> None:
    duplicated = df.loc[df["Gene"].duplicated(keep=False), "Gene"].astype(str).unique()
    if len(duplicated) > 0:
        preview = ", ".join(duplicated[:10])
        raise ValueError(
            f"{label}: expected one row per gene, but duplicate Gene values were found "
            f"(examples: {preview})."
        )


def read_cohort_results(path: Path, label: str) -> pd.DataFrame:
    """Read, validate, and normalize one cohort-level Firth result file."""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(
            f"{path}: missing required column(s): {', '.join(sorted(missing))}"
        )

    df["Gene"] = df["Gene"].astype(str).str.strip()
    _validate_unique_genes(df, label)

    for column in NUMERIC_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    # A valid meta-analysis row requires a positive OR and positive CI bounds.
    valid = (
        np.isfinite(df["Odds_Ratio"])
        & np.isfinite(df["95%_CI_Lower"])
        & np.isfinite(df["95%_CI_Upper"])
        & (df["Odds_Ratio"] > 0)
        & (df["95%_CI_Lower"] > 0)
        & (df["95%_CI_Upper"] > 0)
        & df["P_value"].between(0, 1, inclusive="both")
    )
    df = df.loc[valid].copy()

    df["logOR"] = np.log(df["Odds_Ratio"])
    df["SE"] = (
        np.log(df["95%_CI_Upper"]) - np.log(df["95%_CI_Lower"])
    ) / (2 * Z975)

    df = df[np.isfinite(df["SE"]) & (df["SE"] > 0)].copy()
    return df


def meta_analysis(ukb: pd.DataFrame, aou: pd.DataFrame) -> pd.DataFrame:
    """Perform two-cohort fixed-effect IVW meta-analysis by gene."""
    merged = pd.merge(ukb, aou, on="Gene", suffixes=("_UKB", "_AOU"), how="inner")

    if merged.empty:
        return pd.DataFrame()

    merged["weight_UKB"] = 1.0 / np.square(merged["SE_UKB"])
    merged["weight_AOU"] = 1.0 / np.square(merged["SE_AOU"])
    weight_sum = merged["weight_UKB"] + merged["weight_AOU"]

    merged["meta_logOR"] = (
        merged["logOR_UKB"] * merged["weight_UKB"]
        + merged["logOR_AOU"] * merged["weight_AOU"]
    ) / weight_sum
    merged["meta_SE"] = np.sqrt(1.0 / weight_sum)
    merged["meta_OR"] = np.exp(merged["meta_logOR"])
    merged["meta_CI_lower"] = np.exp(merged["meta_logOR"] - Z975 * merged["meta_SE"])
    merged["meta_CI_upper"] = np.exp(merged["meta_logOR"] + Z975 * merged["meta_SE"])
    merged["meta_Z"] = merged["meta_logOR"] / merged["meta_SE"]
    merged["meta_P_value"] = merged["meta_Z"].apply(compute_p_value_mpmath)

    merged["Q"] = (
        merged["weight_UKB"] * np.square(merged["logOR_UKB"] - merged["meta_logOR"])
        + merged["weight_AOU"] * np.square(merged["logOR_AOU"] - merged["meta_logOR"])
    )
    merged["p_het"] = 1.0 - chi2.cdf(merged["Q"], df=1)

    merged["UKB_case_carriers"] = merged["Cases_with_variant_UKB"]
    merged["UKB_case_noncarriers"] = merged["Cases_without_variant_UKB"]
    merged["UKB_control_carriers"] = merged["Controls_with_variant_UKB"]
    merged["UKB_control_noncarriers"] = merged["Controls_without_variant_UKB"]
    merged["AOU_case_carriers"] = merged["Cases_with_variant_AOU"]
    merged["AOU_case_noncarriers"] = merged["Cases_without_variant_AOU"]
    merged["AOU_control_carriers"] = merged["Controls_with_variant_AOU"]
    merged["AOU_control_noncarriers"] = merged["Controls_without_variant_AOU"]

    merged["Total_case_carriers"] = merged["UKB_case_carriers"] + merged["AOU_case_carriers"]
    merged["Total_case_noncarriers"] = merged["UKB_case_noncarriers"] + merged["AOU_case_noncarriers"]
    merged["Total_control_carriers"] = merged["UKB_control_carriers"] + merged["AOU_control_carriers"]
    merged["Total_control_noncarriers"] = merged["UKB_control_noncarriers"] + merged["AOU_control_noncarriers"]
    merged["Dataset"] = "Discovery_Meta"
    merged["Method"] = "Fixed"

    return merged[[
        "Gene", "Dataset",
        "meta_OR", "meta_CI_lower", "meta_CI_upper", "meta_P_value",
        "meta_logOR", "meta_SE", "Q", "p_het", "Method",
        "UKB_case_carriers", "UKB_case_noncarriers",
        "UKB_control_carriers", "UKB_control_noncarriers",
        "AOU_case_carriers", "AOU_case_noncarriers",
        "AOU_control_carriers", "AOU_control_noncarriers",
        "Total_case_carriers", "Total_case_noncarriers",
        "Total_control_carriers", "Total_control_noncarriers",
    ]].copy()


def create_cohort_table(
    meta_df: pd.DataFrame,
    ukb: pd.DataFrame,
    aou: pd.DataFrame,
) -> pd.DataFrame:
    """Create a long table with UKB, AoU, and meta rows for each analyzed gene."""
    ukb_by_gene = ukb.set_index("Gene")
    aou_by_gene = aou.set_index("Gene")
    rows: list[dict[str, object]] = []

    for _, meta_row in meta_df.iterrows():
        gene = str(meta_row["Gene"])
        u = ukb_by_gene.loc[gene]
        a = aou_by_gene.loc[gene]

        rows.extend([
            {
                "Gene": gene,
                "Dataset": "UKB",
                "Cases_with_variant": u["Cases_with_variant"],
                "Cases_without_variant": u["Cases_without_variant"],
                "Controls_with_variant": u["Controls_with_variant"],
                "Controls_without_variant": u["Controls_without_variant"],
                "OR": u["Odds_Ratio"],
                "CI_lower": u["95%_CI_Lower"],
                "CI_upper": u["95%_CI_Upper"],
                "Pvalue": u["P_value"],
                "P_het": "",
                "Method": "",
            },
            {
                "Gene": gene,
                "Dataset": "AOU",
                "Cases_with_variant": a["Cases_with_variant"],
                "Cases_without_variant": a["Cases_without_variant"],
                "Controls_with_variant": a["Controls_with_variant"],
                "Controls_without_variant": a["Controls_without_variant"],
                "OR": a["Odds_Ratio"],
                "CI_lower": a["95%_CI_Lower"],
                "CI_upper": a["95%_CI_Upper"],
                "Pvalue": a["P_value"],
                "P_het": "",
                "Method": "",
            },
            {
                "Gene": gene,
                "Dataset": "Discovery_Meta",
                "Cases_with_variant": meta_row["Total_case_carriers"],
                "Cases_without_variant": meta_row["Total_case_noncarriers"],
                "Controls_with_variant": meta_row["Total_control_carriers"],
                "Controls_without_variant": meta_row["Total_control_noncarriers"],
                "OR": meta_row["meta_OR"],
                "CI_lower": meta_row["meta_CI_lower"],
                "CI_upper": meta_row["meta_CI_upper"],
                "Pvalue": meta_row["meta_P_value"],
                "P_het": meta_row["p_het"],
                "Method": meta_row["Method"],
            },
        ])

    return pd.DataFrame(rows)


def _format_output_pvalues(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = out[column].apply(
                lambda x: format_p_value(x) if not (isinstance(x, str) and x == "") else ""
            )
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fixed-effect inverse-variance-weighted discovery meta-analysis of "
            "UKB and AoU Firth logistic-regression results."
        )
    )
    parser.add_argument("--ukb", required=True, type=Path, help="UKB cohort result TSV/TXT.")
    parser.add_argument("--aou", required=True, type=Path, help="AoU cohort result TSV/TXT.")
    parser.add_argument(
        "--meta-output",
        type=Path,
        default=Path("discovery_meta_summary.tsv"),
        help="Output path for one-row-per-gene meta-analysis results.",
    )
    parser.add_argument(
        "--cohort-output",
        type=Path,
        default=Path("discovery_meta_with_cohorts.tsv"),
        help="Output path for the UKB/AoU/META long table.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        ukb = read_cohort_results(args.ukb, "UKB")
        aou = read_cohort_results(args.aou, "AoU")
        meta = meta_analysis(ukb, aou)

        if meta.empty:
            raise ValueError("No genes had valid results in both UKB and AoU.")

        # Sort using numeric high-precision P values before formatting for output.
        meta["_p_sort"] = meta["meta_P_value"].apply(float)
        meta = meta.sort_values(["_p_sort", "Gene"]).drop(columns="_p_sort").reset_index(drop=True)

        cohort_table = create_cohort_table(meta, ukb, aou)

        meta_out = meta.copy()
        for column in ["meta_OR", "meta_CI_lower", "meta_CI_upper", "meta_logOR", "meta_SE", "Q"]:
            meta_out[column] = pd.to_numeric(meta_out[column], errors="coerce").round(4)
        meta_out = _format_output_pvalues(meta_out, ["meta_P_value", "p_het"])
        cohort_out = _format_output_pvalues(cohort_table, ["Pvalue", "P_het"])

        meta_out.to_csv(args.meta_output, sep="\t", index=False)
        cohort_out.to_csv(args.cohort_output, sep="\t", index=False)

    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1

    print(f"Wrote discovery meta-analysis: {args.meta_output}")
    print(f"Wrote discovery cohort table: {args.cohort_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
