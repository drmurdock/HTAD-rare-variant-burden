#!/usr/bin/env python3
"""
Fixed-effect inverse-variance-weighted and DerSimonian-Laird random-effects
meta-analysis across discovery and replication cohorts.

Inputs are cohort-level TSV/TXT files with at least:
    gene, beta, se

A cohort label should be supplied in a ``Cohort`` or ``cohort`` column. If it
is absent, the script derives the cohort label from the filename prefix.
Optional columns used for annotation and carrier counts include:
    Phenotype, Variant_Class, p,
    carriers_case, carriers_ctrl, cases_total, controls_total

The analysis is restricted to genes represented in the discovery cohorts
(UKB and AoU). Overall meta-analysis results are reported only when UKB, AoU,
and at least one replication cohort contribute after any requested filtering.
The fixed-effect IVW analysis is the primary overall meta-analysis.
DerSimonian-Laird random-effects results are generated as a sensitivity analysis.

Replication support is evaluated among non-UKB/AoU cohorts. A gene has
independent replication support when at least one replication cohort has
P < 0.05 and an effect direction concordant with the fixed-effect UKB+AoU
discovery estimate.

Study-wide significance thresholds:
    missense:      2.6e-6
    LOF:           1.3e-5
    missense+LOF:  1.3e-5
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Optional: for heterogeneity p-value
try:
    from scipy.stats import chi2  # type: ignore
    _HAVE_SCIPY = True
except Exception:
    chi2 = None  # type: ignore
    _HAVE_SCIPY = False

Z975 = 1.96

# Significance thresholds
MISSENSE_SIG_THRESHOLD = 2.6e-6
LOF_SIG_THRESHOLD = 1.3e-5

# -------------------- Gene categorization helpers --------------------
ESTABLISHED_HTAD_GENES = {
    "FBN1", "ACTA2", "COL3A1", "MYH11", "SMAD3",
    "TGFB2", "TGFBR1", "TGFBR2", "MYLK", "LOX", "PRKG1"
}

PUTATIVE_HTAD_GENES = {
    "FBN2", "MFAP5", "COL5A1", "COL5A2", "BGN", "EFEMP2", "LTBP3", "ELN",
    "THSD4", "SKI", "SMAD2", "SMAD4", "TGFB3", "ARIH1", "FLNA", "MYH9",
    "HCN4", "MAT2A", "FOXE3", "MTOR", "CBS", "NOTCH1", "SLC2A10", "IPO8",
    "PKD1", "PKD2"
}

def categorize_gene(gene_name: str) -> str:
    g = (gene_name or "").upper()
    if g in ESTABLISHED_HTAD_GENES:
        return "Established"
    elif g in PUTATIVE_HTAD_GENES:
        return "Putative"
    else:
        return "Candidate"


# --------------------------- IO helpers ---------------------------
def _read_one(path: Path) -> pd.DataFrame:
    """Read one TSV/TXT. Validate basics and normalize columns we depend on."""
    df = pd.read_csv(path, sep="\t")

    # Basic validation
    for col in ["gene", "beta", "se"]:
        if col not in df.columns:
            raise ValueError(f"{path}: required column '{col}' not found; columns present: {list(df.columns)}")

    # Cohort column (prefer existing, else derive from filename)
    if "Cohort" in df.columns:
        df["cohort"] = df["Cohort"].astype(str)
    elif "cohort" in df.columns:
        df["cohort"] = df["cohort"].astype(str)
    else:
        stem = path.stem
        df["cohort"] = stem.split("_")[0] if "_" in stem else stem
        print(f"[WARNING] {path}: no Cohort column found; inferred cohort '{df['cohort'].iloc[0]}' from filename.")

    df["cohort"] = df["cohort"].astype(str).str.strip()

    duplicated = df.duplicated(subset=["gene", "cohort"], keep=False) if "gene" in df.columns else pd.Series(False, index=df.index)
    if duplicated.any():
        examples = df.loc[duplicated, ["gene", "cohort"]].drop_duplicates().head(10)
        pairs = ", ".join(f"{r.gene}/{r.cohort}" for r in examples.itertuples(index=False))
        raise ValueError(f"{path}: duplicate gene/cohort rows found (examples: {pairs}).")

    # Clean types
    df["gene"] = df["gene"].astype(str).str.strip()
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    if "p" in df.columns:
        df["p"] = pd.to_numeric(df["p"], errors="coerce")

    # Drop rows invalid for meta computation
    df = df[np.isfinite(df["beta"]) & np.isfinite(df["se"]) & (df["se"] > 0)].copy()

    # Keep all useful columns if present
    keep_cols = [
        "gene", "cohort", "Cohort", "Phenotype", "Variant_Class", "Variant Class",
        "beta", "se", "or", "ci_l", "ci_u", "p",
        "carriers_case", "carriers_ctrl", "cases_total", "controls_total",
        "var", "weight",
    ]
    cols = [c for c in keep_cols if c in df.columns]

    base = df[["gene", "cohort", "beta", "se"]].copy()
    extras = [c for c in cols if c not in base.columns]
    if extras:
        base = base.join(df[extras])
    return base


def _validate_discovery_annotations(all_df: pd.DataFrame) -> None:
    """Require UKB and AoU annotations to agree for each gene when supplied."""
    if "cohort" not in all_df.columns:
        return

    disc_mask = all_df["cohort"].astype(str).str.upper().isin(["UKB", "AOU"])
    discovery = all_df.loc[disc_mask].copy()

    for column in ["Phenotype", "Variant_Class", "Variant Class"]:
        if column not in discovery.columns:
            continue
        nonmissing = discovery.loc[discovery[column].notna(), ["gene", column]].copy()
        if nonmissing.empty:
            continue
        counts = nonmissing.groupby("gene")[column].nunique(dropna=True)
        bad_genes = counts[counts > 1].index.tolist()
        if bad_genes:
            preview = ", ".join(map(str, bad_genes[:10]))
            raise ValueError(
                f"UKB and AoU {column} annotations disagree for {len(bad_genes)} gene(s) "
                f"(examples: {preview})."
            )


def _fixed_effect_beta(block: pd.DataFrame) -> float:
    """Return an inverse-variance-weighted fixed-effect beta for a cohort block."""
    block = block[np.isfinite(block["beta"]) & np.isfinite(block["se"]) & (block["se"] > 0)]
    if block.empty:
        return np.nan
    weights = 1.0 / np.square(block["se"].to_numpy())
    if weights.sum() <= 0:
        return np.nan
    return float(np.sum(weights * block["beta"].to_numpy()) / weights.sum())


# ----------------------- statistical helpers ----------------------
def _two_tailed_p_from_z(z: np.ndarray) -> np.ndarray:
    """Two-sided p-values using erfc:
       p = 2 * (1 - Phi(|z|)) = erfc(|z|/sqrt(2))
    """
    absz = np.abs(z)
    erfc_vec = np.vectorize(lambda x: math.erfc(x), otypes=[float])
    return erfc_vec(absz / math.sqrt(2.0))


def _meta_one_gene_fixed(block: pd.DataFrame) -> pd.Series:
    """Compute fixed-effects IVW meta for a single gene block."""
    w = 1.0 / np.square(block["se"].to_numpy())
    beta = block["beta"].to_numpy()
    wsum = w.sum()
    if wsum <= 0:
        return pd.Series({})
    beta_hat = float(np.sum(w * beta) / wsum)
    se_hat = float(1.0 / math.sqrt(wsum))
    z = beta_hat / se_hat if se_hat > 0 else np.nan
    p = float(_two_tailed_p_from_z(np.array([z]))[0]) if np.isfinite(z) else np.nan

    # OR and CI on OR scale
    or_meta = math.exp(beta_hat)
    ci_l = math.exp(beta_hat - Z975 * se_hat)
    ci_u = math.exp(beta_hat + Z975 * se_hat)

    k = int(len(block))
    # Heterogeneity Q
    diff = beta - beta_hat
    Q = float(np.sum(w * np.square(diff))) if k > 1 else 0.0
    df_ = max(k - 1, 0)
    if _HAVE_SCIPY and df_ > 0:
        p_het = float(1.0 - chi2.cdf(Q, df_))
    else:
        p_het = np.nan
    I2 = 0.0
    if k > 1 and Q > 0:
        I2 = max(0.0, (Q - df_) / Q) * 100.0

    cohorts_list = ",".join(block["cohort"].astype(str).tolist())

    return pd.Series({
        "k": k,
        "beta_meta": beta_hat,
        "se_meta": se_hat,
        "z": z,
        "p": p,                # meta p
        "or_meta": or_meta,
        "ci_l": ci_l,          # meta CI (OR scale)
        "ci_u": ci_u,
        "Q": Q,
        "I2": I2,
        "p_het": p_het,
        "contributing_cohorts": cohorts_list,
        "Model": "Fixed",
    })


def _meta_one_gene_random(block: pd.DataFrame) -> pd.Series:
    """Compute random-effects (DerSimonian–Laird) meta for a single gene block."""
    beta = block["beta"].to_numpy()
    se = block["se"].to_numpy()
    w_fixed = 1.0 / np.square(se)
    wsum = w_fixed.sum()
    if wsum <= 0:
        return pd.Series({})

    # Fixed-effects pooled effect for Q
    beta_fixed = float(np.sum(w_fixed * beta) / wsum)
    k = int(len(block))
    df_ = max(k - 1, 0)
    # Cochran's Q under fixed
    Q = float(np.sum(w_fixed * np.square(beta - beta_fixed))) if k > 1 else 0.0
    if _HAVE_SCIPY and df_ > 0:
        p_het = float(1.0 - chi2.cdf(Q, df_))
    else:
        p_het = np.nan
    I2 = 0.0
    if k > 1 and Q > 0:
        I2 = max(0.0, (Q - df_) / Q) * 100.0

    # DerSimonian–Laird tau^2
    c = wsum - (np.square(w_fixed).sum() / wsum) if wsum > 0 else 0.0
    tau2 = max(0.0, (Q - df_) / c) if (k > 1 and c > 0) else 0.0

    # Random-effects weights
    w_star = 1.0 / (np.square(se) + tau2)
    w_star_sum = w_star.sum()
    beta_hat = float(np.sum(w_star * beta) / w_star_sum) if w_star_sum > 0 else np.nan
    se_hat = float(1.0 / math.sqrt(w_star_sum)) if w_star_sum > 0 else np.nan
    z = beta_hat / se_hat if se_hat and se_hat > 0 else np.nan
    p = float(_two_tailed_p_from_z(np.array([z]))[0]) if np.isfinite(z) else np.nan

    or_meta = math.exp(beta_hat) if np.isfinite(beta_hat) else np.nan
    ci_l = math.exp(beta_hat - Z975 * se_hat) if np.isfinite(beta_hat) and np.isfinite(se_hat) else np.nan
    ci_u = math.exp(beta_hat + Z975 * se_hat) if np.isfinite(beta_hat) and np.isfinite(se_hat) else np.nan

    cohorts_list = ",".join(block["cohort"].astype(str).tolist())

    return pd.Series({
        "k": k,
        "beta_meta": beta_hat,
        "se_meta": se_hat,
        "z": z,
        "p": p,
        "or_meta": or_meta,
        "ci_l": ci_l,
        "ci_u": ci_u,
        "Q": Q,
        "I2": I2,
        "p_het": p_het,
        "contributing_cohorts": cohorts_list,
        "Model": "Random",
    })


# ------------------------------ engine ------------------------------
def run_meta(
    inputs: List[Path],
    exclude_zero_case_carriers: bool,
    bonf_m: int,
    model: str = "Fixed"
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[str]]:
    """
    Compute meta-analysis (Fixed or Random) and return:
      meta_df, long_df, vclass_for_naming
    """
    frames: List[pd.DataFrame] = [_read_one(p) for p in inputs]
    all_df = pd.concat(frames, ignore_index=True)
    _validate_discovery_annotations(all_df)

    # Build a filtered copy for meta-analysis (sensitivity option)
    df_for_meta = all_df.copy()
    if exclude_zero_case_carriers and "carriers_case" in df_for_meta.columns:
        df_for_meta["carriers_case"] = pd.to_numeric(df_for_meta["carriers_case"], errors="coerce")
        # Exclude ONLY rows where carriers_case == 0 (keep NaN or >0)
        df_for_meta = df_for_meta[~(df_for_meta["carriers_case"].fillna(-1) == 0)].copy()

    # Choose per-gene function
    _per_gene = _meta_one_gene_fixed if str(model).lower().startswith("fix") else _meta_one_gene_random

    # Cross-version compatible GroupBy.apply (pandas ≥2.2 supports include_groups)
    try:
        meta = (
            df_for_meta
            .groupby("gene", sort=False)
            .apply(_per_gene, include_groups=False)
            .reset_index()
        )
    except TypeError:
        cols_for_apply = [c for c in df_for_meta.columns if c != "gene"]
        meta = (
            df_for_meta
            .groupby("gene", sort=False)[cols_for_apply]
            .apply(_per_gene)
            .reset_index()
        )

    # Bonferroni computed from numeric p where possible
    m_tests = float(bonf_m)
    if "p" in meta.columns:
        p_num = pd.to_numeric(meta["p"], errors="coerce")
        p_bonf = np.minimum(1.0, p_num * m_tests)
        meta["Bonferroni_P"] = p_bonf
        meta["Bonferroni_Significant"] = p_num.le(0.05 / m_tests)
    else:
        meta["Bonferroni_P"] = np.nan
        meta["Bonferroni_Significant"] = np.nan

    # Discovery phenotype & variant class from UKB/AoU (they should match)
    disc_mask = all_df["cohort"].astype(str).str.upper().isin(["UKB", "AOU"]) if "cohort" in all_df.columns else pd.Series(False, index=all_df.index)
    if "Phenotype" in all_df.columns:
        phen_map = (
            all_df.loc[disc_mask & all_df["Phenotype"].notna(), ["gene", "Phenotype"]]
                  .drop_duplicates(subset=["gene"]).set_index("gene")["Phenotype"]
        )
        meta["Discovery_Phenotype"] = meta["gene"].map(phen_map)
    else:
        meta["Discovery_Phenotype"] = np.nan

    varcol = "Variant_Class" if "Variant_Class" in all_df.columns else ("Variant Class" if "Variant Class" in all_df.columns else None)
    if varcol is not None:
        vclass_map = (
            all_df.loc[disc_mask & all_df[varcol].notna(), ["gene", varcol]]
                  .drop_duplicates(subset=["gene"]).set_index("gene")[varcol]
        )
        meta["Discovery_Variant_Class"] = meta["gene"].map(vclass_map)
    else:
        meta["Discovery_Variant_Class"] = np.nan

    # Aggregate total carriers/non-carriers (per gene) from cohorts,
    # respecting the exclude-zero-case-carriers flag
    have_counts = all(c in all_df.columns for c in ["carriers_case", "cases_total", "carriers_ctrl", "controls_total"])
    if have_counts:
        # When the flag is enabled, compute totals from the filtered df_for_meta
        use_filtered = (exclude_zero_case_carriers and "carriers_case" in all_df.columns)
        counts_source = df_for_meta.copy() if use_filtered else all_df.copy()

        # Ensure numeric
        for c in ["carriers_case", "cases_total", "carriers_ctrl", "controls_total"]:
            counts_source[c] = pd.to_numeric(counts_source[c], errors="coerce")

        sums = (
            counts_source
            .groupby("gene", as_index=False)
            .agg({
                "carriers_case": "sum",
                "cases_total": "sum",
                "carriers_ctrl": "sum",
                "controls_total": "sum",
            })
        )
        sums["total_noncarriers_case"] = sums["cases_total"] - sums["carriers_case"]
        sums["total_noncarriers_ctrl"] = sums["controls_total"] - sums["carriers_ctrl"]
        sums = sums.rename(columns={
            "carriers_case": "total_carriers_case",
            "carriers_ctrl": "total_carriers_ctrl",
        })
        meta = meta.merge(
            sums[[
                "gene",
                "total_carriers_case",
                "total_noncarriers_case",
                "total_carriers_ctrl",
                "total_noncarriers_ctrl",
            ]],
            on="gene",
            how="left",
        )
    else:
        meta["total_carriers_case"] = np.nan
        meta["total_noncarriers_case"] = np.nan
        meta["total_carriers_ctrl"] = np.nan
        meta["total_noncarriers_ctrl"] = np.nan


    # ---------------- Discovery-only inclusion & NA logic ----------------
    # Keep only discovery genes; compute meta only if UKB + AOU + ≥1 additional cohort present
    disc = {"UKB", "AOU"}
    disc_mask2_all = all_df["cohort"].astype(str).str.upper().isin(list(disc)) if "cohort" in all_df.columns else pd.Series(False, index=all_df.index)
    disc_genes = set(all_df.loc[disc_mask2_all, "gene"].dropna().astype(str).unique())

    # Restrict meta to discovery genes only
    meta = meta[meta["gene"].isin(disc_genes)].copy()

    # Determine cohort universe from df_for_meta (after exclusions, if any)
    all_cohorts_post = set(df_for_meta["cohort"].astype(str).str.upper().unique()) if "cohort" in df_for_meta.columns else set()
    addl_cohorts = sorted(all_cohorts_post - disc)

    def _split(s):
        if pd.isna(s) or not str(s).strip():
            return set()
        return {x.strip().upper() for x in str(s).split(",") if x.strip()}

    present_sets = meta["contributing_cohorts"].map(_split) if "contributing_cohorts" in meta.columns else pd.Series([set()] * len(meta), index=meta.index)
    has_both_disc = present_sets.map(lambda s: {"UKB", "AOU"}.issubset(s))
    n_addl_present = present_sets.map(lambda s: len(set(addl_cohorts) & s) if addl_cohorts else 0)
    ok_for_meta = has_both_disc & (n_addl_present >= 1)

    # For discovery genes that are NOT ok_for_meta, set meta fields to "NA" and print a message
    need_na = ~ok_for_meta
    if need_na.any():
        for g, cohorts in meta.loc[need_na, ["gene", "contributing_cohorts"]].itertuples(index=False):
            print(f"[INFO] Meta set to NA for gene {g}: no additional contributing cohort after filters (contributing: {cohorts})")

    for c in ["beta_meta","se_meta","z","p","or_meta","ci_l","ci_u","Q","I2","p_het"]:
        if c in meta.columns:
            meta[c] = meta[c].astype("object")
            meta.loc[need_na, c] = "NA"

    # Set Bonferroni fields to NA when the pooled meta-analysis is unavailable.
    if "Bonferroni_P" in meta.columns:
        meta["Bonferroni_P"] = meta["Bonferroni_P"].astype("object")
        mask_p_na_str = meta["p"].astype(str).eq("NA") if "p" in meta.columns else False
        meta.loc[need_na | mask_p_na_str, "Bonferroni_P"] = "NA"
    if "Bonferroni_Significant" in meta.columns:
        meta["Bonferroni_Significant"] = meta["Bonferroni_Significant"].astype("object")
        mask_p_na_str = meta["p"].astype(str).eq("NA") if "p" in meta.columns else False
        meta.loc[need_na | mask_p_na_str, "Bonferroni_Significant"] = "NA"

    # Sort meta by numeric p (coerce), then gene
    p_sort = pd.to_numeric(meta["p"], errors="coerce") if "p" in meta.columns else pd.Series([], dtype=float)
    meta = meta.assign(_p_sort=p_sort).sort_values(["_p_sort", "gene"], na_position="last").drop(columns=["_p_sort"]).reset_index(drop=True)

    # ---------- Replication support ----------
    # Discovery direction is the fixed-effect UKB+AoU beta for each gene.
    disc = {"UKB", "AOU"}
    cohorts_upper = df_for_meta["cohort"].astype(str).str.upper()
    discovery_rows = df_for_meta.loc[cohorts_upper.isin(list(disc)), ["gene", "beta", "se"]].copy()
    discovery_beta = discovery_rows.groupby("gene", sort=False)[["beta", "se"]].apply(_fixed_effect_beta)

    if "p" in df_for_meta.columns:
        rep = df_for_meta.loc[~cohorts_upper.isin(list(disc)), ["gene", "beta", "p"]].copy()
        rep["p"] = pd.to_numeric(rep["p"], errors="coerce")
        rep["beta"] = pd.to_numeric(rep["beta"], errors="coerce")
        rep["Discovery_Beta"] = rep["gene"].map(discovery_beta)
        rep["Nominal"] = rep["p"].lt(0.05)
        rep["Direction_Concordant"] = (
            np.sign(rep["beta"]) == np.sign(rep["Discovery_Beta"])
        ) & rep["beta"].notna() & rep["Discovery_Beta"].notna() & rep["Discovery_Beta"].ne(0)
        rep["Independent_Support"] = rep["Nominal"] & rep["Direction_Concordant"]

        any_rep_sig = rep.groupby("gene")["Nominal"].any() if not rep.empty else pd.Series(dtype=bool)
        any_rep_concordant = rep.groupby("gene")["Direction_Concordant"].any() if not rep.empty else pd.Series(dtype=bool)
        any_independent_support = rep.groupby("gene")["Independent_Support"].any() if not rep.empty else pd.Series(dtype=bool)

        meta["Nominal_Sig_in_Replication"] = meta["gene"].map(any_rep_sig).fillna(False)
        meta["Directionally_Concordant_Replication"] = meta["gene"].map(any_rep_concordant).fillna(False)
        meta["Independent_Replication_Support"] = meta["gene"].map(any_independent_support).fillna(False)
    else:
        meta["Nominal_Sig_in_Replication"] = False
        meta["Directionally_Concordant_Replication"] = False
        meta["Independent_Replication_Support"] = False

    # Arrange the meta-analysis summary columns.
    desired_front = [
        "gene",
        "Discovery_Phenotype",
        "Discovery_Variant_Class",
        "total_carriers_case",
        "total_noncarriers_case",
        "total_carriers_ctrl",
        "total_noncarriers_ctrl",
        "k",
        "contributing_cohorts",
        "meta_OR",
        "meta_CI_lower",
        "meta_CI_upper",
        "meta_P_value",
        "meta_SE",
        "Model",
        "Q",
        "I2",
        "p_het",
        "Bonferroni_P",
        "Bonferroni_Significant",
        "Nominal_Sig_in_Replication",
        "Directionally_Concordant_Replication",
        "Independent_Replication_Support",
    ]
    # Rename core meta-analysis fields before reordering.
    meta = meta.rename(columns={
        "or_meta": "meta_OR",
        "ci_l": "meta_CI_lower",
        "ci_u": "meta_CI_upper",
        "p": "meta_P_value",
        "se_meta": "meta_SE",
        "nominal_sig_in_replication": "Nominal_Sig_in_Replication",
    })
    meta = meta[[c for c in desired_front if c in meta.columns]]

    # Add gene category to the meta summary.
    meta.insert(1, "Gene_Category", meta["gene"].map(categorize_gene))

    # ---------------- Build the long table (cohort file) ----------------
    eligible_genes = set(meta["gene"].tolist())

    if exclude_zero_case_carriers and "carriers_case" in all_df.columns:
        cc_num = pd.to_numeric(all_df["carriers_case"], errors="coerce")
        long_source = all_df[~(cc_num.fillna(-1) == 0)].copy()
    else:
        long_source = all_df.copy()

    # Join meta annotations needed to construct the long table.
    meta_for_join = meta.drop(
        columns=[
            "meta_P_value", "meta_CI_lower", "meta_CI_upper",
            "Bonferroni_P", "Bonferroni_Significant",
            "Nominal_Sig_in_Replication",
            "Directionally_Concordant_Replication",
            "Independent_Replication_Support",
        ],
        errors="ignore",
    )
    long_individual = (
        long_source[long_source["gene"].isin(eligible_genes)]
        .merge(meta_for_join, on="gene", how="left")
    )

    # Construct one META row per gene with meta values copied
    long_cols = list(long_individual.columns)
    meta_rows = pd.DataFrame({c: np.nan for c in long_cols}, index=meta.index)
    meta_rows["gene"] = meta["gene"]

    # cohort labels
    if "cohort" in long_cols:
        meta_rows["cohort"] = "META"
    if "Cohort" in long_cols:
        meta_rows["Cohort"] = "META"

    # discovery phenotype/class (force META Phenotype = "TAD" if present)
    if "Phenotype" in long_cols:
        meta_rows["Phenotype"] = "TAD"
    if "Variant_Class" in long_cols and "Discovery_Variant_Class" in meta.columns:
        meta_rows["Variant_Class"] = meta["Discovery_Variant_Class"]

    # effect columns from meta (use renamed fields)
    if "se" in long_cols and "meta_SE" in meta.columns:
        meta_rows["se"] = meta["meta_SE"]
    if "or" in long_cols and "meta_OR" in meta.columns:
        meta_rows["or"] = meta["meta_OR"]
    if "ci_l" in long_cols and "meta_CI_lower" in meta.columns:
        meta_rows["ci_l"] = meta["meta_CI_lower"]
    if "ci_u" in long_cols and "meta_CI_upper" in meta.columns:
        meta_rows["ci_u"] = meta["meta_CI_upper"]
    if "p" in long_cols and "meta_P_value" in meta.columns:
        meta_rows["p"] = meta["meta_P_value"]

    # counts (use totals aggregated across cohorts)
    if "carriers_case" in long_cols and "total_carriers_case" in meta.columns:
        meta_rows["carriers_case"] = meta["total_carriers_case"]
    if "carriers_ctrl" in long_cols and "total_carriers_ctrl" in meta.columns:
        meta_rows["carriers_ctrl"] = meta["total_carriers_ctrl"]
    if "cases_total" in long_cols and {"total_carriers_case", "total_noncarriers_case"}.issubset(set(meta.columns)):
        meta_rows["cases_total"] = meta["total_carriers_case"] + meta["total_noncarriers_case"]
    if "controls_total" in long_cols and {"total_carriers_ctrl", "total_noncarriers_ctrl"}.issubset(set(meta.columns)):
        meta_rows["controls_total"] = meta["total_carriers_ctrl"] + meta["total_noncarriers_ctrl"]

    # derived var/weight if needed — be robust to "NA" strings
    if "se" in meta_rows.columns and "var" in long_cols:
        se_numeric = pd.to_numeric(meta_rows["se"], errors="coerce")
        meta_rows["var"] = np.square(se_numeric)
    if "var" in meta_rows.columns and "weight" in long_cols:
        denom = pd.to_numeric(meta_rows["var"], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            meta_rows["weight"] = 1.0 / denom

    # Include meta-only fields on META rows (populate values explicitly)
    if "k" in long_cols and "k" in meta.columns:
        meta_rows["k"] = meta["k"]
    if "Q" in long_cols and "Q" in meta.columns:
        meta_rows["Q"] = meta["Q"]
    if "I2" in long_cols and "I2" in meta.columns:
        meta_rows["I2"] = meta["I2"]
    if "p_het" in long_cols and "p_het" in meta.columns:
        meta_rows["p_het"] = meta["p_het"]
    if "contributing_cohorts" in long_cols and "contributing_cohorts" in meta.columns:
        meta_rows["contributing_cohorts"] = meta["contributing_cohorts"]
    if "Model" in meta.columns and "Model" in long_cols:
        meta_rows["Model"] = meta["Model"]

    # Combine cohort rows + meta rows
    long = pd.concat([long_individual, meta_rows], ignore_index=True, sort=False)

    # Single 'cohort' column (drop capitalized duplicate if present)
    if "Cohort" in long.columns and "cohort" in long.columns:
        long = long.drop(columns=["Cohort"])

    # Order rows within gene: UKB, AOU, others, META (last)
    cser = long["cohort"].astype(str).str.upper()
    prio = cser.map({"UKB": 0, "AOU": 1}).fillna(2)
    prio = prio.where(~cser.eq("META"), 99)
    long = long.assign(_prio=prio).sort_values(["gene", "_prio", "cohort"]).drop(columns=["_prio"]).reset_index(drop=True)

    # Propagate the pooled P value across long-table rows for convenient filtering.
    if "p" in long.columns and "meta_P_value" in meta.columns:
        meta_p_map = meta.set_index("gene")["meta_P_value"].to_dict()
        long["p_META"] = long["gene"].map(meta_p_map)

    # Cast some meta-only cols to object type, then blank for non-META rows
    if "cohort" in long.columns:
        nonmeta_mask = ~long["cohort"].astype(str).str.upper().eq("META")
        meta_only_cols = [
            "Discovery_Phenotype","Discovery_Variant_Class",
            "total_carriers_case","total_noncarriers_case","total_carriers_ctrl","total_noncarriers_ctrl",
            "k","Model","Q","I2","p_het","contributing_cohorts",
        ]
        cols_existing = [c for c in meta_only_cols if c in long.columns]
        if cols_existing:
            long[cols_existing] = long[cols_existing].astype("object")
            long.loc[nonmeta_mask, cols_existing] = ""

    # Final column order for cohort file
    desired_long_order = [
        "gene","cohort","beta","se","Phenotype","Variant_Class","or","ci_l","ci_u","p",
        "carriers_case","carriers_ctrl","cases_total","controls_total","var","weight",
        "Discovery_Phenotype","Discovery_Variant_Class",
        "total_carriers_case","total_noncarriers_case","total_carriers_ctrl","total_noncarriers_ctrl",
        "k","Model","Q","I2","p_het","contributing_cohorts","p_META"
    ]
    for c in desired_long_order:
        if c not in long.columns:
            long[c] = np.nan if c != "p_META" else ""
    long = long[desired_long_order]

    # ---------------- Meta-only output formatting ----------------
    # Use concise, publication-facing column names in the summary table.
    meta = meta.rename(columns={
        "Gene_Category": "Category",                 # Gene_Category -> Category
        "Discovery_Phenotype": "Phenotype",          # Discovery_Phenotype -> Phenotype
        "Discovery_Variant_Class": "Variant_Class",  # Discovery_Variant_Class -> Variant_Class
        "k": "k_Cohorts",                             # k -> k_Cohorts
        "contributing_cohorts": "Contributing_Cohorts",
        "meta_OR": "OR",                              # meta_OR -> OR
        "meta_CI_lower": "CI_lower",                  # meta_CI_lower -> CI_lower
        "meta_CI_upper": "CI_upper",                  # meta_CI_upper -> CI_upper
        "meta_P_value": "P",                          # meta_P_value -> P
    })

    # Add aggregate carrier and non-carrier totals.
    if {"total_carriers_case","total_carriers_ctrl","total_noncarriers_case","total_noncarriers_ctrl"}.issubset(set(meta.columns)):
        carriers = pd.to_numeric(meta["total_carriers_case"], errors="coerce") + pd.to_numeric(meta["total_carriers_ctrl"], errors="coerce")
        noncarriers = pd.to_numeric(meta["total_noncarriers_case"], errors="coerce") + pd.to_numeric(meta["total_noncarriers_ctrl"], errors="coerce")
        # Insert aggregate counts after Variant_Class.
        insert_idx = meta.columns.get_loc("Variant_Class") + 1 if "Variant_Class" in meta.columns else len(meta.columns)
        meta.insert(insert_idx, "Carriers", carriers)
        meta.insert(insert_idx + 1, "Non-Carriers", noncarriers)
    else:
        # Preserve the output schema when aggregate counts are unavailable.
        insert_idx = meta.columns.get_loc("Variant_Class") + 1 if "Variant_Class" in meta.columns else len(meta.columns)
        meta.insert(insert_idx, "Carriers", np.nan)
        meta.insert(insert_idx + 1, "Non-Carriers", np.nan)

    # Add formatted OR (95% CI) display column.
    def _fmt_sig(x):
        # Display OR and CI values to one decimal place; keep "NA" as-is.
        try:
            if pd.isna(x): return "NA"
            if isinstance(x, str) and x.upper() == "NA": return "NA"
            return f"{float(x):.1f}"
        except Exception:
            return "NA"

    or_idx = meta.columns.get_loc("OR") + 1 if "OR" in meta.columns else len(meta.columns)
    # Build the label
    or_str = meta["OR"].apply(_fmt_sig) if "OR" in meta.columns else pd.Series(["NA"]*len(meta))
    lo_str = meta["CI_lower"].apply(_fmt_sig) if "CI_lower" in meta.columns else pd.Series(["NA"]*len(meta))
    hi_str = meta["CI_upper"].apply(_fmt_sig) if "CI_upper" in meta.columns else pd.Series(["NA"]*len(meta))
    meta.insert(or_idx + 2, "OR (95% CI)", or_str + " (" + lo_str + "–" + hi_str + ")")

    # Keep p_het numeric (no string formatting) so Excel sees numbers, not text
    if "p_het" in meta.columns:
        meta["p_het"] = pd.to_numeric(meta["p_het"], errors="coerce")

    # Final column ordering for the meta summary.
    meta_order = [
        "gene",
        "Category",
        "Phenotype",
        "Variant_Class",
        "Carriers",
        "Non-Carriers",
        "total_carriers_case",
        "total_noncarriers_case",
        "total_carriers_ctrl",
        "total_noncarriers_ctrl",
        "k_Cohorts",
        "Contributing_Cohorts",
        "OR",
        "CI_lower",
        "CI_upper",
        "OR (95% CI)",
        "P",
        "meta_SE",
        "Model",
        "Q",
        "I2",
        "p_het",
        "Bonferroni_P",
        "Bonferroni_Significant",
        "Nominal_Sig_in_Replication",
        "Directionally_Concordant_Replication",
        "Independent_Replication_Support",
    ]
    meta = meta[[c for c in meta_order if c in meta.columns]]

    # Rename 'gene' to 'Gene' in both outputs
    meta = meta.rename(columns={"gene": "Gene"})
    long = long.rename(columns={"gene": "Gene"})

    # Determine a single discovery Variant_Class for naming (from UKB/AoU)
    vclass_for_naming: Optional[str] = None
    if "Variant_Class" in meta.columns:
        vals = meta["Variant_Class"].dropna().astype(str).unique().tolist()
        if len(vals) == 1:
            vclass_for_naming = vals[0]

    return meta, long, vclass_for_naming


def _apply_study_wide_sig(meta_df: pd.DataFrame) -> pd.DataFrame:
    """Add the manuscript study-wide significance flag from variant class and meta-analysis P value."""
    df = meta_df.copy()
    vclass_series = df.get("Variant_Class", pd.Series([np.nan] * len(df)))
    pvals = pd.to_numeric(df.get("P", pd.Series([np.nan] * len(df))), errors="coerce")

    # Normalize Variant Class
    vclass_norm = vclass_series.astype(str).str.strip().str.lower()
    missense_mask = vclass_norm.eq("missense")
    lof_or_combo_mask = vclass_norm.eq("lof") | vclass_norm.eq("missense+lof")

    sig = pd.Series(False, index=df.index)
    sig = sig | (missense_mask & pvals.le(MISSENSE_SIG_THRESHOLD))
    sig = sig | (lof_or_combo_mask & pvals.le(LOF_SIG_THRESHOLD))

    # If P is "NA" (string) from insufficient cohorts, leave as False
    df["Study-wide Significant"] = sig.fillna(False)
    return df


def _apply_cross_cohort_validation(meta_df: pd.DataFrame) -> pd.DataFrame:
    """Flag genes meeting both replication-support and study-wide significance criteria."""
    df = meta_df.copy()
    replication = df.get(
        "Independent_Replication_Support",
        pd.Series(False, index=df.index),
    ).fillna(False).astype(bool)
    study_wide = df.get(
        "Study-wide Significant",
        pd.Series(False, index=df.index),
    ).fillna(False).astype(bool)
    df["Cross_Cohort_Validation"] = replication & study_wide
    return df


# ------------------------------- CLI -------------------------------
def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fixed-effect IVW and DerSimonian-Laird random-effects meta-analysis "
            "across UKB, AoU, and replication cohort result files."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path, help="Input TSV/TXT files (>=2) with at least: gene,beta,se. 'Cohort' used if present.")
    parser.add_argument("-o", "--output", type=Path, default=Path("meta_results.tsv"), help="Base path for the meta-only TSV (model + variant class + bonf_m + label will be inserted).")
    parser.add_argument("--label", type=str, default=None, help="Suffix label appended before extension (e.g., 'candidates' → meta_results_missense+lof_candidates.tsv)")
    parser.add_argument(
        "--exclude-zero-case-carriers",
        action="store_true",
        help=(
            "Exclude cohort-gene rows with carriers_case == 0 from meta-analysis and "
            "with_cohorts output. Use this for the manuscript replication/meta-analysis workflow."
        ),
    )
    parser.add_argument("--bonf_m", type=int, default=72, help="Number of tests for Bonferroni correction (default: 72)."
)

    args = parser.parse_args(argv)

    # Print which Bonferroni correction value is being used
    print(f"[INFO] Using Bonferroni correction with {args.bonf_m} tests")

    if len(args.inputs) < 2:
        print("ERROR: Provide at least two input files.", file=sys.stderr)
        return 2

    try:
        # Run Fixed
        meta_fixed, long_fixed, vclass_for_naming_fixed = run_meta(
            args.inputs, args.exclude_zero_case_carriers, args.bonf_m, model="Fixed"
        )
        # Run Random
        meta_random, long_random, vclass_for_naming_random = run_meta(
            args.inputs, args.exclude_zero_case_carriers, args.bonf_m, model="Random"
        )

        # Apply manuscript significance and cross-cohort validation criteria.
        meta_fixed = _apply_cross_cohort_validation(_apply_study_wide_sig(meta_fixed))
        meta_random = _apply_cross_cohort_validation(_apply_study_wide_sig(meta_random))

    except Exception as e:
        print(f"ERROR while processing: {e}", file=sys.stderr)
        return 3

    # Helper to build path with model, vclass, bonf_m, and optional label
    def build_paths(base: Path, model_name: str, vclass_for_naming: Optional[str]) -> tuple[Path, Path]:
        meta_path = base
        # Insert variant class and Bonferroni test count when available.
        if vclass_for_naming and vclass_for_naming not in meta_path.stem:
            meta_path = meta_path.with_name(f"{meta_path.stem}_{vclass_for_naming}_{args.bonf_m}{meta_path.suffix}")
        else:
            # ensure bonf_m inclusion if not already there and vclass absent
            if str(args.bonf_m) not in meta_path.stem:
                meta_path = meta_path.with_name(f"{meta_path.stem}_{args.bonf_m}{meta_path.suffix}")
        if args.label:
            meta_path = meta_path.with_name(f"{meta_path.stem}_{args.label}{meta_path.suffix}")
        # Insert the model token before the extension.
        meta_path = meta_path.with_name(f"{meta_path.stem}_{model_name}{meta_path.suffix}")

        long_path = meta_path.with_suffix("")
        long_path = long_path.with_name(long_path.name + ".with_cohorts.tsv")
        return meta_path, long_path

    # Build paths
    meta_path_fixed, long_path_fixed = build_paths(args.output, "fixed", vclass_for_naming_fixed)
    meta_path_random, long_path_random = build_paths(args.output, "random_DL", vclass_for_naming_random)

    # Write TSV + XLSX for Fixed
    try:
        meta_fixed.to_csv(meta_path_fixed, sep="\t", index=False)
        long_fixed.to_csv(long_path_fixed, sep="\t", index=False)
        meta_fixed.to_excel(meta_path_fixed.with_suffix(".xlsx"), index=False)
        long_fixed.to_excel(long_path_fixed.with_suffix(".xlsx"), index=False)
    except Exception as e:
        print(f"WARNING: Failed to write one or more FIXED output files: {e}", file=sys.stderr)

    # Write TSV + XLSX for Random
    try:
        meta_random.to_csv(meta_path_random, sep="\t", index=False)
        long_random.to_csv(long_path_random, sep="\t", index=False)
        meta_random.to_excel(meta_path_random.with_suffix(".xlsx"), index=False)
        long_random.to_excel(long_path_random.with_suffix(".xlsx"), index=False)
    except Exception as e:
        print(f"WARNING: Failed to write one or more RANDOM output files: {e}", file=sys.stderr)

    print(f"Wrote FIXED meta-only: {len(meta_fixed):,} genes → {meta_path_fixed}")
    print(f"Wrote FIXED long table: {len(long_fixed):,} rows → {long_path_fixed}")
    print(f"Wrote RANDOM meta-only: {len(meta_random):,} genes → {meta_path_random}")
    print(f"Wrote RANDOM long table: {len(long_random):,} rows → {long_path_random}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())




