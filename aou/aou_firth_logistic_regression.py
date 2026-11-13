#!/usr/bin/env python3
"""
All of Us gene-based rare-variant burden analysis using Firth penalized logistic regression.

The script is designed to run as a worker job in the All of Us Researcher
Workbench. File paths are provided through environment variables by the job
submission wrapper.

Required environment variables
------------------------------
GENE_FILE
    Tab-delimited gene carrier file containing: gene, carrier_samples.
CASES_FILE
    Tab-delimited case cohort file containing: person_id.
CONTROLS_FILE
    Tab-delimited control cohort file containing: person_id.
COVARIATES_FILE
    Tab-delimited covariate file containing person_id, age_at_enrollment,
    sex_at_birth, pca_1, pca_2, pca_3, pca_4, and pca_5. If ancestry filtering
    is requested, ancestry_pred is also required.
OUTPUT_FILE
    Local path for the regression results.
ANCESTRY_PRED
    An ancestry code or "all".

Optional environment variables
------------------------------
RELATED_FILE
    Tab-delimited file containing person_id values to exclude.
FIRTHLOGIST_ZIP
    Localized path to firthlogist.zip. This is used only if firthlogist is not
    already importable in the runtime environment.
"""

from __future__ import annotations

import os
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_ENV = (
    "GENE_FILE",
    "CASES_FILE",
    "CONTROLS_FILE",
    "COVARIATES_FILE",
    "OUTPUT_FILE",
    "ANCESTRY_PRED",
)

INDEPENDENT_VARS = [
    "carrier",
    "age_at_enrollment",
    "sex_binary",
    "pca_1",
    "pca_2",
    "pca_3",
    "pca_4",
    "pca_5",
]

OUTPUT_HEADER = (
    "Gene\tCases_with_variant\tCases_without_variant\t"
    "Controls_with_variant\tControls_without_variant\t"
    "Odds_Ratio\t95%_CI_Lower\t95%_CI_Upper\tP_value\n"
)


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        fail(f"{label} is missing required column(s): {', '.join(missing)}")


def load_firthlogist():
    """Import firthlogist, falling back to a localized ZIP if necessary."""
    try:
        from firthlogist import FirthLogisticRegression

        return FirthLogisticRegression
    except ImportError:
        pass

    zip_value = os.environ.get("FIRTHLOGIST_ZIP")
    if not zip_value:
        fail(
            "firthlogist is not installed and FIRTHLOGIST_ZIP was not provided. "
            "Install firthlogist or localize firthlogist.zip with the job launcher."
        )

    zip_path = Path(zip_value)
    if not zip_path.is_file():
        fail(f"FIRTHLOGIST_ZIP does not point to a readable file: {zip_path}")

    extract_dir = Path(tempfile.mkdtemp(prefix="firthlogist_"))
    try:
        with zipfile.ZipFile(zip_path, "r") as archive:
            archive.extractall(extract_dir)
    except (zipfile.BadZipFile, OSError) as exc:
        fail(f"Unable to extract FIRTHLOGIST_ZIP: {exc}")

    # Locate the package without assuming a particular archive root directory.
    package_dirs = list(extract_dir.rglob("firthlogist/__init__.py"))
    if package_dirs:
        package_parent = package_dirs[0].parent.parent
        sys.path.insert(0, str(package_parent))
    else:
        # Some distributions may provide a module file instead of a package.
        module_files = list(extract_dir.rglob("firthlogist.py"))
        if not module_files:
            fail("Could not locate the firthlogist package inside FIRTHLOGIST_ZIP.")
        sys.path.insert(0, str(module_files[0].parent))

    try:
        from firthlogist import FirthLogisticRegression

        return FirthLogisticRegression
    except ImportError as exc:
        fail(f"Unable to import firthlogist after extracting FIRTHLOGIST_ZIP: {exc}")


def format_p_value(p_value: float) -> str:
    return f"{p_value:.2e}" if p_value < 0.001 else f"{p_value:.4f}"


def no_case_carrier_line(
    gene: str,
    cases_with_variant: int,
    cases_without_variant: int,
    controls_with_variant: int,
    controls_without_variant: int,
) -> str:
    return (
        f"{gene}\t{cases_with_variant}\t{cases_without_variant}\t"
        f"{controls_with_variant}\t{controls_without_variant}\t"
        "NA\tNA\tNA\tNA\n"
    )


def count_group(df: pd.DataFrame, phenotype: int, carrier: int) -> int:
    return int(((df["phenotype"] == phenotype) & (df["carrier"] == carrier)).sum())


def main() -> int:
    missing_env = [name for name in REQUIRED_ENV if not os.environ.get(name)]
    if missing_env:
        fail(f"Missing required environment variable(s): {', '.join(missing_env)}")

    gene_file = Path(os.environ["GENE_FILE"])
    cases_file = Path(os.environ["CASES_FILE"])
    controls_file = Path(os.environ["CONTROLS_FILE"])
    covariates_file = Path(os.environ["COVARIATES_FILE"])
    output_file = Path(os.environ["OUTPUT_FILE"])
    related_value = os.environ.get("RELATED_FILE")
    related_file = Path(related_value) if related_value else None
    ancestry_filter = os.environ["ANCESTRY_PRED"].strip()

    FirthLogisticRegression = load_firthlogist()

    gene_df = pd.read_csv(gene_file, sep="\t")
    cases = pd.read_csv(cases_file, sep="\t")
    controls = pd.read_csv(controls_file, sep="\t")
    covariates = pd.read_csv(covariates_file, sep="\t")

    require_columns(gene_df, ["gene", "carrier_samples"], "GENE_FILE")
    require_columns(cases, ["person_id"], "CASES_FILE")
    require_columns(controls, ["person_id"], "CONTROLS_FILE")

    required_covariates = [
        "person_id",
        "age_at_enrollment",
        "sex_at_birth",
        "pca_1",
        "pca_2",
        "pca_3",
        "pca_4",
        "pca_5",
    ]
    if ancestry_filter.lower() != "all":
        required_covariates.append("ancestry_pred")
    require_columns(covariates, required_covariates, "COVARIATES_FILE")

    cases = cases.copy()
    controls = controls.copy()
    cases["phenotype"] = 1
    controls["phenotype"] = 0
    pheno = pd.concat([cases, controls], ignore_index=True)

    # Retain the analysis behavior of excluding IDs listed in RELATED_FILE.
    if related_file is not None:
        related = pd.read_csv(related_file, sep="\t")
        require_columns(related, ["person_id"], "RELATED_FILE")
        related_ids = set(related["person_id"])
        pheno = pheno[~pheno["person_id"].isin(related_ids)].copy()

    df = pd.merge(pheno, covariates, on="person_id", how="inner")
    df["person_id"] = pd.to_numeric(df["person_id"], errors="coerce")
    df["sex_binary"] = df["sex_at_birth"].map({"male": 1, "female": 0})

    if ancestry_filter.lower() != "all":
        df = df[
            df["ancestry_pred"].astype(str).str.lower() == ancestry_filter.lower()
        ].copy()

    results: list[str] = []

    for _, row in gene_df.iterrows():
        gene_name = str(row["gene"])
        carrier_samples = row["carrier_samples"]

        if pd.isna(carrier_samples) or not str(carrier_samples).strip():
            print(f"Skipping {gene_name}: no carrier data")
            continue

        try:
            carrier_set = {
                int(sample.strip())
                for sample in str(carrier_samples).split(",")
                if sample.strip()
            }
        except ValueError:
            print(f"Skipping {gene_name}: unable to parse carrier samples")
            continue

        df["carrier"] = df["person_id"].isin(carrier_set).astype(int)

        cases_with_variant = count_group(df, phenotype=1, carrier=1)
        cases_without_variant = count_group(df, phenotype=1, carrier=0)
        controls_with_variant = count_group(df, phenotype=0, carrier=1)
        controls_without_variant = count_group(df, phenotype=0, carrier=0)

        if cases_with_variant == 0:
            line = no_case_carrier_line(
                gene_name,
                cases_with_variant,
                cases_without_variant,
                controls_with_variant,
                controls_without_variant,
            )
            results.append(line)
            print(line, end="")
            continue

        X = df[INDEPENDENT_VARS]
        y = df["phenotype"]

        try:
            model = FirthLogisticRegression(
                skip_ci=False,
                skip_pvals=False,
                max_iter=50,
            )
            model.fit(X, y)
        except Exception as exc:
            print(f"Error fitting model for {gene_name}: {exc}")
            continue

        coefficients = model.coef_
        p_values = model.pvals_
        ci_lower_values = (
            np.exp(model.ci_[:, 0])
            if hasattr(model, "ci_")
            else ["NA"] * len(coefficients)
        )
        ci_upper_values = (
            np.exp(model.ci_[:, 1])
            if hasattr(model, "ci_")
            else ["NA"] * len(coefficients)
        )

        carrier_index = INDEPENDENT_VARS.index("carrier")
        odds_ratio = np.exp(coefficients[carrier_index])
        ci_lower = ci_lower_values[carrier_index]
        ci_upper = ci_upper_values[carrier_index]
        p_value = p_values[carrier_index]

        result_line = (
            f"{gene_name}\t"
            f"{cases_with_variant}\t"
            f"{cases_without_variant}\t"
            f"{controls_with_variant}\t"
            f"{controls_without_variant}\t"
            f"{odds_ratio:.4f}\t"
            f"{ci_lower:.4f}\t"
            f"{ci_upper:.4f}\t"
            f"{format_p_value(p_value)}\n"
        )
        results.append(result_line)
        print(result_line, end="")

    # Keep genes with NA results at the end, otherwise sort by P value.
    results.sort(
        key=lambda line: (
            float(line.rstrip("\n").split("\t")[-1])
            if line.rstrip("\n").split("\t")[-1] != "NA"
            else float("inf")
        )
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as handle:
        handle.write(OUTPUT_HEADER)
        handle.writelines(results)

    print(f"Wrote {len(results)} gene results to {output_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
