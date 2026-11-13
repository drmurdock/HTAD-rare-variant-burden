#!/usr/bin/env python3
"""
UK Biobank gene-based rare-variant burden analysis using Firth penalized logistic regression.

Analysis supporting the manuscript:
"Gene-Based Rare Variant Burden Analyses Across Biobanks Identify Novel
High-Risk Genes for Thoracic Aortic Disease."

This script is designed for use in UKB Research Analysis Platform / DNAnexus
environments. The statistical model is:

    case/control status ~ Sex + Age + High_risk_variant + PC1-PC4 [+ HTN]

HTN is included by default and is removed only when --no_htn is specified.

Phenotype codes:
    A = Dissection
    B = TAA
    C = TAA-surgery
    D = Dissection + TAA-surgery

The input variant file is expected to have already been filtered to the desired
variant class/mask (for example missense, LOF, or missense + LOF). This script
does not redefine variant masks.
"""

import argparse
import logging
import os
import sys
from multiprocessing import Manager, Pool, cpu_count

import numpy as np
import pandas as pd
from firthlogist import FirthLogisticRegression
from tqdm import tqdm


# UKB phenotype filenames used for the manuscript analysis. The default root is
# a DNAnexus-mounted project directory and can be changed with --data-dir.
CASE_FILENAMES = {
    'A': os.path.join(
        'eid_latest',
        'eid_dissection_510.demographics_icd10.marfan_euro.486.txt',
    ),
    'B': os.path.join(
        'eid_latest',
        'eid_TAA_no_rupture_no_dissection_1152.demographics_icd10.marfan_euro_1104.txt',
    ),
    'C': os.path.join(
        'eid_latest',
        'eid_TAA_no_rupture_no_dissection_open_surgery_271.demographics_icd10.marfan_euro_262.txt',
    ),
    # D = combined Dissection + TAA-surgery phenotype.
    'D': os.path.join(
        'eid_latest',
        'eid_dissection_taa_open_surgery.demographics_icd10.748.txt',
    ),
}

PHENOTYPE_PREFIX = {
    'A': 'dissection',
    'B': 'TAA',
    'C': 'TAA+surgery',
    'D': 'dissection+TAA+surgery',
}

DEFAULT_CONTROLS_FILENAME = (
    'eid_UKB_full_dataset_490541.filtered_447570.demographics_icd10.marfan.txt'
)
DEFAULT_PCS_FILENAME = 'UKB.PRS_principal_components.csv'
DEFAULT_EXCLUSION_FILENAME = os.path.join('eid_latest', 'ukb_related_to_remove.txt')


def setup_logging(include_htn, prefix, output_dir):
    """Set up file logging."""
    log_filename = f"{prefix}_logistic_regression_firths_PC"
    if include_htn:
        log_filename += "_HTN"
    log_filename += ".log"
    log_path = os.path.join(output_dir, log_filename)

    # Important when the function is reused in notebooks or interactive jobs.
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        filename=os.path.abspath(log_path),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        force=True,
    )
    logging.info("Starting analysis with logging to %s", os.path.abspath(log_path))


def load_exclusion_eids(exclusion_file):
    """Load EIDs to exclude from the analysis."""
    try:
        with open(exclusion_file, 'r') as f:
            return {line.strip() for line in f if line.strip()}
    except Exception as exc:
        logging.error("Error loading exclusion EIDs file: %s", exc)
        raise


def require_columns(df, required_columns, context):
    """Fail early with a clear error if a required input column is missing."""
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s) in {context}: {missing}")


def parse_variant_eids(values):
    """Return cleaned carrier EIDs from a Series of comma-separated EID fields."""
    eids = set()
    for value in values.dropna():
        eids.update(
            eid.strip()
            for eid in str(value).split(',')
            if eid and eid.strip()
        )
    return eids


def prepare_logistic_regression_input(
    variant_file,
    gene,
    cases_file,
    controls_file,
    pcs_file,
    include_htn=False,
    exclusion_eids=None,
):
    """Construct the per-gene case/control dataframe used by the Firth model."""
    try:
        controls = pd.read_csv(controls_file, sep='\t', dtype={'EID': str})
        cases = pd.read_csv(cases_file, sep='\t', dtype={'EID': str})
        pcs = pd.read_csv(pcs_file, dtype={'EID': str})

        controls.columns = controls.columns.str.strip()
        cases.columns = cases.columns.str.strip()
        pcs.columns = pcs.columns.str.strip()

        require_columns(controls, ['EID', 'Sex', 'Age'], 'controls file')
        require_columns(cases, ['EID', 'Sex', 'Age'], 'cases file')
        require_columns(pcs, ['EID', 'PC_1', 'PC_2', 'PC_3', 'PC_4'], 'PC file')

        if include_htn:
            require_columns(controls, ['HTN_dx'], 'controls file')
            require_columns(cases, ['HTN_dx'], 'cases file')

        # Normalize IDs before exclusions/merges.
        controls['EID'] = controls['EID'].astype(str).str.strip()
        cases['EID'] = cases['EID'].astype(str).str.strip()
        pcs['EID'] = pcs['EID'].astype(str).str.strip()

        # Apply the optional relatedness exclusion.
        if exclusion_eids:
            initial_cases_count = len(cases)
            initial_controls_count = len(controls)

            cases = cases[~cases['EID'].isin(exclusion_eids)]
            controls = controls[~controls['EID'].isin(exclusion_eids)]

            logging.info(
                "Gene: %s - Excluded %s cases and %s controls based on exclusion list",
                gene,
                initial_cases_count - len(cases),
                initial_controls_count - len(controls),
            )

        controls['Dissection_dx'] = 0
        cases['Dissection_dx'] = 1
        cohort_data = pd.concat([controls, cases], ignore_index=True)

        variant_data = pd.read_csv(variant_file, sep='\t', dtype={'EIDs': str})
        variant_data.columns = variant_data.columns.str.strip()
        require_columns(variant_data, ['Gene', 'EIDs'], 'variant file')
        variant_data = variant_data[variant_data['Gene'] == gene]


        eids_with_variant = parse_variant_eids(variant_data['EIDs'])

        # Individuals without PC data are not retained in the analytic cohort.
        cohort_data = pd.merge(cohort_data, pcs, on='EID', how='inner')

        cohort_data['High_risk_variant'] = cohort_data['EID'].isin(eids_with_variant).astype(int)
        cohort_data['Sex'] = cohort_data['Sex'].map({'female': 0, 'male': 1})

        if include_htn:
            cohort_data['HTN_dx'] = cohort_data['HTN_dx'].apply(
                lambda value: 1 if value == 'HTN' else 0
            )

        initial_sample_count = cohort_data.shape[0]
        initial_cohort_data = cohort_data.copy()

        required_columns = [
            'Sex',
            'Age',
            'High_risk_variant',
            'PC_1',
            'PC_2',
            'PC_3',
            'PC_4',
        ]
        if include_htn:
            required_columns.append('HTN_dx')

        cohort_data = cohort_data.dropna(subset=required_columns)

        final_eids = cohort_data['EID']
        excluded_samples = initial_cohort_data[
            ~initial_cohort_data['EID'].isin(final_eids)
        ]
        excluded_cases = excluded_samples[excluded_samples['Dissection_dx'] == 1]
        excluded_controls = excluded_samples[excluded_samples['Dissection_dx'] == 0]

        final_data_columns = [
            'EID',
            'Sex',
            'Age',
            'Dissection_dx',
            'High_risk_variant',
            'PC_1',
            'PC_2',
            'PC_3',
            'PC_4',
        ]
        if include_htn:
            final_data_columns.append('HTN_dx')
        final_data = cohort_data[final_data_columns]

        cases_with_variant = final_data[
            (final_data['Dissection_dx'] == 1)
            & (final_data['High_risk_variant'] == 1)
        ].shape[0]
        cases_without_variant = final_data[
            (final_data['Dissection_dx'] == 1)
            & (final_data['High_risk_variant'] == 0)
        ].shape[0]
        controls_with_variant = final_data[
            (final_data['Dissection_dx'] == 0)
            & (final_data['High_risk_variant'] == 1)
        ].shape[0]
        controls_without_variant = final_data[
            (final_data['Dissection_dx'] == 0)
            & (final_data['High_risk_variant'] == 0)
        ].shape[0]

        logging.info("Gene: %s", gene)
        logging.info("Total initial samples: %s", initial_sample_count)
        logging.info("Samples excluded due to missing data: %s", excluded_samples.shape[0])
        logging.info("Excluded cases: %s", excluded_cases.shape[0])
        logging.info("Excluded controls: %s", excluded_controls.shape[0])
        logging.info("EIDs of excluded cases: %s", excluded_cases['EID'].tolist())
        logging.info(
            "Total cases after exclusion: %s",
            final_data[final_data['Dissection_dx'] == 1].shape[0],
        )
        logging.info(
            "Total controls after exclusion: %s",
            final_data[final_data['Dissection_dx'] == 0].shape[0],
        )
        logging.info("Cases with high-risk variant: %s", cases_with_variant)
        logging.info("Cases without high-risk variant: %s", cases_without_variant)
        logging.info("Controls with high-risk variant: %s", controls_with_variant)
        logging.info("Controls without high-risk variant: %s", controls_without_variant)

        if include_htn:
            htn_cases = final_data[
                (final_data['Dissection_dx'] == 1) & (final_data['HTN_dx'] == 1)
            ].shape[0]
            htn_controls = final_data[
                (final_data['Dissection_dx'] == 0) & (final_data['HTN_dx'] == 1)
            ].shape[0]
            logging.info("Hypertensive cases: %s", htn_cases)
            logging.info("Hypertensive controls: %s", htn_controls)

        return (
            final_data,
            cases_with_variant,
            cases_without_variant,
            controls_with_variant,
            controls_without_variant,
        )

    except Exception as exc:
        logging.exception("Error processing gene %s: %s", gene, exc)
        return None, 0, 0, 0, 0


def perform_firth_logistic_regression(data, gene, include_htn=False, max_iter=50):
    """Fit the same Firth logistic-regression model used in the manuscript."""
    independent_vars = [
        'Age',
        'High_risk_variant',
        'PC_1',
        'PC_2',
        'PC_3',
        'PC_4',
    ]
    if 'Sex' in data.columns:
        independent_vars.insert(0, 'Sex')
    if include_htn:
        independent_vars.append('HTN_dx')

    X = data[independent_vars].values
    y = data['Dissection_dx'].values

    flr_model = FirthLogisticRegression(
        skip_ci=False,
        skip_pvals=False,
        max_iter=max_iter,
    )
    flr_model.fit(X, y)

    coef = flr_model.coef_
    pvals = flr_model.pvals_
    ci_lower = (
        np.exp(flr_model.ci_[:, 0])
        if hasattr(flr_model, 'ci_')
        else ['NA'] * len(coef)
    )
    ci_upper = (
        np.exp(flr_model.ci_[:, 1])
        if hasattr(flr_model, 'ci_')
        else ['NA'] * len(coef)
    )

    summary = f"\n===== Firth Logistic Regression Summary for {gene} =====\n"
    summary += "Variable\tCoefficient\tP-value\t95% CI Lower\t95% CI Upper\n"
    for i, var in enumerate(independent_vars):
        summary += (
            f"{var}\t{coef[i]:.4f}\t{format_p_value(pvals[i])}\t"
            f"{ci_lower[i]:.4f}\t{ci_upper[i]:.4f}\n"
        )
    summary += "\n"

    hrv_index = independent_vars.index('High_risk_variant')
    return (
        np.exp(coef[hrv_index]),
        ci_lower[hrv_index],
        ci_upper[hrv_index],
        pvals[hrv_index],
        summary,
    )


def handle_no_variants(
    gene,
    cases_with_variant,
    cases_without_variant,
    controls_with_variant,
    controls_without_variant,
):
    return (
        f"{gene}\t{cases_with_variant}\t{cases_without_variant}\t"
        f"{controls_with_variant}\t{controls_without_variant}\t"
        "NA\tNA\tNA\tNA\n"
    )


def format_p_value(p_value):
    """Use scientific notation only for P < 0.001."""
    if p_value < 0.001:
        return f"{p_value:.2e}"
    return f"{p_value:.4f}"


def process_gene(args):
    (
        gene,
        variant_file,
        cases_file,
        controls_file,
        pcs_file,
        include_htn,
        summary_queue,
        high_risk_queue,
        exclusion_eids,
    ) = args

    try:
        (
            final_data,
            cases_with_variant,
            cases_without_variant,
            controls_with_variant,
            controls_without_variant,
        ) = prepare_logistic_regression_input(
            variant_file,
            gene,
            cases_file,
            controls_file,
            pcs_file,
            include_htn,
            exclusion_eids,
        )

        if final_data is None:
            return handle_no_variants(gene, 0, 0, 0, 0)

        if cases_with_variant == 0:
            result = handle_no_variants(
                gene,
                cases_with_variant,
                cases_without_variant,
                controls_with_variant,
                controls_without_variant,
            )
            summary_queue.put("")
        else:
            odds_ratio, ci_lower, ci_upper, p_value, summary = (
                perform_firth_logistic_regression(final_data, gene, include_htn)
            )
            result = (
                f"{gene}\t{cases_with_variant}\t{cases_without_variant}\t"
                f"{controls_with_variant}\t{controls_without_variant}\t"
                f"{odds_ratio:.4f}\t{ci_lower:.4f}\t{ci_upper:.4f}\t"
                f"{format_p_value(p_value)}\n"
            )
            summary_queue.put(summary)

        # Write high-risk variant details for carriers in the analytic cohort.
        gene_variants = pd.read_csv(variant_file, sep='\t', dtype={'EIDs': str})
        gene_variants.columns = gene_variants.columns.str.strip()
        gene_variants = gene_variants[gene_variants['Gene'] == gene]

        # Variant annotations are supplementary output fields rather than model
        # inputs. Missing annotation columns are reported as NA and do not prevent
        # the gene-level Firth analysis from completing.
        require_columns(gene_variants, ['EIDs'], 'variant file')

        case_eid_set = set(
            final_data.loc[final_data['Dissection_dx'] == 1, 'EID'].astype(str)
        )
        control_eid_set = set(
            final_data.loc[final_data['Dissection_dx'] == 0, 'EID'].astype(str)
        )

        for _, variant_row in gene_variants.iterrows():
            eid_value = variant_row['EIDs']
            if pd.isna(eid_value):
                continue

            eids_with_variant = {
                eid.strip()
                for eid in str(eid_value).split(',')
                if eid and eid.strip()
            }

            cases_eids = sorted(eids_with_variant.intersection(case_eid_set))
            controls_eids = sorted(eids_with_variant.intersection(control_eid_set))

            if not cases_eids and not controls_eids:
                continue

            cases_count = len(cases_eids)
            controls_count = len(controls_eids)

            if cases_count > 0 and controls_count > 0:
                cohort_flag = 'case,control'
            elif cases_count > 0:
                cohort_flag = 'case'
            else:
                cohort_flag = 'control'

            def annotation_value(column):
                value = variant_row.get(column, 'NA')
                return 'NA' if pd.isna(value) else str(value)

            high_risk_queue.put(
                f"{gene}\t{annotation_value('Variant')}\t{annotation_value('HGVS.c')}\t"
                f"{annotation_value('HGVS.p')}\t{annotation_value('gnomad_exome_af')}\t"
                f"{annotation_value('gnomad_genome_af')}\t"
                f"{annotation_value('alpha_missense_score')}\t{annotation_value('revel_score')}\t"
                f"{cohort_flag}\t{cases_count}\t{controls_count}\t"
                f"{','.join(cases_eids)}\t{','.join(controls_eids)}\n"
            )

        return result

    except Exception as exc:
        logging.exception("Error processing gene %s: %s", gene, exc)
        return handle_no_variants(gene, 0, 0, 0, 0)


def process_gene_with_progress(args):
    (
        gene,
        variant_file,
        cases_file,
        controls_file,
        pcs_file,
        include_htn,
        summary_queue,
        high_risk_queue,
        progress_queue,
        exclusion_eids,
    ) = args

    result = process_gene(
        (
            gene,
            variant_file,
            cases_file,
            controls_file,
            pcs_file,
            include_htn,
            summary_queue,
            high_risk_queue,
            exclusion_eids,
        )
    )
    progress_queue.put(1)
    return result


def resolve_path(explicit_path, data_dir, relative_filename):
    """Use an explicit path when provided, otherwise resolve under data_dir."""
    if explicit_path:
        return explicit_path
    return os.path.join(data_dir, relative_filename)


def validate_input_file(path, description):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{description} not found: {path}\n"
            "Use the corresponding command-line override if your UKB RAP project "
            "uses a different directory layout."
        )


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run the manuscript UK Biobank gene-based Firth logistic regression "
            "analysis using one pre-filtered qualifying-variant file."
        )
    )

    # Positional arguments.
    parser.add_argument('label', help='Prefix used for output files')
    parser.add_argument('variant_file', help='Pre-filtered qualifying-variant TSV')
    parser.add_argument(
        'phenotype',
        choices=sorted(CASE_FILENAMES),
        help=(
            'A=Dissection, B=TAA, C=TAA-surgery, '
            'D=Dissection+TAA-surgery'
        ),
    )

    parser.add_argument(
        '--no_htn',
        action='store_true',
        help=(
            'Exclude HTN from the model. By default HTN is included.'
        ),
    )
    parser.add_argument(
        '--exclude_eids',
        action='store_true',
        help='Exclude EIDs from the relatedness exclusion file',
    )

    parser.add_argument(
        '--data-dir',
        default='/mnt/project/scripts/lr_data',
        help=(
            'Root directory containing manuscript UKB phenotype/covariate files '
            '(default: /mnt/project/scripts/lr_data)'
        ),
    )
    parser.add_argument(
        '--cases-file',
        help='Override the phenotype-specific case file',
    )
    parser.add_argument(
        '--controls-file',
        help='Override the manuscript control file',
    )
    parser.add_argument(
        '--pcs-file',
        help='Override the manuscript PC file (must contain PC_1 through PC_4)',
    )
    parser.add_argument(
        '--exclusion-file',
        help='Override the relatedness exclusion file used with --exclude_eids',
    )
    parser.add_argument(
        '--output-dir',
        default='.',
        help='Directory for result, summary, variant-detail, and log files',
    )

    return parser


def main():
    args = build_parser().parse_args()

    include_htn = not args.no_htn

    os.makedirs(args.output_dir, exist_ok=True)

    cases_file = resolve_path(
        args.cases_file,
        args.data_dir,
        CASE_FILENAMES[args.phenotype],
    )
    controls_file = resolve_path(
        args.controls_file,
        args.data_dir,
        DEFAULT_CONTROLS_FILENAME,
    )
    pcs_file = resolve_path(
        args.pcs_file,
        args.data_dir,
        DEFAULT_PCS_FILENAME,
    )
    exclusion_file = resolve_path(
        args.exclusion_file,
        args.data_dir,
        DEFAULT_EXCLUSION_FILENAME,
    )

    validate_input_file(args.variant_file, 'Variant file')
    validate_input_file(cases_file, 'Cases file')
    validate_input_file(controls_file, 'Controls file')
    validate_input_file(pcs_file, 'PC file')
    if args.exclude_eids:
        validate_input_file(exclusion_file, 'Exclusion file')

    prefix_name = f"{args.label}_{PHENOTYPE_PREFIX[args.phenotype]}"
    setup_logging(include_htn, prefix_name, args.output_dir)

    logging.info("Phenotype code: %s", args.phenotype)
    logging.info("Cases file: %s", os.path.abspath(cases_file))
    logging.info("Controls file: %s", os.path.abspath(controls_file))
    logging.info("PC file: %s", os.path.abspath(pcs_file))
    logging.info("Variant file: %s", os.path.abspath(args.variant_file))
    logging.info("HTN included: %s", include_htn)

    exclusion_eids = None
    if args.exclude_eids:
        try:
            exclusion_eids = load_exclusion_eids(exclusion_file)
        except Exception as exc:
            print(f"Error loading exclusion EIDs file: {exc}", file=sys.stderr)
            return 1

        print(
            f"Loaded {len(exclusion_eids)} EIDs to exclude from analysis "
            f"from {exclusion_file}"
        )
        logging.info(
            "Loaded %s EIDs to exclude from analysis from %s",
            len(exclusion_eids),
            exclusion_file,
        )

    variant_data = pd.read_csv(args.variant_file, sep='\t', dtype={'EIDs': str})
    variant_data.columns = variant_data.columns.str.strip()
    require_columns(variant_data, ['Gene', 'EIDs'], 'variant file')
    unique_genes = variant_data['Gene'].dropna().unique()

    if len(unique_genes) == 0:
        print("No genes were found in the variant file.", file=sys.stderr)
        return 1

    output_filename = f"{prefix_name}_log_regression_results_firths_PC"
    summary_filename = f"{prefix_name}_log_regression_analysis_firths_PC"
    high_risk_variants_filename = f"{prefix_name}_high_risk_variants_details"

    if include_htn:
        output_filename += "_HTN"
        summary_filename += "_HTN"
        high_risk_variants_filename += "_HTN"

    output_filename = os.path.join(args.output_dir, output_filename + '.txt')
    summary_filename = os.path.join(args.output_dir, summary_filename + '.txt')
    high_risk_variants_filename = os.path.join(
        args.output_dir,
        high_risk_variants_filename + '.txt',
    )

    manager = Manager()
    summary_queue = manager.Queue()
    high_risk_queue = manager.Queue()
    progress_queue = manager.Queue()

    args_list = [
        (
            gene,
            args.variant_file,
            cases_file,
            controls_file,
            pcs_file,
            include_htn,
            summary_queue,
            high_risk_queue,
            progress_queue,
            exclusion_eids,
        )
        for gene in unique_genes
    ]

    num_processes = min(cpu_count(), len(unique_genes))
    results = []

    with Pool(processes=num_processes) as pool:
        with tqdm(
            total=len(unique_genes),
            desc=f"Processing genes ({num_processes} processes)",
        ) as pbar:
            result_async = pool.map_async(process_gene_with_progress, args_list)

            while not result_async.ready():
                for _ in range(progress_queue.qsize()):
                    progress_queue.get()
                    pbar.update()
                result_async.wait(0.1)

            pbar.update(pbar.total - pbar.n)
        results = result_async.get()

    results.sort(
        key=lambda row: (
            float(row.split('\t')[-1])
            if 'NA' not in row.split('\t')[-1]
            else float('inf')
        )
    )

    with open(output_filename, 'w') as f_out:
        f_out.write(
            "Gene\tCases_with_variant\tCases_without_variant\t"
            "Controls_with_variant\tControls_without_variant\tOdds_Ratio\t"
            "95%_CI_Lower\t95%_CI_Upper\tP_value\n"
        )
        f_out.writelines(results)

    with open(summary_filename, 'w') as f_summary:
        while not summary_queue.empty():
            f_summary.write(summary_queue.get())

    with open(high_risk_variants_filename, 'w') as f_high_risk:
        f_high_risk.write(
            "Gene\tVariant\tHGVS.c\tHGVS.p\tgnomad_exome_af\t"
            "gnomad_genome_af\talpha_missense_score\trevel_score\t"
            "Cohort_Flag\tCases_Count\tControls_Count\tCases_EIDs\tControls_EIDs\n"
        )
        while not high_risk_queue.empty():
            f_high_risk.write(high_risk_queue.get())

    print(f"\nFirth logistic regression results saved to '{output_filename}'.")
    print(f"Firth logistic regression summaries saved to '{summary_filename}'.")
    print(
        "Case and control variant details saved to "
        f"'{high_risk_variants_filename}'."
    )

    return 0


if __name__ == '__main__':
    sys.exit(main())
