# HTAD Rare-Variant Burden Analysis

Code supporting the primary gene-based rare-variant burden and meta-analysis results reported in:

**Gene-Based Rare Variant Burden Analyses Across Biobanks Identify Novel High-Risk Genes for Thoracic Aortic Disease**

This repository contains the statistical analysis code used for UK Biobank (UKB) and All of Us (AoU) gene-based burden testing, the UKB+AoU discovery meta-analysis, and the subsequent multi-cohort meta-analysis incorporating replication cohorts. The UKB burden analyses were run within the **UK Biobank Research Analysis Platform (RAP), hosted by DNAnexus**, while the AoU analyses were run within the **All of Us Researcher Workbench**. The discovery and multi-cohort meta-analysis scripts are intended to be run **locally** on prepared **summary-level cohort result files** after the cohort-specific burden analyses have been completed; individual-level participant data are not required for these meta-analysis stages.

## Repository structure

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── ukb/
│   └── ukb_firth_logistic_regression.py
├── aou/
│   ├── run_aou_firth.sh
│   └── aou_firth_logistic_regression.py
└── meta_analysis/
    ├── meta_analysis_discovery.py
    └── meta_analysis_multicohort.py
```

## Analysis workflow

```text
Pre-filtered qualifying variants / gene carrier lists
                    │
            ┌───────┴────────┐
            │                │
            ▼                ▼
           UKB              AoU
    UKB RAP / DNAnexus  Researcher Workbench
     Firth regression   Firth regression
            │                │
            └───────┬────────┘
                    ▼
        UKB + AoU discovery meta-analysis
          local execution; fixed-effect IVW
                    │
                    ▼
          Replication cohort results
                    │
                    ▼
          Multi-cohort meta-analysis
             local execution
          ├── fixed-effect IVW (primary)
          └── DerSimonian-Laird random effects (sensitivity)
```

## Scope

These scripts begin with **pre-filtered qualifying-variant or gene-carrier files**. Variant annotation, pathogenicity filtering, allele-frequency filtering, and construction of the carrier lists occur upstream and are not implemented here.

The primary variant classes analyzed are:

- `missense`
- `lof`
- `missense+lof`

The manuscript phenotype codes used in the public analysis scripts are:

| Code | Phenotype |
|---|---|
| `A` | Dissection |
| `B` | TAA |
| `C` | TAA-surgery |
| `D` | Dissection + TAA-surgery |

## 1. UK Biobank burden analysis

The UK Biobank analysis was performed within the **UK Biobank Research Analysis Platform (RAP), hosted by DNAnexus**. The `ukb_firth_logistic_regression.py` script is intended to be run in that environment and assumes access to project-mounted UKB phenotype, covariate, principal-component, relatedness, and pre-filtered variant files. Alternative file locations can be supplied through the command-line options described below.

Script:

```text
ukb/ukb_firth_logistic_regression.py
```

Within the RAP/DNAnexus environment, the UKB analysis uses Firth penalized logistic regression with a binary gene-level qualifying-variant carrier indicator. The model includes age, sex, and PC1-PC4.

The public script preserves the historical command-line behavior in which hypertension (`HTN`) is included by default. To reproduce the **primary manuscript model**, use `--no_htn`. Related individuals are removed when `--exclude_eids` is supplied.

### Example

```bash
python ukb/ukb_firth_logistic_regression.py \
    missense \
    /path/to/prefiltered_missense_variants.tsv \
    A \
    --no_htn \
    --exclude_eids \
    --data-dir /path/to/ukb_analysis_files \
    --output-dir results/ukb
```

The first positional argument is an output label, the second is the pre-filtered qualifying-variant file, and the third is the phenotype code.

The script also accepts explicit overrides for the case, control, principal-component, and relatedness-exclusion files:

```text
--cases-file
--controls-file
--pcs-file
--exclusion-file
```

### Expected UKB variant input

The qualifying-variant file must contain at least:

```text
Gene
EIDs
```

where `EIDs` contains comma-separated UKB participant identifiers for carriers of each variant. Variant annotation columns are used when available for variant-detail output but are not required for the regression itself.

## 2. All of Us burden analysis

Scripts:

```text
aou/run_aou_firth.sh
aou/aou_firth_logistic_regression.py
```

The shell script is the job-submission wrapper for the All of Us Researcher Workbench. It submits chromosome- and phenotype-specific jobs and passes the required files to the Python worker.

The AoU Firth model includes:

```text
carrier + age_at_enrollment + sex_binary + PC1 + PC2 + PC3 + PC4 + PC5
```

Related individuals listed in the supplied exclusion file are removed before model fitting. The manuscript analysis uses `ANCESTRY_PRED=all` in the launcher.

### Required environment variables

Before using the launcher, define:

```bash
export WORKSPACE_BUCKET="gs://your-workspace-bucket"
export FIRTHLOGIST_URI="gs://path/to/firthlogist.zip"
```

The launcher uses the current AoU dsub helper by default:

```text
$HOME/aou_dsub.gcb.bash
```

This can be changed with `AOU_DSUB_HELPER`. The worker-script location can be changed with `AOU_WORKER_SCRIPT`.

### Examples

Run all four phenotypes and all chromosomes for missense variants:

```bash
bash aou/run_aou_firth.sh --variant-class missense
```

Run selected phenotypes:

```bash
bash aou/run_aou_firth.sh \
    --variant-class lof \
    --phenotypes A,D
```

Run selected chromosomes:

```bash
bash aou/run_aou_firth.sh \
    --variant-class missense+lof \
    --chromosomes 1,2,3,X
```

The launcher accepts chromosomes `1-22`, `X`, and `Y`.

### Expected AoU gene-carrier input

Each chromosome-specific carrier file must contain at least:

```text
gene
carrier_samples
```

where `carrier_samples` is a comma-separated list of AoU `person_id` values carrying a qualifying variant in that gene.

## 3. UKB + AoU discovery meta-analysis

Script:

```text
meta_analysis/meta_analysis_discovery.py
```

This script is intended to be run **locally** and performs a two-cohort fixed-effect inverse-variance-weighted meta-analysis of prepared **summary-level UKB and AoU Firth result files**. It does not require individual-level participant data.

For each cohort, the script reconstructs the log odds ratio and standard error from the reported odds ratio and 95% confidence interval:

```text
beta = log(OR)
SE = [log(CI_upper) - log(CI_lower)] / (2 x 1.96)
```

Only genes with valid estimates in both UKB and AoU are included. The output includes the pooled OR, 95% CI, P value, Cochran's Q, and heterogeneity P value.

### Example

```bash
python meta_analysis/meta_analysis_discovery.py \
    --ukb results/ukb_dissection_missense.txt \
    --aou results/aou_dissection_missense.txt \
    --meta-output results/discovery_dissection_missense.tsv \
    --cohort-output results/discovery_dissection_missense_with_cohorts.tsv
```

The UKB and AoU input files are expected to contain:

```text
Gene
Cases_with_variant
Cases_without_variant
Controls_with_variant
Controls_without_variant
Odds_Ratio
95%_CI_Lower
95%_CI_Upper
P_value
```

## 4. Multi-cohort meta-analysis

Script:

```text
meta_analysis/meta_analysis_multicohort.py
```

This script is intended to be run **locally** and performs the overall meta-analysis across prepared **summary-level UKB, AoU, and replication-cohort result files**. It does not require individual-level participant data.

The **fixed-effect inverse-variance-weighted model is the primary analysis**. A **DerSimonian-Laird random-effects model** is generated as a sensitivity analysis.

For the manuscript replication workflow, cohort-gene rows with zero case carriers can be excluded with:

```text
--exclude-zero-case-carriers
```

Overall meta-analysis results are reported only when UKB, AoU, and at least one replication cohort contribute after filtering.

### Minimum input columns

Each cohort-level file must contain:

```text
gene
beta
se
```

A `Cohort` or `cohort` column is strongly recommended. Optional columns used for annotations and carrier counts include:

```text
Phenotype
Variant_Class
p
carriers_case
carriers_ctrl
cases_total
controls_total
```

The conversion of cohort-specific regression results into this standardized `beta`/`se` schema is an upstream formatting step and is not performed by this script.

### Example

```bash
python meta_analysis/meta_analysis_multicohort.py \
    UKB_results.tsv \
    AOU_results.tsv \
    PMBB_results.tsv \
    MGBB_results.tsv \
    MVP_results.tsv \
    --exclude-zero-case-carriers \
    -o results/overall_meta.tsv
```

Both fixed-effect and DerSimonian-Laird random-effects TSV and XLSX outputs are generated.

### Replication support

Independent replication support is defined in the script as at least one replication cohort with:

1. nominal `P < 0.05`, and
2. an effect direction concordant with the fixed-effect UKB+AoU discovery estimate.

The same replication cohort must satisfy both conditions.

### Study-wide significance thresholds

The analysis uses the manuscript thresholds:

| Variant class | Threshold |
|---|---:|
| Missense | `2.6e-6` |
| LOF | `1.3e-5` |
| Missense+LOF | `1.3e-5` |

`Cross_Cohort_Validation` requires both independent replication support and study-wide significance in the overall meta-analysis.

## Software requirements

Python dependencies are listed in `requirements.txt`.

Install with:

```bash
python -m pip install -r requirements.txt
```

The Firth regression analyses use the [`firthlogist`](https://github.com/jzluo/firthlogist) package. In AoU, the worker can either import an installed copy of `firthlogist` or use the ZIP file localized by `run_aou_firth.sh`.

The UKB script is intended for execution within the UK Biobank RAP/DNAnexus environment. The AoU launcher additionally requires access to the All of Us Researcher Workbench job-submission environment. The two meta-analysis scripts do not require either cloud platform and are intended for local execution once the required cohort-level result files have been prepared.

## Data access

This repository does not contain individual-level UK Biobank or All of Us data, participant identifiers, carrier lists, or controlled-access cohort data. Researchers must obtain access to the relevant resources through the corresponding data-access procedures and construct the required input files within the approved computing environments, including the UK Biobank RAP/DNAnexus environment for UKB analyses and the All of Us Researcher Workbench for AoU analyses.

## Reproducibility notes

- The scripts are intended to reproduce the statistical analysis stages described above once the required pre-filtered cohort inputs have been prepared.
- The UKB burden-analysis script is designed for execution within the **UK Biobank Research Analysis Platform (RAP), hosted by DNAnexus**, and uses paths consistent with a RAP project-mounted filesystem unless overridden.
- The AoU burden-analysis scripts are designed for execution within the **All of Us Researcher Workbench**. `run_aou_firth.sh` submits chromosome-, phenotype-, and variant-class-specific cloud jobs, and `aou_firth_logistic_regression.py` performs the corresponding gene-based Firth logistic-regression analyses.
- The discovery and multi-cohort meta-analysis scripts are intended for **local execution using prepared summary-level cohort result files**. Individual-level participant data are not required for these meta-analysis stages.
- Cloud job-submission interfaces may evolve independently of the statistical analysis code. `run_aou_firth.sh` uses the current All of Us dsub helper interface used when this repository was prepared.
- Filenames and directory layouts can be changed using the command-line or environment-variable overrides documented in the scripts.

## Citation

If you use this code, please cite the current preprint:

> Murdock DR, Guan P, Guo D, et al. **Gene-Based Rare Variant Burden Analyses Across Biobanks Identify Novel High-Risk Genes for Thoracic Aortic Disease.** *medRxiv*. 2026. doi: [10.64898/2026.09.16.26363254](https://doi.org/10.64898/2026.09.16.26363254).

## License

See the repository `LICENSE` file.
