#!/usr/bin/env bash
set -euo pipefail

# Submit All of Us gene-based Firth logistic-regression jobs by phenotype,
# chromosome, and variant class.
#
# Required environment variables:
#   WORKSPACE_BUCKET   All of Us workspace bucket (for example, gs://...)
#   FIRTHLOGIST_URI    Cloud URI to firthlogist.zip
#
# Optional environment variables:
#   AOU_DSUB_HELPER    Path to the current AoU dsub helper script
#                      (default: $HOME/aou_dsub.gcb.bash)
#   AOU_DSUB_IMAGE     Container image used for submitted jobs
#   AOU_OUTPUT_DIR     Output directory in the workspace bucket
#   AOU_LOG_DIR        dsub logging directory in the workspace bucket
#   AOU_WORKER_SCRIPT  Cloud path to aou_firth_logistic_regression.py
#
# Usage:
#   ./run_aou_firth.sh --variant-class missense
#   ./run_aou_firth.sh --variant-class lof --phenotypes A,C
#   ./run_aou_firth.sh --variant-class missense+lof --phenotypes D --chromosomes 1,2,3,X
#
# Manuscript phenotype codes:
#   A = Dissection
#   B = TAA
#   C = TAA-surgery
#   D = Dissection + TAA-surgery

VARIANT_CLASS=""
PHENOTYPES="A,B,C,D"
CHROMOSOMES="all"

usage() {
    cat <<'USAGE'
Usage:
  run_aou_firth.sh --variant-class <missense|lof|missense+lof> [options]

Options:
  --variant-class CLASS   Required. missense, lof, or missense+lof
  --phenotypes LIST       Comma-separated phenotype codes (default: A,B,C,D)
  --chromosomes LIST      Comma-separated chromosomes, or 'all' (default: all)
  -h, --help              Show this help message

Phenotypes:
  A  Dissection
  B  TAA
  C  TAA-surgery
  D  Dissection + TAA-surgery

Examples:
  ./run_aou_firth.sh --variant-class missense
  ./run_aou_firth.sh --variant-class lof --phenotypes A,D
  ./run_aou_firth.sh --variant-class missense+lof --chromosomes 1,2,3,X
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --variant-class)
            VARIANT_CLASS="${2:-}"
            shift 2
            ;;
        --phenotypes)
            PHENOTYPES="${2:-}"
            shift 2
            ;;
        --chromosomes)
            CHROMOSOMES="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$VARIANT_CLASS" ]]; then
    echo "ERROR: --variant-class is required." >&2
    usage >&2
    exit 2
fi

case "$VARIANT_CLASS" in
    missense|lof|missense+lof) ;;
    *)
        echo "ERROR: --variant-class must be missense, lof, or missense+lof." >&2
        exit 2
        ;;
esac

: "${WORKSPACE_BUCKET:?Set WORKSPACE_BUCKET before running this script.}"
: "${FIRTHLOGIST_URI:?Set FIRTHLOGIST_URI to the cloud URI for firthlogist.zip.}"

AOU_DSUB_HELPER="${AOU_DSUB_HELPER:-${HOME}/aou_dsub.gcb.bash}"
AOU_DSUB_IMAGE="${AOU_DSUB_IMAGE:-us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14}"
AOU_OUTPUT_DIR="${AOU_OUTPUT_DIR:-${WORKSPACE_BUCKET}/aou_firth_results}"
AOU_LOG_DIR="${AOU_LOG_DIR:-${WORKSPACE_BUCKET}/dsub/logs}"

if [[ ! -f "$AOU_DSUB_HELPER" ]]; then
    echo "ERROR: AoU dsub helper not found: $AOU_DSUB_HELPER" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "$AOU_DSUB_HELPER"

# Input files shared across analyses.
CONTROLS_FILE="${WORKSPACE_BUCKET}/cohorts/non_taad_controls_v3_no_mfs_lds_veds_ids.tsv"
COVARIATES_FILE="${WORKSPACE_BUCKET}/PCA/AOUv8_pca_covariates.ancestry.tsv"
RELATED_FILE="${WORKSPACE_BUCKET}/cohorts/aou_related_to_remove.v4.txt"
AOU_WORKER_SCRIPT="${AOU_WORKER_SCRIPT:-${WORKSPACE_BUCKET}/dsub/aou_firth_logistic_regression.py}"

# Manuscript phenotypes.
declare -A CASE_FILES=(
    [A]="${WORKSPACE_BUCKET}/cohorts/aortic_dissection_v2_ids.tsv"
    [B]="${WORKSPACE_BUCKET}/cohorts/taa_no_rupture_no_dissection_ids_v2.tsv"
    [C]="${WORKSPACE_BUCKET}/cohorts/taa_no_rupture_no_dissection_ids_v2.surgery.tsv"
    [D]="${WORKSPACE_BUCKET}/cohorts/aortic_dissection_v2_taa_no_rupture_open_surgery_ids.tsv"
)

declare -A PHENOTYPE_PREFIXES=(
    [A]="dissection"
    [B]="taa"
    [C]="taa_surgery"
    [D]="dissection_taa_surgery"
)

IFS=',' read -r -a PHENOTYPE_ARRAY <<< "$PHENOTYPES"

if [[ "$CHROMOSOMES" == "all" ]]; then
    CHROMOSOME_ARRAY=({1..22} X Y)
else
    IFS=',' read -r -a CHROMOSOME_ARRAY <<< "$CHROMOSOMES"
fi

for phenotype_code in "${PHENOTYPE_ARRAY[@]}"; do
    phenotype_code="${phenotype_code//[[:space:]]/}"

    if [[ -z "${CASE_FILES[$phenotype_code]+x}" ]]; then
        echo "ERROR: Invalid phenotype code '$phenotype_code'. Use A, B, C, or D." >&2
        exit 2
    fi

    case_file="${CASE_FILES[$phenotype_code]}"
    phenotype_prefix="${PHENOTYPE_PREFIXES[$phenotype_code]}"

    for chrom in "${CHROMOSOME_ARRAY[@]}"; do
        chrom="${chrom//[[:space:]]/}"

        if [[ ! "$chrom" =~ ^([1-9]|1[0-9]|2[0-2]|X|Y)$ ]]; then
            echo "ERROR: Invalid chromosome '$chrom'. Use 1-22, X, or Y." >&2
            exit 2
        fi

        gene_file="${WORKSPACE_BUCKET}/carrier_lists_latest/${VARIANT_CLASS}/chr${chrom}_${VARIANT_CLASS}_gene_carrier_counts.QC.tsv"
        local_output="/tmp/logistic_regression_results_chr${chrom}_${phenotype_prefix}_${VARIANT_CLASS}.txt"
        cloud_output="${AOU_OUTPUT_DIR}/chr${chrom}_${phenotype_prefix}_${VARIANT_CLASS}_firths_LR_AOU.all.txt"

        echo "Submitting ${VARIANT_CLASS}: phenotype ${phenotype_code} (${phenotype_prefix}), chr${chrom}"

        aou_dsub \
            --logging "$AOU_LOG_DIR" \
            --input FIRTHLOGIST_ZIP="$FIRTHLOGIST_URI" \
            --image "$AOU_DSUB_IMAGE" \
            --env GENE_FILE="$gene_file" \
            --env CASES_FILE="$case_file" \
            --env CONTROLS_FILE="$CONTROLS_FILE" \
            --env COVARIATES_FILE="$COVARIATES_FILE" \
            --env RELATED_FILE="$RELATED_FILE" \
            --env OUTPUT_FILE="$local_output" \
            --env ANCESTRY_PRED="all" \
            --output OUTPUT_FILE="$cloud_output" \
            --script "$AOU_WORKER_SCRIPT" \
            --machine-type "e2-standard-32" \
            --boot-disk-size 150
    done
done
