#!/bin/bash
#SBATCH --job-name=t2d_prs
#SBATCH --array=1-154
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=t2d_prs_%A_%a.out
#SBATCH --error=t2d_prs_%A_%a.err

set -euo pipefail

# ============================================================
# T2D chromosome × P-value-threshold PRS scoring
#
# 22 chromosomes × 7 thresholds = 154 array tasks.
#
# This reproduces the submitted threshold-only PRS workflow.
# No LD clumping is performed.
# ============================================================

# ------------------------------------------------------------
# Repository root
# ------------------------------------------------------------

SCRIPT_DIR="$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)"

REPO_DIR="$(
  cd "${SCRIPT_DIR}/.."
  pwd
)"

cd "${REPO_DIR}"

# ------------------------------------------------------------
# Arguments / configuration
# ------------------------------------------------------------

# Directory containing:
#
#   Testing_Chr1.pgen
#   Testing_Chr1.pvar
#   Testing_Chr1.psam
#   ...
#   Testing_Chr22.*
#
GENO_DIR="${1:-}"

if [[ -z "${GENO_DIR}" ]]; then
    echo "Usage:"
    echo "  sbatch scripts/02_generate_prs.sbatch /path/to/plink_files"
    exit 1
fi

SCORE_DIR="outputs/sumstats/t2d/score_files"
OUT_DIR="outputs/prs/t2d"

PLINK2_BIN="${PLINK2_BIN:-plink2}"

# Historical filename labels:
#
# p1e5  = P <= 1e-5
# p1e3  = P <= 1e-3
# p0.01 = P <= 0.01
# ...
# p1    = P <= 1
THRESHOLD_LABELS=(
  "1e5"
  "1e3"
  "0.01"
  "0.05"
  "0.1"
  "0.5"
  "1"
)

# ------------------------------------------------------------
# Check PLINK
# ------------------------------------------------------------

if ! command -v "${PLINK2_BIN}" >/dev/null 2>&1; then
    echo "ERROR: plink2 not found in PATH." >&2
    echo "Load or install PLINK 2 before submitting this job." >&2
    exit 1
fi

echo "PLINK version:"
"${PLINK2_BIN}" --version

# ------------------------------------------------------------
# Convert array index -> chromosome + threshold
# ------------------------------------------------------------

N_THRESHOLDS="${#THRESHOLD_LABELS[@]}"
EXPECTED_TASKS=$((22 * N_THRESHOLDS))

if (( SLURM_ARRAY_TASK_ID < 1 ||
      SLURM_ARRAY_TASK_ID > EXPECTED_TASKS )); then

    echo \
      "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and ${EXPECTED_TASKS}." \
      >&2

    exit 1
fi

TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

CHR=$(
  (
    TASK_ID / N_THRESHOLDS
  ) + 1
)

THRESHOLD="${
  THRESHOLD_LABELS[
    TASK_ID % N_THRESHOLDS
  ]
}"

# ------------------------------------------------------------
# Input / output files
# ------------------------------------------------------------

PFILE="${GENO_DIR}/Testing_Chr${CHR}"

SCORE_FILE="${
  SCORE_DIR
}/t2d_chr${CHR}_p${THRESHOLD}.score.tsv"

OUT_PREFIX="${
  OUT_DIR
}/t2d_chr${CHR}_p${THRESHOLD}_prs"

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------
# Input validation
# ------------------------------------------------------------

echo "Chromosome: ${CHR}"
echo "Threshold:  ${THRESHOLD}"
echo "Genotype:   ${PFILE}"
echo "Score file: ${SCORE_FILE}"
echo "Output:     ${OUT_PREFIX}"

if [[ ! -s "${SCORE_FILE}" ]]; then

    echo \
      "ERROR: score file missing or empty: ${SCORE_FILE}" \
      >&2

    exit 1
fi

for ext in pgen pvar psam; do

    if [[ ! -f "${PFILE}.${ext}" ]]; then

        echo \
          "ERROR: missing genotype file: ${PFILE}.${ext}" \
          >&2

        exit 1
    fi

done

# ------------------------------------------------------------
# PLINK 2 scoring
#
# Score-file columns:
#
#   1 = variant_id
#   2 = effect_allele
#   3 = beta = log(odds_ratio)
#
# The resulting SCORE1_AVG and ALLELE_CT fields are later used
# to reconstruct chromosome-weighted genome-wide PRSs.
# ------------------------------------------------------------

"${PLINK2_BIN}" \
  --pfile "${PFILE}" \
  --score "${SCORE_FILE}" 1 2 3 header \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory 16000 \
  --out "${OUT_PREFIX}"

echo \
  "Finished chromosome ${CHR}, threshold ${THRESHOLD}"
