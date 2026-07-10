#!/bin/bash

#SBATCH --job-name=coloc_gtex_v10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=70:00:00
#SBATCH --array=2-51
# Per-task stdout/stderr are renamed to logs/coloc_<tissue>.out/.err below via `exec`,
# so no static --output/--error here. If a task fails before that `exec` line runs,
# SLURM falls back to its default slurm-%A_%a.out/.err in the submission directory.

set -euo pipefail

# CONST: edit these paths for each project.
PROJECT_DIR="/projects/cbmr_shared/people/wkq953/non-GDPR/project_share/Saliva/coloc"
REGION_INPUT="${PROJECT_DIR}/result/all_genes_within_range_with_ensemblID.csv"
GTEX_EQTL_DIR="/datasets/cbmr_shared/resources/gtex/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations"
OUTPUT_DIR="${PROJECT_DIR}/result/coloc_pqtl_eqtl_gtex_v10_w_region"
SUBMIT_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${SUBMIT_SCRIPT_DIR}/.." && pwd)"
GTEX_TISSUE_TABLE="${CODE_DIR}/data/gtex_eqtl_sample_size_formatted.csv"
FUNCTION_FILE="${CODE_DIR}/R/coloc_functions.R"
RUN_SCRIPT="${CODE_DIR}/scripts/run_coloc_gtex_v10.R"

# Optional: point to a project renv/library if needed.
COLOC_R_LIB_PATH="/home/wkq953/segment/pipeline/multiome-pipeline_sh_29012025/renv/library/R-4.3/x86_64-pc-linux-gnu"

mkdir -p logs "${OUTPUT_DIR}"

# for GTEx v10 tissue format; array starts from line 2 to skip the header.
LINE_NUM="${SLURM_ARRAY_TASK_ID}"
TASK="$(sed -n "${LINE_NUM}p" "${GTEX_TISSUE_TABLE}")"
eqtl_tissue="$(echo "${TASK}" | cut -d',' -f1)"

exec > "logs/coloc_${eqtl_tissue}.out" 2> "logs/coloc_${eqtl_tissue}.err"

module load --auto R/4.3.3

echo "Running coloc for GTEx v10 tissue: ${eqtl_tissue}"
echo "Started at: $(date)"

Rscript "${RUN_SCRIPT}" \
  --tissue "${eqtl_tissue}" \
  --project-dir "${PROJECT_DIR}" \
  --input "${REGION_INPUT}" \
  --output-dir "${OUTPUT_DIR}" \
  --eqtl-dir "${GTEX_EQTL_DIR}" \
  --function-file "${FUNCTION_FILE}" \
  --lib-path "${COLOC_R_LIB_PATH}" \
  --input-format auto \
  --dist 500000

echo "Finished at: $(date)"
