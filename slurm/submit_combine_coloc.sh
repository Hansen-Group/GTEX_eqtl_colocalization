#!/bin/bash

#SBATCH --job-name=coloc_combine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/coloc_combine_%j.out
#SBATCH --error=logs/coloc_combine_%j.err

set -euo pipefail

# CONST: keep these in sync with slurm/submit_coloc_gtex_v10.sh
PROJECT_DIR="target_project_pth"
REGION_INPUT="${PROJECT_DIR}/result/all_genes_within_range_with_ensemblID.csv"
OUTPUT_DIR="${PROJECT_DIR}/result/coloc_pqtl_eqtl_gtex_v10_w_region"
SUBMIT_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${SUBMIT_SCRIPT_DIR}/.." && pwd)"
FUNCTION_FILE="${CODE_DIR}/R/coloc_functions.R"
RUN_SCRIPT="${CODE_DIR}/scripts/combine_coloc_summaries.R"

# Optional: point to a project renv/library if needed.
COLOC_R_LIB_PATH="/home/wkq953/segment/pipeline/multiome-pipeline_sh_29012025/renv/library/R-4.3/x86_64-pc-linux-gnu"

mkdir -p logs

module load --auto R/4.3.3

echo "Combining coloc summaries in: ${OUTPUT_DIR}"
echo "Started at: $(date)"

Rscript "${RUN_SCRIPT}" \
  --project-dir "${PROJECT_DIR}" \
  --input "${REGION_INPUT}" \
  --output-dir "${OUTPUT_DIR}" \
  --function-file "${FUNCTION_FILE}" \
  --lib-path "${COLOC_R_LIB_PATH}" \
  --input-format auto \
  --dist 500000 \
  --h4-threshold 0.7

echo "Finished at: $(date)"
