# GTEx v10 eQTL colocalization workflow

## 1. Aim and Description

This project runs GTEx v10 eQTL colocalization on the CBMR Esrum HPC. Jobs are submitted with SLURM, with one GTEx tissue per array task, followed by one combine job after all tissues finish.

The current code supports:

- pQTL-eQTL coloc: protein/QTL summary statistics against GTEx eQTL.
- GWAS-eQTL coloc: GWAS summary statistics can be used through the same `file_path` input, as long as the file is formatted with the required association columns.

Main files:

- `R/coloc_functions.R`: shared input parsing, GTEx formatting, allele matching, coloc, and output helpers.
- `scripts/run_coloc_gtex_v10.R`: runs one GTEx tissue.
- `scripts/combine_coloc_summaries.R`: combines per-tissue summaries after the array finishes.
- `scripts/check_coloc_env.R`: checks or installs required R packages and writes package versions.
- `slurm/submit_coloc_gtex_v10.sh`: SLURM array submission.
- `slurm/submit_combine_coloc.sh`: SLURM combine submission.

### 1.1 GTEx eQTL Data

GTEx v10 eQTL all-pairs data comes from:

https://gtexportal.org/home/downloads/adult-gtex/qtl

On Esrum, the data has already been downloaded here:

```text
/datasets/cbmr_shared/resources/gtex/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/
```

The expected file naming pattern is:

```text
<GTEx_tissue>.v10.allpairs.chr<CHR>.parquet
```

Example:

```text
Liver.v10.allpairs.chr1.parquet
```

To use a custom GTEx/eQTL download, change `--eqtl-dir` in the R command or `GTEX_EQTL_DIR` in `slurm/submit_coloc_gtex_v10.sh`.

### 1.2 Environment Setup

Recommended Esrum module:

```bash
module load --auto R/4.3.3
```

Required R packages:

```text
arrow
BSgenome
coloc
data.table
dplyr
GenomeInfoDb
GenomicRanges
IRanges
S4Vectors
SNPlocs.Hsapiens.dbSNP155.GRCh38
stringr
tibble
```

Optional packages, only needed for `get_ensemblID()`:

```text
httr
jsonlite
```

Run the environment check/install script:

```bash
Rscript scripts/check_coloc_env.R \
  --lib-path /home/wkq953/segment/pipeline/multiome-pipeline_sh_29012025/renv/library/R-4.3/x86_64-pc-linux-gnu \
  --output coloc_env_versions.tsv
```

By default, `scripts/check_coloc_env.R` installs missing CRAN/Bioconductor packages and then writes the R/package version table. Use `--check-only` to fail instead of installing.

The current Esrum library path used by the SLURM scripts is:

```text
/home/wkq953/segment/pipeline/multiome-pipeline_sh_29012025/renv/library/R-4.3/x86_64-pc-linux-gnu
```

After connecting to Esrum, run the script above against the library used for this project and keep `coloc_env_versions.tsv` with the project outputs so exact package versions are recorded.

### 1.3 Notes

- Coordinates are assumed to be GRCh38 because GTEx eQTL coordinates use hg38.
- GTEx variant IDs are converted to `chr`, `position`, `effect_allele`, and `other_allele`.
- rsIDs for GTEx variants are annotated with `SNPlocs.Hsapiens.dbSNP155.GRCh38`.
- Rows with missing rsID, beta, varbeta, p-value, or allele frequency are removed before coloc.
- GTEx eQTL is passed to `coloc.abf()` as a quantitative dataset with `sdY = 1` by default. This follows the coloc discussion in [chr1swallace/coloc#201](https://github.com/chr1swallace/coloc/issues/201), which asks whether "`sdY = 1` be used for GTEx eQTL data in coloc". Override with `--sdY` only if you intentionally need a different eQTL phenotype standard deviation.
- `match_coloc_region()` aligns pQTL/GWAS and eQTL effect alleles by literal match or full flip, with no separate check for palindromic A/T or C/G SNPs. This relies on both inputs reporting effect alleles against the same GRCh38 reference/`+` strand. GTEx does this through `variant_id`; check custom pQTL/GWAS strand conventions before trusting coloc results near palindromic variants.

## 2. How to Run

Run one tissue:

```bash
module load --auto R/4.3.3

Rscript scripts/run_coloc_gtex_v10.R \
  --tissue Liver \
  --project-dir /projects/cbmr_shared/people/wkq953/non-GDPR/project_share/Saliva/coloc \
  --input result/all_genes_within_range_with_ensemblID.csv \
  --output-dir result/coloc_pqtl_eqtl_gtex_v10_w_region \
  --eqtl-dir /datasets/cbmr_shared/resources/gtex/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations \
  --lib-path /home/wkq953/segment/pipeline/multiome-pipeline_sh_29012025/renv/library/R-4.3/x86_64-pc-linux-gnu \
  --input-format auto \
  --dist 500000
```

`--input-format auto` detects either fixed-region input or lead-SNP input.

Fixed-region input requires:

```text
range_bp, CHR, start, end, Gene, ensembl, file_path, N
```

Example file:

```text
examples/region_input_example.csv
```

Lead-SNP input requires:

```text
exposure, lead_cispQTL, chr, pos, PG.Genes, file_path, N
```

and one of:

```text
ensemblID, ensembl_gene_id, ensembl
```

Example file:

```text
examples/lead_input_example.tsv
```

For lead-SNP input, the tested region is:

```text
start = max(1, pos - dist)
end   = pos + dist
```

The association file pointed to by `file_path` is used as the pQTL or GWAS dataset. Required columns are:

```text
chr, position, rsid, effect_allele, other_allele, beta, varbeta, pvalue, maf
```

Example file:

```text
examples/association_sumstats_example.tsv
```

Accepted alternative columns:

- `chromosome_hg38` + `base_pair_location_hg38`: converted to `chr`, `position`, and annotated with rsID.
- `p_value`: renamed to `pvalue`.
- `MAF`: renamed to `maf`.
- `eaf`: converted to `maf = ifelse(eaf > 0.5, 1 - eaf, eaf)` if `maf` is missing.
- `se`, `SE`, `StdErr`, or `standard_error`: converted to `varbeta = se^2` if `varbeta` is missing.

## 3. SLURM Submission

Edit the constants at the top of `slurm/submit_coloc_gtex_v10.sh`:

```bash
PROJECT_DIR="/projects/cbmr_shared/people/wkq953/non-GDPR/project_share/Saliva/coloc"
REGION_INPUT="${PROJECT_DIR}/result/all_genes_within_range_with_ensemblID.csv"
GTEX_TISSUE_TABLE="/projects/holbaek-AUDIT/people/wkq953/project_proteomics_holbaek/data/eqtl_gtex_v10_hg38/gtex_eqtl_sample_size_formatted.csv"
GTEX_EQTL_DIR="/datasets/cbmr_shared/resources/gtex/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations"
OUTPUT_DIR="${PROJECT_DIR}/result/coloc_pqtl_eqtl_gtex_v10_w_region"
COLOC_R_LIB_PATH="/home/wkq953/segment/pipeline/multiome-pipeline_sh_29012025/renv/library/R-4.3/x86_64-pc-linux-gnu"
```

Submit from the repo root:

```bash
sbatch slurm/submit_coloc_gtex_v10.sh
```

The default array is `#SBATCH --array=2-51`, assuming line 1 of `GTEX_TISSUE_TABLE` is a header and column 1 is the GTEx tissue name. Logs are written to:

```text
logs/coloc_<tissue>.out
logs/coloc_<tissue>.err
```

Combine results once all array jobs finish:

```bash
sbatch slurm/submit_combine_coloc.sh
```

Or chain the combine job:

```bash
job_id=$(sbatch --parsable slurm/submit_coloc_gtex_v10.sh)
sbatch --dependency=afterok:${job_id} slurm/submit_combine_coloc.sh
```

Do not run the combine script from inside each tissue task; parallel combine jobs can race while writing the same summary files.

## 4. Output

For each successful gene-region-tissue coloc:

```text
<gene>_<region>_<tissue>_coloc_results_summary.tsv
<gene>_<region>_<tissue>_coloc_results_full.tsv
```

Per tissue:

```text
<tissue>_all_coloc_results.rds
```

After running `scripts/combine_coloc_summaries.R` or `slurm/submit_combine_coloc.sh`:

```text
coloc_all_summary_w_region.tsv
coloc_all_summary_w_region.rds
coloc_sig_pp07_summary_w_region.tsv
coloc_sig_pp07_summary_w_region.rds
```

Important output columns:

- `gene_name`, `gene_symbol`, `ensembl`, `region`, `tissue`, `input_format`
- `nsnps`
- `PP.H0.abf`, `PP.H1.abf`, `PP.H2.abf`, `PP.H3.abf`, `PP.H4.abf`
- `pp4_cond = PP.H4.abf / (PP.H3.abf + PP.H4.abf)`
- `CHR`, `start`, `end`

Use `readRDS()` for `.rds` files, not `data.table::fread()`.
