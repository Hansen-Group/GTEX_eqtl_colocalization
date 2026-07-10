#!/usr/bin/env Rscript

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
}

parse_args <- function(args) {
  parsed <- list()
  positional <- character()
  i <- 1

  while (i <= length(args)) {
    arg <- args[[i]]

    if (startsWith(arg, "--")) {
      key_value <- sub("^--", "", arg)
      if (grepl("=", key_value, fixed = TRUE)) {
        pieces <- strsplit(key_value, "=", fixed = TRUE)[[1]]
        parsed[[pieces[[1]]]] <- paste(pieces[-1], collapse = "=")
      } else {
        key <- key_value
        next_arg <- if (i < length(args)) args[[i + 1]] else NA_character_
        if (!is.na(next_arg) && !startsWith(next_arg, "--")) {
          parsed[[key]] <- next_arg
          i <- i + 1
        } else {
          parsed[[key]] <- TRUE
        }
      }
    } else {
      positional <- c(positional, arg)
    }

    i <- i + 1
  }

  parsed$positional <- positional
  parsed
}

arg_or_default <- function(parsed, name, default = NULL) {
  value <- parsed[[name]]
  if (is.null(value) || identical(value, "")) {
    return(default)
  }
  value
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/combine_coloc_summaries.R --project-dir /path/to/project --input result/all_genes_within_range_with_ensemblID.csv --output-dir result/coloc_pqtl_eqtl_gtex_v10_w_region",
      "",
      "Combines every *_coloc_results_summary.tsv file in --output-dir into:",
      "  coloc_all_summary_w_region.tsv/.rds",
      "  coloc_sig_pp<threshold>_summary_w_region.tsv/.rds",
      "",
      "Run this once, after every SLURM array task in submit_coloc_gtex_v10.sh has finished.",
      "Do not run it from inside a per-tissue task: multiple array tasks writing this file",
      "concurrently can race and silently drop tissues from the combined output.",
      "",
      "Optional:",
      "  --output-dir      Output directory [project-dir/result/coloc_pqtl_eqtl_gtex_v10_w_region]",
      "  --h4-threshold    PP.H4.abf threshold for significant RDS/TSV [0.7]",
      "  --input-format    auto, region, or lead [auto]",
      "  --dist            Window size around lead SNP for lead input [500000]",
      "  --function-file   Path to R/coloc_functions.R",
      "  --lib-path        Optional .libPaths() entry for a project renv/library",
      sep = "\n"
    )
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!is.null(args$help)) {
  usage()
  quit(status = 0)
}

lib_path <- arg_or_default(args, "lib-path", Sys.getenv("COLOC_R_LIB_PATH", unset = ""))
if (!identical(lib_path, "") && dir.exists(lib_path)) {
  .libPaths(lib_path)
}

script_dir <- get_script_dir()
default_function_file <- normalizePath(file.path(script_dir, "..", "R", "coloc_functions.R"), mustWork = FALSE)
function_file <- arg_or_default(args, "function-file", default_function_file)
source(function_file)
load_coloc_packages()

project_dir <- arg_or_default(
  args,
  "project-dir",
  Sys.getenv("COLOC_PROJECT_DIR", unset = getwd())
)
project_dir <- normalizePath(project_dir, mustWork = FALSE)

input_path <- arg_or_default(
  args,
  "input",
  Sys.getenv("COLOC_INPUT", unset = file.path(project_dir, "result", "all_genes_within_range_with_ensemblID.csv"))
)
if (!grepl("^/", input_path)) {
  input_path <- file.path(project_dir, input_path)
}

output_dir <- arg_or_default(
  args,
  "output-dir",
  Sys.getenv("COLOC_OUTPUT_DIR", unset = file.path(project_dir, "result", "coloc_pqtl_eqtl_gtex_v10_w_region"))
)
if (!grepl("^/", output_dir)) {
  output_dir <- file.path(project_dir, output_dir)
}

h4_threshold <- as.numeric(arg_or_default(args, "h4-threshold", Sys.getenv("COLOC_H4_THRESHOLD", unset = "0.7")))
input_format <- arg_or_default(args, "input-format", Sys.getenv("COLOC_INPUT_FORMAT", unset = "auto"))
dist <- as.integer(arg_or_default(args, "dist", Sys.getenv("COLOC_LEAD_DIST", unset = "500000")))

cat("Project dir:", project_dir, "\n")
cat("Input table:", input_path, "\n")
cat("Output dir:", output_dir, "\n")

raw_input <- read_coloc_table(input_path)
region_input <- standardize_coloc_input(raw_input, input_format = input_format, dist = dist)

combined_summary <- combine_coloc_summaries(output_dir)
if (nrow(combined_summary) == 0) {
  cat("No coloc summary files were found in", output_dir, "\n")
  quit(status = 0)
}

combined_summary <- combined_summary |>
  dplyr::left_join(
    region_input |>
      dplyr::select(result_name, result_region, range_bp, CHR, start, end) |>
      dplyr::distinct(),
    by = c("gene_name" = "result_name", "region" = "result_region")
  )

all_summary_tsv <- file.path(output_dir, "coloc_all_summary_w_region.tsv")
all_summary_rds <- file.path(output_dir, "coloc_all_summary_w_region.rds")
data.table::fwrite(combined_summary, all_summary_tsv, sep = "\t")
saveRDS(combined_summary, all_summary_rds)

sig_summary <- combined_summary |>
  dplyr::filter(PP.H4.abf >= h4_threshold)
sig_tsv <- file.path(output_dir, paste0("coloc_sig_pp", gsub("\\.", "", h4_threshold), "_summary_w_region.tsv"))
sig_rds <- file.path(output_dir, paste0("coloc_sig_pp", gsub("\\.", "", h4_threshold), "_summary_w_region.rds"))
data.table::fwrite(sig_summary, sig_tsv, sep = "\t")
saveRDS(sig_summary, sig_rds)

cat("Combined", nrow(combined_summary), "rows across",
  length(list.files(output_dir, pattern = "_coloc_results_summary.tsv$")),
  "summary files\n")
cat("Saved combined summary:", all_summary_tsv, "\n")
cat("Saved significant summary:", sig_tsv, "\n")
