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

as_flag <- function(value) {
  if (is.null(value)) {
    return(FALSE)
  }
  if (is.logical(value)) {
    return(isTRUE(value))
  }
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/run_coloc_gtex_v10.R --tissue Liver --project-dir /path/to/project --input result/all_genes_within_range_with_ensemblID.csv",
      "",
      "Required:",
      "  --tissue       GTEx v10 tissue name, e.g. Liver",
      "  --project-dir  Project root directory",
      "  --input        Region/gene input table. Relative paths are resolved under project-dir.",
      "",
      "Optional:",
      "  --output-dir      Output directory [project-dir/result/coloc_pqtl_eqtl_gtex_v10_w_region]",
      "  --eqtl-dir        GTEx v10 allpairs parquet directory",
      "  --sdY             eQTL sdY for coloc quant dataset [1]",
      "  --input-format    auto, region, or lead [auto]",
      "  --dist            Window size around lead SNP for lead input [500000]",
      "  --function-file   Path to R/coloc_functions.R",
      "  --lib-path        Optional .libPaths() entry for a project renv/library",
      "  --env-check-only  Check R/package versions and exit",
      "",
      "Backward compatible positional use:",
      "  Rscript scripts/run_coloc_gtex_v10.R Liver",
      sep = "\n"
    )
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (as_flag(args$help)) {
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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

eqtl_dir <- arg_or_default(
  args,
  "eqtl-dir",
  Sys.getenv(
    "GTEX_V10_EQTL_DIR",
    unset = "/datasets/cbmr_shared/resources/gtex/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/"
  )
)

sdY <- as.numeric(arg_or_default(args, "sdY", Sys.getenv("COLOC_EQTL_SDY", unset = "1")))
input_format <- arg_or_default(args, "input-format", Sys.getenv("COLOC_INPUT_FORMAT", unset = "auto"))
dist <- as.integer(arg_or_default(args, "dist", Sys.getenv("COLOC_LEAD_DIST", unset = "500000")))

tissue <- arg_or_default(args, "tissue", NULL)
if (is.null(tissue) && length(args$positional) >= 1) {
  tissue <- args$positional[[1]]
}

env_file <- file.path(output_dir, "coloc_env_versions.tsv")
env_info <- check_coloc_env(output_file = env_file)
print(env_info)
load_coloc_packages()

if (as_flag(args$`env-check-only`)) {
  cat("Environment check completed. Version table written to:", env_file, "\n")
  quit(status = 0)
}

if (is.null(tissue) || identical(tissue, "")) {
  usage()
  stop("Missing required --tissue argument.", call. = FALSE)
}

cat("Project dir:", project_dir, "\n")
cat("Input table:", input_path, "\n")
cat("Output dir:", output_dir, "\n")
cat("GTEx v10 eQTL dir:", eqtl_dir, "\n")
cat("Processing eQTL tissue:", tissue, "\n")
cat("Using eQTL sdY:", sdY, "\n")
cat("Input format:", input_format, "\n")
cat("Lead SNP window distance:", dist, "\n")

raw_input <- read_coloc_table(input_path)
region_input <- standardize_coloc_input(raw_input, input_format = input_format, dist = dist)

row_names <- region_input$row_label

results <- lapply(
  seq_len(nrow(region_input)),
  function(i) {
    row_df <- region_input[i, , drop = FALSE]

    tryCatch(
      {
        message("Running row ", i, ": ", row_df$result_name, " | region: ", row_df$result_region)
        run_coloc(
          pqtl_info = row_df,
          eQTL_tissue = tissue,
          coloc_result_dir = output_dir,
          eqtl_dir = eqtl_dir,
          sdY = sdY
        )
      },
      error = function(e) {
        message("[ERROR] row ", i, ": ", conditionMessage(e))
        NULL
      }
    )
  }
)
names(results) <- row_names

result_rds <- file.path(output_dir, paste0(tissue, "_all_coloc_results.rds"))
saveRDS(results, result_rds)
cat("Saved per-tissue coloc result list:", result_rds, "\n")

cat(
  "Run scripts/combine_coloc_summaries.R once all tissues have finished",
  "to build the combined summary tables.\n"
)
