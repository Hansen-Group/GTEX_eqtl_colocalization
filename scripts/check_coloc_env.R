#!/usr/bin/env Rscript

parse_args <- function(args) {
  parsed <- list()
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
    }

    i <- i + 1
  }

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

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/check_coloc_env.R --lib-path /path/to/R/library --output coloc_env_versions.tsv",
      "",
      "By default this script installs missing CRAN/Bioconductor packages, then writes",
      "the R/package version table. Use --check-only to fail instead of installing.",
      "",
      "Optional:",
      "  --lib-path       Optional .libPaths() entry for a project renv/library",
      "  --output         Version table path [coloc_env_versions.tsv]",
      "  --function-file  Path to R/coloc_functions.R",
      "  --check-only     Do not install missing packages",
      sep = "\n"
    )
  )
}

install_missing_coloc_packages <- function(required_packages, optional_packages) {
  packages <- unique(c(required_packages, optional_packages))
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) == 0) {
    return(invisible(character()))
  }

  bioc_packages <- c(
    "BSgenome",
    "GenomeInfoDb",
    "GenomicRanges",
    "IRanges",
    "S4Vectors",
    "SNPlocs.Hsapiens.dbSNP155.GRCh38"
  )
  cran_packages <- setdiff(missing_packages, bioc_packages)
  missing_bioc_packages <- intersect(missing_packages, bioc_packages)

  if (length(cran_packages) > 0) {
    install.packages(cran_packages, repos = "https://cloud.r-project.org")
  }

  if (length(missing_bioc_packages) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_bioc_packages, ask = FALSE, update = FALSE)
  }

  invisible(missing_packages)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (as_flag(args$help)) {
  usage()
  quit(status = 0)
}

lib_path <- arg_or_default(args, "lib-path", Sys.getenv("COLOC_R_LIB_PATH", unset = ""))
if (!identical(lib_path, "")) {
  dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
  .libPaths(lib_path)
}

script_dir <- get_script_dir()
default_function_file <- normalizePath(file.path(script_dir, "..", "R", "coloc_functions.R"), mustWork = FALSE)
function_file <- arg_or_default(args, "function-file", default_function_file)
source(function_file)

output_file <- arg_or_default(args, "output", "coloc_env_versions.tsv")
if (!as_flag(args$`check-only`)) {
  install_missing_coloc_packages(required_coloc_packages, optional_coloc_packages)
}

env_info <- check_coloc_env(output_file = output_file)
print(env_info)
cat("Environment version table written to:", output_file, "\n")
