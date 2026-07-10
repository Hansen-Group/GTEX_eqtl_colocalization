required_coloc_packages <- c(
  "arrow",
  "BSgenome",
  "coloc",
  "data.table",
  "dplyr",
  "GenomeInfoDb",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "stringr",
  "tibble"
)

optional_coloc_packages <- c("httr", "jsonlite")

load_coloc_packages <- function(required_packages = required_coloc_packages) {
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

check_coloc_env <- function(
  required_packages = required_coloc_packages,
  optional_packages = optional_coloc_packages,
  output_file = NULL
) {
  packages <- unique(c(required_packages, optional_packages))
  env_info <- data.frame(
    package = packages,
    required = packages %in% required_packages,
    installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    version = NA_character_,
    stringsAsFactors = FALSE
  )

  installed_idx <- env_info$installed
  env_info$version[installed_idx] <- vapply(
    env_info$package[installed_idx],
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  )

  r_info <- data.frame(
    package = "R",
    required = TRUE,
    installed = TRUE,
    version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )

  env_info <- rbind(r_info, env_info)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      env_info,
      file = output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  missing_required <- env_info$package[env_info$required & !env_info$installed]
  if (length(missing_required) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  env_info
}

assert_required_columns <- function(df, required_cols, object_name) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

read_coloc_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  sep <- if (ext == "csv") "," else "\t"
  data.table::fread(path, sep = sep, data.table = FALSE)
}

detect_coloc_input_format <- function(input_df) {
  region_cols <- c("range_bp", "CHR", "start", "end", "Gene", "ensembl", "file_path", "N")
  lead_cols <- c("exposure", "lead_cispQTL", "chr", "pos", "PG.Genes", "file_path", "N")

  if (all(region_cols %in% colnames(input_df))) {
    return("region")
  }
  if (all(lead_cols %in% colnames(input_df)) &&
      any(c("ensemblID", "ensembl_gene_id", "ensembl") %in% colnames(input_df))) {
    return("lead")
  }

  stop(
    "Cannot detect coloc input format. Expected either region columns: ",
    paste(region_cols, collapse = ", "),
    "; or lead-SNP columns: ",
    paste(c(lead_cols, "ensemblID or ensembl_gene_id"), collapse = ", "),
    call. = FALSE
  )
}

standardize_coloc_input <- function(input_df, input_format = "auto", dist = 500000) {
  if (identical(input_format, "auto")) {
    input_format <- detect_coloc_input_format(input_df)
  }

  if (identical(input_format, "region")) {
    required_region_cols <- c("range_bp", "CHR", "start", "end", "Gene", "ensembl", "file_path", "N")
    assert_required_columns(input_df, required_region_cols, "region input")

    input_df |>
      dplyr::mutate(
        input_format = "region",
        CHR = as.character(CHR),
        start = as.integer(start),
        end = as.integer(end),
        Gene = as.character(Gene),
        gene_symbol = as.character(Gene),
        result_name = as.character(Gene),
        result_region = as.character(range_bp),
        row_label = paste(result_name, result_region, sep = "_")
      )
  } else if (identical(input_format, "lead")) {
    required_lead_cols <- c("exposure", "lead_cispQTL", "chr", "pos", "PG.Genes", "file_path", "N")
    assert_required_columns(input_df, required_lead_cols, "lead-SNP input")

    ensembl_col <- intersect(c("ensemblID", "ensembl_gene_id", "ensembl"), colnames(input_df))
    if (length(ensembl_col) == 0) {
      stop("lead-SNP input requires one of: ensemblID, ensembl_gene_id, ensembl", call. = FALSE)
    }
    ensembl_col <- ensembl_col[[1]]

    input_df |>
      dplyr::mutate(
        input_format = "lead",
        range_start = pmax(1L, as.integer(pos) - as.integer(dist)),
        range_end = as.integer(pos) + as.integer(dist),
        range_bp = paste0(range_start, "-", range_end),
        CHR = as.character(chr),
        start = range_start,
        end = range_end,
        Gene = as.character(exposure),
        gene_symbol = as.character(PG.Genes),
        ensembl = as.character(.data[[ensembl_col]]),
        result_name = as.character(exposure),
        result_region = as.character(lead_cispQTL),
        row_label = paste(result_name, result_region, sep = "_")
      )
  } else {
    stop("--input-format must be one of: auto, region, lead", call. = FALSE)
  }
}

get_ensemblID <- function(gene_symbol, entity = "target") {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("get_ensemblID() requires optional packages httr and jsonlite.", call. = FALSE)
  }

  query_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  request_body <- list(
    operationName = "EnsemblIDSearch",
    variables = list(queryString = gene_symbol),
    query = "
      query EnsemblIDSearch($queryString: String!) {
        search(queryString: $queryString) {
          hits {
            id
            name
            entity
          }
        }
      }
    "
  )

  response <- httr::POST(query_url, body = request_body, encode = "json")
  if (httr::http_error(response)) {
    warning(sprintf("Query failed for %s", gene_symbol))
    return(data.frame(gene = gene_symbol, ensembl = NA_character_))
  }

  response_text <- rawToChar(response$content)
  response_data <- jsonlite::fromJSON(response_text)
  hits <- response_data$data$search$hits

  if (is.null(hits) || length(hits) == 0) {
    return(data.frame(gene = gene_symbol, ensembl = NA_character_))
  }

  hits_df <- as.data.frame(hits)
  hits_ensembl_df <- hits_df[hits_df$entity == entity, , drop = FALSE]
  data.frame(gene = gene_symbol, ensembl = hits_ensembl_df$id[1])
}

annotate_snp_with_rsid_hg38 <- function(hg38_data, chromosome = "chr", position = "position") {
  flag_class_gr <- is(hg38_data, "GRanges")

  if (flag_class_gr) {
    hg38_gr <- hg38_data
  } else {
    assert_required_columns(hg38_data, c(chromosome, position), "hg38_data")
    hg38_gr <- GenomicRanges::GRanges(
      seqnames = as.character(hg38_data[[chromosome]]),
      ranges = IRanges::IRanges(
        start = as.numeric(hg38_data[[position]]),
        end = as.numeric(hg38_data[[position]])
      ),
      mcols = hg38_data[, base::setdiff(colnames(hg38_data), c(chromosome, position))]
    )
    GenomeInfoDb::seqlevelsStyle(hg38_gr) <- "NCBI"
  }

  snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
  rsid <- BSgenome::snpsByOverlaps(snps, hg38_gr)
  rsid <- as(rsid, "GRanges")

  S4Vectors::mcols(hg38_gr)$rsid <- NA_character_
  hits <- GenomicRanges::findOverlaps(hg38_gr, rsid, ignore.strand = TRUE)
  S4Vectors::mcols(hg38_gr)$rsid[S4Vectors::queryHits(hits)] <-
    S4Vectors::mcols(rsid)$RefSNP_id[S4Vectors::subjectHits(hits)]

  hg38_gr |>
    sort() |>
    tibble::as_tibble() |>
    dplyr::rename_with(~ gsub("^mcols\\.", "", .x)) |>
    dplyr::rename(chr = seqnames, position = start)
}

format_eqtl_df <- function(eqtl_file_path, ensembl_id, sdY = 1) {
  if (!file.exists(eqtl_file_path)) {
    stop("GTEx eQTL parquet file does not exist: ", eqtl_file_path, call. = FALSE)
  }

  eqtl_df <- arrow::read_parquet(eqtl_file_path) |>
    dplyr::mutate(
      chr = stringr::str_split_fixed(variant_id, "_", 5)[, 1] |>
        stringr::str_remove("^chr") |>
        as.character(),
      position = stringr::str_split_fixed(variant_id, "_", 5)[, 2] |>
        as.integer(),
      other_allele = stringr::str_split_fixed(variant_id, "_", 5)[, 3] |>
        as.character(),
      effect_allele = stringr::str_split_fixed(variant_id, "_", 5)[, 4] |>
        as.character(),
      gene_id = stringr::str_remove(gene_id, "\\.\\d+$"),
      eaf = af,
      maf = ifelse(af > 0.5, 1 - af, af),
      beta = slope,
      varbeta = slope_se^2,
      se = slope_se,
      pvalue = pval_nominal,
      sdY = sdY
    )

  eqtl_gene <- eqtl_df |>
    dplyr::filter(gene_id == ensembl_id)

  if (nrow(eqtl_gene) == 0) {
    message(sprintf("[SKIP] No eQTL found for gene %s", ensembl_id))
    return(NULL)
  }

  eqtl_gene |>
    annotate_snp_with_rsid_hg38() |>
    dplyr::select(
      gene_id, rsid, chr, position, effect_allele, other_allele,
      eaf, maf, beta, varbeta, pvalue, sdY, variant_id
    ) |>
    dplyr::filter(!is.na(rsid), !is.na(pvalue), !is.na(beta), !is.na(varbeta), !is.na(eaf))
}

standardize_pqtl_sumstats <- function(pqtl_df) {
  if (!("chr" %in% colnames(pqtl_df)) || !("position" %in% colnames(pqtl_df))) {
    if (all(c("chromosome_hg38", "base_pair_location_hg38") %in% colnames(pqtl_df))) {
      pqtl_df <- annotate_snp_with_rsid_hg38(
        pqtl_df,
        chromosome = "chromosome_hg38",
        position = "base_pair_location_hg38"
      )
    }
  }

  if (!("pvalue" %in% colnames(pqtl_df)) && "p_value" %in% colnames(pqtl_df)) {
    pqtl_df <- dplyr::rename(pqtl_df, pvalue = p_value)
  }
  if (!("maf" %in% colnames(pqtl_df)) && "MAF" %in% colnames(pqtl_df)) {
    pqtl_df <- dplyr::rename(pqtl_df, maf = MAF)
  }
  if (!("maf" %in% colnames(pqtl_df)) && "eaf" %in% colnames(pqtl_df)) {
    pqtl_df <- dplyr::mutate(pqtl_df, maf = ifelse(eaf > 0.5, 1 - eaf, eaf))
  }
  if (!("varbeta" %in% colnames(pqtl_df))) {
    se_col <- intersect(c("se", "SE", "StdErr", "standard_error"), colnames(pqtl_df))
    if (length(se_col) > 0) {
      pqtl_df <- dplyr::mutate(pqtl_df, varbeta = .data[[se_col[[1]]]]^2)
    }
  }

  pqtl_df
}

prepare_pqtl_region <- function(pqtl_path, region_chr, region_start, region_end) {
  pqtl_df <- data.table::fread(pqtl_path, data.table = FALSE)
  pqtl_df <- standardize_pqtl_sumstats(pqtl_df)

  required_cols <- c(
    "chr", "position", "rsid", "effect_allele", "other_allele",
    "beta", "varbeta", "pvalue", "maf"
  )
  assert_required_columns(pqtl_df, required_cols, "pQTL summary statistics")

  pqtl_df |>
    dplyr::mutate(
      chr = as.character(chr),
      position = as.integer(position),
      rsid = dplyr::na_if(as.character(rsid), "")
    ) |>
    dplyr::filter(
      chr == as.character(region_chr),
      position >= as.integer(region_start),
      position <= as.integer(region_end),
      !is.na(rsid),
      !is.na(beta),
      !is.na(varbeta),
      !is.na(pvalue),
      !is.na(maf)
    )
}

match_coloc_region <- function(dataset1, dataset2, suffix1 = ".eqtl", suffix2 = ".pqtl") {
  ea1 <- paste0("effect_allele", suffix1)
  oa1 <- paste0("other_allele", suffix1)
  beta2 <- paste0("beta", suffix2)
  eaf2 <- paste0("eaf", suffix2)
  ea2 <- paste0("effect_allele", suffix2)
  oa2 <- paste0("other_allele", suffix2)

  matched <- dataset1 |>
    dplyr::filter(rsid %in% dataset2$rsid) |>
    dplyr::left_join(dataset2, by = "rsid", suffix = c(suffix1, suffix2)) |>
    dplyr::filter(
      !is.na(.data[[ea1]]),
      !is.na(.data[[oa1]]),
      !is.na(.data[[ea2]]),
      !is.na(.data[[oa2]])
    ) |>
    dplyr::mutate(
      same_alleles = .data[[ea1]] == .data[[ea2]] & .data[[oa1]] == .data[[oa2]],
      flipped_alleles = .data[[ea1]] == .data[[oa2]] & .data[[oa1]] == .data[[ea2]],
      alleles_match = same_alleles | flipped_alleles
    ) |>
    dplyr::filter(alleles_match)

  if (nrow(matched) == 0) {
    return(matched)
  }

  matched |>
    dplyr::mutate(
      !!beta2 := ifelse(flipped_alleles, -.data[[beta2]], .data[[beta2]]),
      !!eaf2 := if (eaf2 %in% colnames(matched)) {
        ifelse(flipped_alleles, 1 - .data[[eaf2]], .data[[eaf2]])
      } else {
        NA_real_
      }
    )
}

fix_duplicate_rsid <- function(df, rsid_col = "rsid") {
  df |>
    dplyr::group_by(.data[[rsid_col]]) |>
    dplyr::mutate(dup_index = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      !!rsid_col := ifelse(
        duplicated(.data[[rsid_col]]) | duplicated(.data[[rsid_col]], fromLast = TRUE),
        paste0(.data[[rsid_col]], ".", dup_index),
        .data[[rsid_col]]
      )
    ) |>
    dplyr::select(-dup_index)
}

write_coloc_result <- function(
  coloc_results,
  output_dir,
  gene,
  region,
  eQTL_tissue,
  gene_symbol = gene,
  ensembl = NA_character_,
  input_format = NA_character_
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  coloc_results_summary <- as.data.frame(t(coloc_results$summary)) |>
    dplyr::mutate(
      gene_name = gene,
      gene_symbol = gene_symbol,
      ensembl = ensembl,
      region = region,
      tissue = eQTL_tissue,
      input_format = input_format,
      pp4_cond = PP.H4.abf / (PP.H3.abf + PP.H4.abf)
    ) |>
    dplyr::select(gene_name, gene_symbol, ensembl, region, tissue, input_format, dplyr::everything())

  coloc_results_full <- coloc_results$results
  output_prefix <- paste(gene, region, eQTL_tissue, sep = "_")

  data.table::fwrite(
    coloc_results_summary,
    file = file.path(output_dir, paste0(output_prefix, "_coloc_results_summary.tsv")),
    sep = "\t"
  )
  data.table::fwrite(
    coloc_results_full,
    file = file.path(output_dir, paste0(output_prefix, "_coloc_results_full.tsv")),
    sep = "\t"
  )

  invisible(coloc_results_summary)
}

run_coloc <- function(
  pqtl_info,
  eQTL_tissue,
  coloc_result_dir,
  eqtl_dir = "/datasets/cbmr_shared/resources/gtex/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/",
  sdY = 1
) {
  required_region_cols <- c("range_bp", "CHR", "start", "end", "Gene", "ensembl", "file_path", "N")
  assert_required_columns(pqtl_info, required_region_cols, "region input")

  region_id <- pqtl_info$range_bp[[1]]
  region_chr <- pqtl_info$CHR[[1]]
  region_start <- pqtl_info$start[[1]]
  region_end <- pqtl_info$end[[1]]
  gene <- if ("result_name" %in% colnames(pqtl_info)) pqtl_info$result_name[[1]] else pqtl_info$Gene[[1]]
  gene_symbol <- if ("gene_symbol" %in% colnames(pqtl_info)) pqtl_info$gene_symbol[[1]] else pqtl_info$Gene[[1]]
  output_region <- if ("result_region" %in% colnames(pqtl_info)) pqtl_info$result_region[[1]] else region_id
  input_format <- if ("input_format" %in% colnames(pqtl_info)) pqtl_info$input_format[[1]] else NA_character_
  ensembl_id <- pqtl_info$ensembl[[1]]
  pqtl_path <- pqtl_info$file_path[[1]]
  pqtl_sample_size <- unique(pqtl_info$N)[[1]]

  cat(sprintf(
    "[region] Process %s | gene %s (%s), region %s\n",
    gene, gene_symbol, ensembl_id, region_id
  ))

  pqtl_target_region_df <- prepare_pqtl_region(
    pqtl_path = pqtl_path,
    region_chr = region_chr,
    region_start = region_start,
    region_end = region_end
  )

  if (nrow(pqtl_target_region_df) == 0) {
    cat(sprintf("[SKIP] No pQTL variants found for %s in region %s\n", gene, region_id))
    return(NULL)
  }

  eqtl_path <- file.path(
    eqtl_dir,
    paste0(eQTL_tissue, ".v10.allpairs.chr", region_chr, ".parquet")
  )
  eqtl_sumstat <- format_eqtl_df(eqtl_path, ensembl_id, sdY = sdY)

  if (is.null(eqtl_sumstat)) {
    cat(sprintf("[SKIP] No eQTL found for gene %s (%s) in tissue %s\n", gene_symbol, ensembl_id, eQTL_tissue))
    return(NULL)
  }

  eqtl_target_region_df <- eqtl_sumstat |>
    dplyr::filter(
      chr == as.character(region_chr),
      position >= as.integer(region_start),
      position <= as.integer(region_end),
      !is.na(rsid)
    )

  if (nrow(eqtl_target_region_df) == 0) {
    cat(sprintf(
      "[SKIP] No eQTL variants found for gene %s (%s) in tissue %s and region %s\n",
      gene_symbol, ensembl_id, eQTL_tissue, region_id
    ))
    return(NULL)
  }

  coloc_region <- match_coloc_region(
    eqtl_target_region_df,
    pqtl_target_region_df,
    suffix1 = ".eqtl",
    suffix2 = ".pqtl"
  ) |>
    fix_duplicate_rsid()

  if (nrow(coloc_region) == 0) {
    cat(sprintf(
      "[SKIP] No allele-aligned shared variants for gene %s (%s), tissue %s, region %s\n",
      gene_symbol, ensembl_id, eQTL_tissue, region_id
    ))
    return(NULL)
  }

  ds_pqtl <- list(
    snp = coloc_region$rsid,
    beta = coloc_region$beta.pqtl,
    varbeta = coloc_region$varbeta.pqtl,
    pvalues = coloc_region$pvalue.pqtl,
    MAF = coloc_region$maf.pqtl,
    N = pqtl_sample_size,
    type = "quant"
  )
  ds_eqtl <- list(
    snp = coloc_region$rsid,
    beta = coloc_region$beta.eqtl,
    varbeta = coloc_region$varbeta.eqtl,
    pvalues = coloc_region$pvalue.eqtl,
    MAF = coloc_region$maf.eqtl,
    sdY = sdY,
    type = "quant"
  )

  cat(sprintf("[Run coloc] %s | %s | %s | nsnps=%s\n", gene, output_region, eQTL_tissue, nrow(coloc_region)))
  coloc_result <- coloc::coloc.abf(dataset1 = ds_pqtl, dataset2 = ds_eqtl)
  write_coloc_result(
    coloc_result,
    coloc_result_dir,
    gene,
    output_region,
    eQTL_tissue,
    gene_symbol = gene_symbol,
    ensembl = ensembl_id,
    input_format = input_format
  )

  coloc_result
}

combine_coloc_summaries <- function(coloc_result_dir) {
  coloc_files <- list.files(
    coloc_result_dir,
    pattern = "_coloc_results_summary.tsv$",
    full.names = TRUE
  )

  if (length(coloc_files) == 0) {
    return(data.frame())
  }

  coloc_files |>
    lapply(data.table::fread) |>
    dplyr::bind_rows()
}
