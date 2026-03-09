# Wrapper API intentionally mirrors expected DEExplorer user entrypoint.

normalize_marralm_input <- function(input_MArrayLM, arg_expr = substitute(input_MArrayLM)) {
  if (is.null(input_MArrayLM)) {
    stop("`input_MArrayLM` is required and must contain at least one MArrayLM fit.", call. = FALSE)
  }

  if (inherits(input_MArrayLM, "MArrayLM")) {
    fit_name <- paste(deparse(arg_expr, width.cutoff = 500L), collapse = "")
    return(list(fits = stats::setNames(list(input_MArrayLM), fit_name)))
  }

  if (!is.list(input_MArrayLM)) {
    stop("`input_MArrayLM` must be a MArrayLM object or a list of MArrayLM objects.", call. = FALSE)
  }

  is_fit <- vapply(input_MArrayLM, inherits, logical(1), what = "MArrayLM")
  if (!all(is_fit)) {
    flattened_hint <- !is.null(names(input_MArrayLM)) &&
      anyDuplicated(names(input_MArrayLM)) &&
      "coefficients" %in% names(input_MArrayLM)
    if (flattened_hint) {
      starts <- which(names(input_MArrayLM) == "coefficients")
      ends <- c(starts[-1L] - 1L, length(input_MArrayLM))
      chunks <- Map(function(start, end) input_MArrayLM[start:end], starts, ends)

      can_recover <- length(chunks) >= 1L &&
        all(vapply(chunks, function(chunk) {
          all(c("coefficients", "p.value", "t", "Amean") %in% names(chunk))
        }, logical(1)))

      if (can_recover) {
        recovered <- lapply(chunks, function(chunk) {
          class(chunk) <- "MArrayLM"
          chunk
        })
        names(recovered) <- paste0("fit", seq_along(recovered))
        warning(
          paste0(
            "`input_MArrayLM` looked like flattened `c(efit, efit2, ...)` input. ",
            "Recovered fits automatically, but `list(efit = efit, efit2 = efit2)` is preferred."
          ),
          call. = FALSE
        )
        return(list(fits = recovered))
      }
    }
    stop("Every element of `input_MArrayLM` must inherit from limma::MArrayLM.", call. = FALSE)
  }

  fit_names <- names(input_MArrayLM)
  if (is.null(fit_names)) {
    fit_names <- rep("", length(input_MArrayLM))
  }
  missing_names <- is_missing_string(fit_names)
  fit_names[missing_names] <- paste0("fit", seq_along(input_MArrayLM))[missing_names]
  fit_names <- make.unique(fit_names)
  names(input_MArrayLM) <- fit_names

  list(fits = input_MArrayLM)
}

normalize_gsea_files_input <- function(input_gsea_files) {
  out <- list(
    fgsea_files = NULL,
    msigdb_genesets = NULL
  )

  if (is.null(input_gsea_files)) {
    return(out)
  }

  if (is.character(input_gsea_files)) {
    out$fgsea_files <- input_gsea_files
    return(out)
  }

  if (is.list(input_gsea_files)) {
    if (!is.null(names(input_gsea_files))) {
      pick <- function(x, key) {
        if (!key %in% names(x)) return(NULL)
        x[[key]]
      }
      out$fgsea_files <- pick(input_gsea_files, "fgsea_files") %||%
        pick(input_gsea_files, "fgsea")
      out$msigdb_genesets <- pick(input_gsea_files, "msigdb_genesets") %||%
        pick(input_gsea_files, "genesets") %||%
        pick(input_gsea_files, "msigdb")
      if (!is.null(out$fgsea_files) || !is.null(out$msigdb_genesets)) {
        return(out)
      }
      if (all(nzchar(names(input_gsea_files)))) {
        out$msigdb_genesets <- input_gsea_files
        return(out)
      }
    }
  }

  stop(
    paste0(
      "`input_gsea_files` must be one of: NULL, a character vector of fGSEA TSV file paths, ",
      "or a named list with `fgsea_files` and/or `msigdb_genesets`."
    ),
    call. = FALSE
  )
}

#' Build a DEExplorer Shiny app bundle and write ui.R/server.R
#'
#' @param input_DGEList `edgeR::DGEList` object.
#' @param input_MArrayLM One `MArrayLM` or a named list of one or more `MArrayLM` fits.
#' @param input_gsea_files Optional character vector of fGSEA TSV file paths, or a named
#'   list with `fgsea_files` (character vector) and/or `msigdb_genesets` (named list).
#' @param msigdb_path Optional path to an `.rds` file containing a named list of gene sets
#'   (gene-set names as keys, character vectors of gene symbols as values).
#' @param out Output folder where `ui.R`, `server.R` and `app-data.rds` are written.
#' @param title App title.
#' @param overwrite Overwrite existing generated files.
#' @param launch Run the generated app immediately.
#' @param sample_id_col Sample id column in `dge$samples`.
#' @param gene_id_col Gene id column in `dge$genes`.
#' @param gene_symbol_col Gene symbol column in `dge$genes`.
#' @param fdr_cutoff Significance threshold for DE coloring.
#' @param prior_count Prior count for log-CPM.
#'
#' @return Invisible list with generated file paths.
deexplorer <- function(
  input_DGEList,
  input_MArrayLM,
  input_gsea_files = NULL,
  msigdb_path = NULL,
  out = ".",
  title = "DEExplorer",
  overwrite = TRUE,
  launch = FALSE,
  sample_id_col = "Sample",
  gene_id_col = "ENSEMBL",
  gene_symbol_col = "SYMBOL",
  fdr_cutoff = 0.05,
  prior_count = 2
) {
  fits <- normalize_marralm_input(input_MArrayLM = input_MArrayLM, arg_expr = substitute(input_MArrayLM))
  gsea <- normalize_gsea_files_input(input_gsea_files)

  app_dir <- normalizePath(out, winslash = "/", mustWork = FALSE)

  call_args <- c(
    list(
      dge = input_DGEList,
      app_dir = app_dir,
      title = title,
      overwrite = overwrite,
      launch = launch,
      msigdb_genesets = gsea$msigdb_genesets,
      msigdb_path = msigdb_path,
      fgsea_files = gsea$fgsea_files,
      sample_id_col = sample_id_col,
      gene_id_col = gene_id_col,
      gene_symbol_col = gene_symbol_col,
      fdr_cutoff = fdr_cutoff,
      prior_count = prior_count
    ),
    fits$fits
  )

  do.call(write_deexplorer_app, call_args)
}
