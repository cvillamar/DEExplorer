# Wrapper API intentionally mirrors expected DEExplorer user entrypoint.

normalize_marralm_input <- function(in_marralm, arg_expr = substitute(in_marralm)) {
  if (is.null(in_marralm)) {
    stop("`in_marralm` is required and must contain at least one MArrayLM fit.", call. = FALSE)
  }

  if (inherits(in_marralm, "MArrayLM")) {
    fit_name <- paste(deparse(arg_expr, width.cutoff = 500L), collapse = "")
    return(list(fits = stats::setNames(list(in_marralm), fit_name)))
  }

  if (!is.list(in_marralm)) {
    stop("`in_marralm` must be a MArrayLM object or a list of MArrayLM objects.", call. = FALSE)
  }

  is_fit <- vapply(in_marralm, inherits, logical(1), what = "MArrayLM")
  if (!all(is_fit)) {
    flattened_hint <- !is.null(names(in_marralm)) &&
      anyDuplicated(names(in_marralm)) &&
      "coefficients" %in% names(in_marralm)
    if (flattened_hint) {
      starts <- which(names(in_marralm) == "coefficients")
      ends <- c(starts[-1L] - 1L, length(in_marralm))
      chunks <- Map(function(start, end) in_marralm[start:end], starts, ends)

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
            "`in_marralm` looked like flattened `c(efit, efit2, ...)` input. ",
            "Recovered fits automatically, but `list(efit = efit, efit2 = efit2)` is preferred."
          ),
          call. = FALSE
        )
        return(list(fits = recovered))
      }
    }
    stop("Every element of `in_marralm` must inherit from limma::MArrayLM.", call. = FALSE)
  }

  fit_names <- names(in_marralm)
  if (is.null(fit_names)) {
    fit_names <- rep("", length(in_marralm))
  }
  missing_names <- is_missing_string(fit_names)
  fit_names[missing_names] <- paste0("fit", seq_along(in_marralm))[missing_names]
  fit_names <- make.unique(fit_names)
  names(in_marralm) <- fit_names

  list(fits = in_marralm)
}

normalize_gsea_input <- function(in_gsea) {
  out <- list(
    fgsea_dir = NULL,
    msigdb_genesets = NULL,
    msigdb_path = NULL
  )

  if (is.null(in_gsea)) {
    return(out)
  }

  if (is.data.frame(in_gsea)) {
    out$msigdb_genesets <- in_gsea
    return(out)
  }

  if (is.character(in_gsea) && length(in_gsea) == 1L) {
    gsea_path <- normalizePath(in_gsea, winslash = "/", mustWork = FALSE)
    if (dir.exists(gsea_path)) {
      out$fgsea_dir <- gsea_path
      return(out)
    }
    if (file.exists(gsea_path) && grepl("\\.rds$", gsea_path, ignore.case = TRUE)) {
      out$msigdb_path <- gsea_path
      return(out)
    }
    stop("`in_gsea` path must point to a directory with fGSEA TSV files or an `.rds` geneset file.", call. = FALSE)
  }

  if (is.list(in_gsea)) {
    pick <- function(x, key) {
      if (is.null(names(x)) || !key %in% names(x)) {
        return(NULL)
      }
      x[[key]]
    }

    has_config_names <- !is.null(names(in_gsea)) &&
      any(names(in_gsea) %in% c("fgsea_dir", "fgsea", "msigdb_path", "msigdb_genesets", "genesets", "msigdb"))
    if (has_config_names) {
      fgsea_dir <- pick(in_gsea, "fgsea_dir") %||% pick(in_gsea, "fgsea")
      msigdb_genesets <- pick(in_gsea, "msigdb_genesets") %||% pick(in_gsea, "genesets") %||% pick(in_gsea, "msigdb")
      msigdb_path <- pick(in_gsea, "msigdb_path")

      if (!is.null(fgsea_dir)) {
        fgsea_dir <- normalizePath(as.character(fgsea_dir)[[1L]], winslash = "/", mustWork = FALSE)
        if (!dir.exists(fgsea_dir)) {
          stop("`in_gsea$fgsea_dir` does not exist.", call. = FALSE)
        }
      }
      if (!is.null(msigdb_path)) {
        msigdb_path <- normalizePath(as.character(msigdb_path)[[1L]], winslash = "/", mustWork = FALSE)
        if (!file.exists(msigdb_path)) {
          stop("`in_gsea$msigdb_path` does not exist.", call. = FALSE)
        }
      }

      out$fgsea_dir <- fgsea_dir
      out$msigdb_genesets <- msigdb_genesets
      out$msigdb_path <- msigdb_path
      return(out)
    }

    if (!is.null(names(in_gsea)) && all(nzchar(names(in_gsea)))) {
      out$msigdb_genesets <- in_gsea
      return(out)
    }
  }

  stop(
    paste0(
      "`in_gsea` must be one of: NULL, a directory path, an .rds path, ",
      "a geneset data.frame, a named geneset list, or a named config list ",
      "with `fgsea_dir`/`msigdb_path`/`msigdb_genesets`."
    ),
    call. = FALSE
  )
}

#' Build a DEExplorer Shiny app bundle and write ui.R/server.R
#'
#' @param in_dgelist `edgeR::DGEList` object.
#' @param in_marralm One `MArrayLM` or a named list of one or more `MArrayLM` fits.
#' @param in_gsea Optional fGSEA / geneset input.
#' @param out Output folder where `ui.R`, `server.R` and `app-data.rds` are written.
#' @param title App title.
#' @param overwrite Overwrite existing generated files.
#' @param launch Run the generated app immediately.
#' @param sample_id_col Sample id column in `dge$samples`.
#' @param gene_id_col Gene id column in `dge$genes`.
#' @param gene_symbol_col Gene symbol column in `dge$genes`.
#' @param fdr_cutoff Significance threshold for DE coloring.
#' @param prior_count Prior count for log-CPM.
#' @param enrichr_databases Enrichr database names offered in the UI.
#'
#' @return Invisible list with generated file paths.
deexplorer <- function(
  in_dgelist,
  in_marralm,
  in_gsea = NULL,
  out = ".",
  title = "DEExplorer",
  overwrite = TRUE,
  launch = FALSE,
  sample_id_col = "Sample",
  gene_id_col = "ENSEMBL",
  gene_symbol_col = "SYMBOL",
  fdr_cutoff = 0.05,
  prior_count = 2,
  enrichr_databases = default_enrichr_databases()
) {
  fits <- normalize_marralm_input(in_marralm = in_marralm, arg_expr = substitute(in_marralm))
  gsea <- normalize_gsea_input(in_gsea)

  app_dir <- normalizePath(out, winslash = "/", mustWork = FALSE)

  call_args <- c(
    list(
      dge = in_dgelist,
      app_dir = app_dir,
      title = title,
      overwrite = overwrite,
      launch = launch,
      msigdb_genesets = gsea$msigdb_genesets,
      msigdb_path = gsea$msigdb_path,
      fgsea_dir = gsea$fgsea_dir,
      sample_id_col = sample_id_col,
      gene_id_col = gene_id_col,
      gene_symbol_col = gene_symbol_col,
      fdr_cutoff = fdr_cutoff,
      prior_count = prior_count,
      enrichr_databases = enrichr_databases
    ),
    fits$fits
  )

  do.call(write_deexplorer_app, call_args)
}
