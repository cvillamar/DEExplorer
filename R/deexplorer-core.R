`%||%` <-
function (x, y) 
{
    if (is.null(x) || length(x) == 0L) {
        y
    }
    else {
        x
    }
}
build_de_payload <-
function (fits, fit_names, gene_df, gene_lookup, fdr_cutoff) 
{
    tables <- list()
    catalog_rows <- list()
    for (fit_index in seq_along(fits)) {
        fit <- fits[[fit_index]]
        fit_name <- fit_names[[fit_index]]
        method <- infer_fit_method(fit)
        method_label <- build_method_label(method, fit)
        for (contrast_name in colnames(fit$coefficients)) {
            tt <- build_single_de_table(fit = fit, fit_name = fit_name, 
                contrast_name = contrast_name, gene_df = gene_df, 
                gene_lookup = gene_lookup, fdr_cutoff = fdr_cutoff)
            contrast_key <- paste(fit_name, contrast_name, sep = "::")
            tables[[contrast_key]] <- tt
            catalog_rows[[length(catalog_rows) + 1L]] <- data.frame(contrast_key = contrast_key, 
                fit_label = fit_name, contrast = contrast_name, 
                method = method, method_label = method_label, 
                display_label = paste0(fit_name, " / ", contrast_name, 
                  " [", method_label, "]"), stringsAsFactors = FALSE)
        }
    }
    contrast_catalog <- do.call(rbind, catalog_rows)
    rownames(contrast_catalog) <- contrast_catalog$contrast_key
    list(tables = tables, contrast_catalog = contrast_catalog)
}
build_de_plot <-
function (de_df, x_col, y_col, title, source_id, highlighted_gene_keys = character(),
    manual_gene_keys = character(), leading_edge_keys = character(),
    active_gene_key = character(), x_title, y_title)
{
    layers <- subset_plot_layers(de_df = de_df, highlighted_gene_keys = highlighted_gene_keys,
        manual_gene_keys = manual_gene_keys, leading_edge_keys = leading_edge_keys,
        active_gene_key = active_gene_key)
    has_highlight <- length(highlighted_gene_keys) > 0L || length(manual_gene_keys) > 0L
    bg_nonsig_color <- if (has_highlight) "rgba(113, 128, 150, 0.15)" else "rgba(113, 128, 150, 0.50)"
    bg_sig_color    <- if (has_highlight) "rgba(197, 48, 48, 0.20)"  else "rgba(197, 48, 48, 0.78)"
    geneset_data <- rbind(layers$geneset, layers$manual)
    xf <- stats::as.formula(paste0("~`", x_col, "`"))
    yf <- stats::as.formula(paste0("~`", y_col, "`"))
    p <- plotly::plot_ly(source = source_id)
    add_fixed_trace <- function(p, data, x_col, y_col, name, marker,
                                legendgroup = name) {
        plotly::add_markers(p, x = data[[x_col]], y = data[[y_col]],
            text = data$hover_text, hoverinfo = "text",
            customdata = data$gene_key, name = name,
            legendgroup = legendgroup, marker = marker, inherit = FALSE)
    }
    p <- add_fixed_trace(p, layers$nonsig, x_col, y_col,
        "FDR >= 0.05", list(color = bg_nonsig_color, size = 6))
    p <- add_fixed_trace(p, layers$sig, x_col, y_col,
        "FDR < 0.05", list(color = bg_sig_color, size = 7))
    p <- add_fixed_trace(p, geneset_data, x_col, y_col,
        "Gene set", list(color = "#e76f51", size = 9, opacity = 0.85))
    p <- add_fixed_trace(p, layers$leading, x_col, y_col,
        "Leading edge", list(color = "rgba(255, 255, 255, 0)",
        size = 12, symbol = "diamond-open",
        line = list(color = "#8d6a00", width = 3)))
    p <- add_fixed_trace(p, layers$active, x_col, y_col,
        "Active gene", list(color = "#0b6e4f", size = 13, symbol = "star"))
    p <- plotly::layout(p, title = list(text = title), xaxis = list(title = x_title),
        yaxis = list(title = y_title), legend = list(orientation = "h",
            y = -0.18))
    p <- plotly::event_register(p, "plotly_click")
    p <- plotly::event_register(p, "plotly_hover")
    p
}
build_expression_matrices <-
function (dge, gene_df, prior_count) 
{
    counts <- dge$counts
    rownames(counts) <- gene_df$gene_key
    cpm <- edgeR::cpm(dge, log = FALSE)
    rownames(cpm) <- gene_df$gene_key
    lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = prior_count)
    rownames(lcpm) <- gene_df$gene_key
    list(counts = counts, cpm = cpm, lcpm = lcpm)
}
build_fgsea_barplot <-
function (fgsea_df, n = 10L)
{
    shiny::validate(shiny::need(
        !is.null(fgsea_df) && nrow(fgsea_df) > 0L,
        "No fGSEA results available for this contrast."
    ))
    fgsea_df <- fgsea_df[order(-abs(fgsea_df$NES)), , drop = FALSE]
    fgsea_df <- utils::head(fgsea_df, n)
    fgsea_df$short_pathway <- ifelse(
        nchar(fgsea_df$pathway) > 50L,
        paste0(substr(fgsea_df$pathway, 1L, 47L), "..."),
        fgsea_df$pathway
    )
    fgsea_df$neg_log10_padj <- -log10(pmax(fgsea_df$padj, .Machine$double.xmin))
    fgsea_df$short_pathway <- factor(fgsea_df$short_pathway,
        levels = rev(fgsea_df$short_pathway))
    plotly::layout(
        plotly::plot_ly(
            data = fgsea_df, x = ~NES, y = ~short_pathway,
            type = "bar", orientation = "h",
            marker = list(
                color = ~neg_log10_padj,
                colorscale = list(c(0, "#b0bec5"), c(0.5, "#e57373"), c(1, "#b71c1c")),
                colorbar = list(title = "-log10(padj)", len = 0.6),
                line = list(width = 0)
            ),
            text = ~paste0("<b>", pathway, "</b><br>",
                           "NES: ", format_number(NES), "<br>",
                           "padj: ", format_pvalue(padj)),
            hoverinfo = "text"
        ),
        title = list(text = "Top fGSEA Enrichment (by |NES|)", font = list(size = 13)),
        xaxis = list(title = "NES", zeroline = TRUE,
            zerolinecolor = "#999999", zerolinewidth = 1),
        yaxis = list(title = ""),
        margin = list(l = 200)
    )
}
build_fgsea_summary <-
function (fgsea_row) 
{
    if (!nrow(fgsea_row)) {
        return(shiny::tagList(shiny::h4("Precomputed fGSEA"), 
            shiny::div(class = "deexplorer-note", "No precomputed fGSEA result was found for the selected gene set and contrast.")))
    }
    leading_edge <- fgsea_row$leading_edge_genes[[1L]]
    shiny::tagList(shiny::h4("Precomputed fGSEA"), shiny::div(class = "deexplorer-metric", 
        shiny::strong("Collection: "), fgsea_row$collection[[1L]]), 
        shiny::div(class = "deexplorer-metric", shiny::strong("NES: "), 
            format_number(fgsea_row$NES[[1L]])), shiny::div(class = "deexplorer-metric", 
            shiny::strong("P.Value: "), format_pvalue(fgsea_row$pval[[1L]])), 
        shiny::div(class = "deexplorer-metric", shiny::strong("adj.P.Val: "), 
            format_pvalue(fgsea_row$padj[[1L]])), shiny::div(class = "deexplorer-metric", 
            shiny::strong("Leading edge: "), paste(leading_edge, 
                collapse = ", ")))
}
build_gene_annotation <-
function (dge, gene_id_col, gene_symbol_col) 
{
    feature_name <- rownames(dge$counts)
    if (is.null(feature_name)) {
        feature_name <- paste0("feature_", seq_len(nrow(dge$counts)))
    }
    raw_genes <- dge$genes
    if (is.null(raw_genes)) {
        raw_genes <- data.frame(row.names = seq_len(nrow(dge$counts)))
    }
    raw_genes <- as.data.frame(raw_genes, stringsAsFactors = FALSE, 
        check.names = FALSE)
    gene_id <- if (gene_id_col %in% colnames(raw_genes)) {
        as.character(raw_genes[[gene_id_col]])
    }
    else {
        feature_name
    }
    symbol <- if (gene_symbol_col %in% colnames(raw_genes)) {
        as.character(raw_genes[[gene_symbol_col]])
    }
    else {
        feature_name
    }
    gene_id <- coalesce_character(gene_id, feature_name)
    symbol <- coalesce_character(symbol, feature_name)
    gene_key <- make.unique(ifelse(is_missing_string(gene_id), 
        feature_name, gene_id))
    gene_df <- data.frame(gene_key = gene_key, gene_id = gene_id, 
        symbol = symbol, feature_name = feature_name, stringsAsFactors = FALSE, 
        check.names = FALSE)
    passthrough_cols <- setdiff(colnames(raw_genes), c(gene_id_col, 
        gene_symbol_col))
    for (col_name in passthrough_cols) {
        gene_df[[col_name]] <- raw_genes[[col_name]]
    }
    gene_df$symbol_upper <- toupper(gene_df$symbol)
    gene_df$gene_id_upper <- toupper(gene_df$gene_id)
    gene_df$feature_upper <- toupper(gene_df$feature_name)
    rownames(gene_df) <- gene_df$gene_key
    gene_df
}
build_gene_boxplot <-
function (gene_expr_df, gene_label, group_col) 
{
    shiny::validate(shiny::need(nrow(gene_expr_df) > 0L, "Hover or select a gene in the table or plots to show its log-CPM profile."))
    plotly::layout(plotly::plot_ly(data = gene_expr_df, x = ~group_label, 
        y = ~log_cpm, color = ~group_label, type = "box", boxpoints = "all", 
        jitter = 0.25, pointpos = 0, text = ~hover_text, hoverinfo = "text", 
        showlegend = FALSE), title = list(text = paste0("log-CPM for ", 
        gene_label)), xaxis = list(title = group_col), yaxis = list(title = "log-CPM"))
}
build_gene_lookup <-
function (gene_df) 
{
    split_lookup <- function(values) {
        keep <- !is_missing_string(values)
        split(gene_df$gene_key[keep], toupper(values[keep]))
    }
    first_entry <- function(x) vapply(x, `[`, character(1), 1L)
    symbol_lookup <- split_lookup(gene_df$symbol)
    gene_id_lookup <- split_lookup(gene_df$gene_id)
    feature_lookup <- split_lookup(gene_df$feature_name)
    list(by_key = stats::setNames(gene_df$gene_key, toupper(gene_df$gene_key)), 
        by_symbol = symbol_lookup, by_gene_id = gene_id_lookup, 
        by_feature = feature_lookup, first_by_symbol = first_entry(symbol_lookup), 
        first_by_gene_id = first_entry(gene_id_lookup), first_by_feature = first_entry(feature_lookup))
}
build_iheatmap <-
function (lcpm_matrix, gene_df, sample_df, col_annotations = NULL,
    row_geneset_flags = NULL, cluster_rows = "none", cluster_cols = "none",
    show_row_labels = TRUE, show_col_labels = TRUE)
{
    mat <- t(scale(t(lcpm_matrix)))
    mat[!is.finite(mat)] <- 0
    row_syms <- gene_df[rownames(mat), "symbol"]
    row_syms <- ifelse(is_missing_string(row_syms), rownames(mat), row_syms)
    rownames(mat) <- make.unique(row_syms)
    colnames(mat) <- sample_df[colnames(mat), "sample_id"]
    hm <- iheatmapr::iheatmap(mat, name = "Z-score",
        cluster_rows = cluster_rows, cluster_cols = cluster_cols,
        colors = grDevices::colorRampPalette(
            c("#2166ac", "#f7f7f7", "#b2182b"))(100),
        layout = list(margin = list(l = 20, b = 20)))
    if (show_row_labels) {
        hm <- iheatmapr::add_row_labels(hm, side = "right")
    }
    if (show_col_labels) {
        hm <- iheatmapr::add_col_labels(hm, side = "bottom")
    }
    if (!is.null(col_annotations) && length(col_annotations)) {
        col_ann_df <- as.data.frame(col_annotations,
            stringsAsFactors = FALSE, check.names = FALSE)
        hm <- iheatmapr::add_col_annotation(hm,
            annotation = col_ann_df, side = "bottom")
    }
    if (!is.null(row_geneset_flags) && length(row_geneset_flags)) {
        gs_names <- names(row_geneset_flags)
        gs_colors <- c("In" = "#000000", "Out" = "#d9d9d9")
        for (i in seq_along(gs_names)) {
            nm <- gs_names[i]
            hm <- iheatmapr::add_row_groups(hm,
                groups = row_geneset_flags[[nm]],
                name = if (i == 1L) "Gene-set" else paste0(".gs_", i),
                title = nm,
                colors = gs_colors,
                show_colorbar = (i == 1L),
                side = "right",
                pname = paste0("gs_", i))
        }
    }
    hm
}
build_method_label <-
function (method, fit)
{
    if (identical(method, "treat")) {
        paste0("treat (|logFC| > ", formatC(fit$treat.lfc[[1]], 
            digits = 3, format = "fg"), ")")
    }
    else {
        "eBayes"
    }
}
build_pca_payload <-
function (lcpm, sample_df) 
{
    pca <- stats::prcomp(t(lcpm), center = TRUE, scale. = FALSE)
    score_df <- as.data.frame(pca$x, stringsAsFactors = FALSE, 
        check.names = FALSE)
    score_df$sample_key <- rownames(score_df)
    score_df <- cbind(score_df, sample_df[score_df$sample_key, 
        , drop = FALSE])
    score_df$hover_text <- make_sample_hover_text(sample_df[score_df$sample_key, 
        , drop = FALSE])
    variance_explained <- (pca$sdev^2)/sum(pca$sdev^2)
    variance_df <- data.frame(pc = paste0("PC", seq_along(variance_explained)), 
        variance_explained = variance_explained, stringsAsFactors = FALSE)
    list(scores = score_df, variance = variance_df)
}
build_pca_plot <-
function (bundle, x_pc, y_pc, color_col = NULL, shape_col = NULL, 
    size_col = NULL) 
{
    pca_df <- bundle$pca_scores
    shiny::validate(shiny::need(x_pc %in% colnames(pca_df), "Selected x-axis PC is unavailable."))
    shiny::validate(shiny::need(y_pc %in% colnames(pca_df), "Selected y-axis PC is unavailable."))
    color_values <- if (!is.null(color_col)) 
        pca_df[[color_col]]
    else rep("Default", nrow(pca_df))
    shape_values <- if (!is.null(shape_col)) 
        pca_df[[shape_col]]
    else rep("circle", nrow(pca_df))
    size_values <- if (!is.null(size_col)) 
        pca_df[[size_col]]
    else rep(14, nrow(pca_df))
    plotly::layout(plotly::plot_ly(data = pca_df, x = pca_df[[x_pc]], 
        y = pca_df[[y_pc]], type = "scatter", 
        mode = "markers+text", text = pca_df$sample_id, textposition = "top center", 
        textfont = list(size = 9, color = "#2f3e46"), hovertext = pca_df$hover_text, 
        hoverinfo = "text", source = "pca_plot", marker = list(color = scale_marker_color(color_values), 
        size = scale_marker_size(size_values), symbol = scale_marker_symbol(shape_values), 
        opacity = 0.88, line = list(color = "#ffffff", width = 1), 
        sizemode = "diameter")), 
        title = list(text = paste0("Sample PCA: ", x_pc, " vs ", 
            y_pc)), xaxis = list(title = x_pc), yaxis = list(title = y_pc), 
        dragmode = "select", showlegend = FALSE)
}
build_pca_variance_plot <-
function (bundle) 
{
    variance_df <- bundle$pca_variance
    pc_order <- suppressWarnings(as.integer(sub("^PC", "", variance_df$pc)))
    if (all(is.na(pc_order))) {
        variance_df <- utils::head(variance_df, 8L)
    }
    else {
        variance_df$pc_order <- pc_order
        variance_df <- variance_df[order(variance_df$pc_order), , 
            drop = FALSE]
        variance_df <- utils::head(variance_df, 8L)
    }
    variance_df$pc <- factor(variance_df$pc, levels = variance_df$pc, 
        ordered = TRUE)
    plotly::layout(plotly::plot_ly(data = variance_df, x = ~pc, 
        y = ~variance_explained * 100, type = "bar", marker = list(color = "#7a4f01")), 
        title = list(text = "Variance explained"), xaxis = list(title = "", 
            categoryorder = "array", categoryarray = as.character(variance_df$pc)), 
        yaxis = list(title = "% variance explained"))
}
build_sample_annotation <-
function (dge, sample_id_col) 
{
    sample_df <- as.data.frame(dge$samples, stringsAsFactors = FALSE, 
        check.names = FALSE)
    sample_df$sample_key <- colnames(dge$counts)
    sample_df$sample_id <- if (!is.null(sample_id_col) && sample_id_col %in% 
        colnames(sample_df)) {
        as.character(sample_df[[sample_id_col]])
    }
    else {
        sample_df$sample_key
    }
    rownames(sample_df) <- sample_df$sample_key
    sample_df
}
build_selected_gene_summary <-
function (gene_row) 
{
    shiny::validate(shiny::need(nrow(gene_row) == 1L, "Hover or select a gene to see statistics here."))
    shiny::tagList(shiny::h4(gene_row$display_symbol), shiny::div(class = "deexplorer-metric", 
        shiny::strong("Gene ID: "), gene_row$gene_id), shiny::div(class = "deexplorer-metric", 
        shiny::strong("Contrast: "), gene_row$contrast), shiny::div(class = "deexplorer-metric", 
        shiny::strong("Method: "), gene_row$method_label), shiny::div(class = "deexplorer-metric", 
        shiny::strong("logFC: "), format_number(gene_row$logFC)), 
        shiny::div(class = "deexplorer-metric", shiny::strong("AveExpr: "), 
            format_number(gene_row$AveExpr)), shiny::div(class = "deexplorer-metric", 
            shiny::strong("P.Value: "), format_pvalue(gene_row$P.Value)), 
        shiny::div(class = "deexplorer-metric", shiny::strong("adj.P.Val: "), 
            format_pvalue(gene_row$adj.P.Val)), shiny::div(class = "deexplorer-metric", 
            shiny::strong("t-statistic: "), format_number(gene_row$t)))
}
build_selection_summary <-
function (contrast_label, manual_gene_count, geneset_name, filtered_row_count)
{
    shiny::tagList(shiny::h4("Current selection"), shiny::div(class = "deexplorer-metric",
        shiny::strong("Contrast: "), contrast_label), shiny::div(class = "deexplorer-metric",
        shiny::strong("Manual highlights: "), manual_gene_count),
        shiny::div(class = "deexplorer-metric", shiny::strong("MSigDB gene set: "),
            geneset_name %||% "None"), shiny::div(class = "deexplorer-metric",
            shiny::strong("Rows in table: "), filtered_row_count),
        shiny::div(class = "deexplorer-note", "Significance colors use adj.P.Val < 0.05."))
}
build_single_de_table <-
function (fit, fit_name, contrast_name, gene_df, gene_lookup, 
    fdr_cutoff) 
{
    method <- infer_fit_method(fit)
    top_fun <- if (identical(method, "treat")) 
        limma::topTreat
    else limma::topTable
    tt <- top_fun(fit, coef = contrast_name, number = Inf, sort.by = "t")
    tt <- as.data.frame(tt, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(tt)) {
        return(tt)
    }
    tt$gene_key <- map_table_gene_keys(tt, gene_lookup)
    annotation <- gene_df[tt$gene_key, , drop = FALSE]
    tt$gene_id <- annotation$gene_id
    tt$symbol <- annotation$symbol
    tt$feature_name <- annotation$feature_name
    tt$display_symbol <- ifelse(is_missing_string(tt$symbol), 
        tt$feature_name, tt$symbol)
    tt$fit_label <- fit_name
    tt$contrast <- contrast_name
    tt$contrast_key <- paste(fit_name, contrast_name, sep = "::")
    tt$method <- method
    tt$method_label <- build_method_label(method, fit)
    tt$rank <- seq_len(nrow(tt))
    if (!"t" %in% colnames(tt)) {
        tt$t <- NA
    }
    if (!"B" %in% colnames(tt)) {
        tt$B <- NA
    }
    tt$significant <- tt$adj.P.Val < fdr_cutoff
    tt$neg_log10_p <- -log10(pmax(tt$P.Value, .Machine$double.xmin))
    tt$hover_text <- make_gene_hover_text(tt)
    tt
}
build_tissue_expression_barplot <-
function (expr_df)
{
    shiny::validate(shiny::need(
        !is.null(expr_df) && nrow(expr_df) > 0L,
        "No tissue expression data available for this gene."
    ))
    expr_df <- expr_df[order(-expr_df$rna_value), , drop = FALSE]
    expr_df <- utils::head(expr_df, 20L)
    expr_df$organ <- factor(expr_df$organ, levels = rev(expr_df$organ))
    plotly::layout(
        plotly::plot_ly(
            data = expr_df, x = ~rna_value, y = ~organ,
            type = "bar", orientation = "h",
            marker = list(color = "#3489ca"),
            text = ~paste0("<b>", organ, "</b><br>",
                           "RNA: ", format_number(rna_value), " TPM"),
            hoverinfo = "text"
        ),
        title = list(text = "Tissue Expression (RNA, TPM)", font = list(size = 13)),
        xaxis = list(title = "TPM"),
        yaxis = list(title = ""),
        margin = list(l = 140)
    )
}
coalesce_character <-
function (primary, fallback) 
{
    primary <- as.character(primary)
    fallback <- as.character(fallback)
    primary[is_missing_string(primary)] <- fallback[is_missing_string(primary)]
    primary
}
collect_fit_inputs <-
function (...) 
{
    fits <- list(...)
    fit_exprs <- as.list(substitute(list(...)))[-1L]
    if (!length(fits)) {
        stop("At least one limma MArrayLM object must be supplied through `...`.", 
            call. = FALSE)
    }
    fit_names <- names(fits)
    missing_names <- is_missing_string(fit_names)
    if (is.null(fit_names)) {
        fit_names <- character(length(fits))
        missing_names <- rep(TRUE, length(fits))
    }
    fit_names[missing_names] <- vapply(fit_exprs[missing_names], 
        function(expr) paste(deparse(expr, width.cutoff = 500L), 
            collapse = ""), character(1))
    fit_names[is_missing_string(fit_names)] <- paste0("fit", 
        seq_along(fits))[is_missing_string(fit_names)]
    fit_names <- make.unique(fit_names)
    list(fits = fits, fit_names = fit_names)
}
create_deexplorer_bundle <-
function (dge, ..., msigdb_genesets = NULL, msigdb_path = NULL,
    fgsea_files = NULL, sample_id_col = "Sample", gene_id_col = "ENSEMBL",
    gene_symbol_col = "SYMBOL", fdr_cutoff = 0.05, prior_count = 2)
{
    validate_dge_input(dge, sample_id_col = sample_id_col)
    fit_info <- collect_fit_inputs(...)
    for (fit_index in seq_along(fit_info$fits)) {
        validate_fit_input(fit = fit_info$fits[[fit_index]],
            fit_name = fit_info$fit_names[[fit_index]], dge = dge,
            gene_id_col = gene_id_col, gene_symbol_col = gene_symbol_col)
    }
    gene_df <- build_gene_annotation(dge = dge, gene_id_col = gene_id_col,
        gene_symbol_col = gene_symbol_col)
    gene_lookup <- build_gene_lookup(gene_df)
    sample_df <- build_sample_annotation(dge = dge, sample_id_col = sample_id_col)
    matrices <- build_expression_matrices(dge = dge, gene_df = gene_df,
        prior_count = prior_count)
    pca <- build_pca_payload(lcpm = matrices$lcpm, sample_df = sample_df)
    de_payload <- build_de_payload(fits = fit_info$fits, fit_names = fit_info$fit_names,
        gene_df = gene_df, gene_lookup = gene_lookup, fdr_cutoff = fdr_cutoff)
    msigdb_genesets <- read_msigdb_genesets(msigdb_genesets = msigdb_genesets,
        msigdb_path = msigdb_path)
    fgsea_results <- read_fgsea_results(fgsea_files = fgsea_files,
        contrast_names = de_payload$contrast_catalog$contrast)
    bundle <- list(counts = matrices$counts, cpm = matrices$cpm, 
        lcpm = matrices$lcpm, sample_df = sample_df, gene_df = gene_df, 
        gene_lookup = gene_lookup, pca_scores = pca$scores, pca_variance = pca$variance, 
        de_tables = de_payload$tables, contrast_catalog = de_payload$contrast_catalog, 
        msigdb = msigdb_genesets, msigdb_names = names(msigdb_genesets), 
        fgsea = fgsea_results, options = list(sample_id_col = sample_id_col,
            gene_id_col = gene_id_col, gene_symbol_col = gene_symbol_col,
            fdr_cutoff = fdr_cutoff, prior_count = prior_count))
    class(bundle) <- "deexplorer_bundle"
    bundle
}
deexplorer_server <-
function (bundle) 
{
    force(bundle)
    function(input, output, session) {
        stopifnot(inherits(bundle, "deexplorer_bundle"))
        first_or <- function(x, default = NULL) {
            if (length(x)) {
                x[[1L]]
            }
            else {
                default
            }
        }
        sample_metadata_choices <- setdiff(colnames(bundle$sample_df), 
            c("sample_key"))
        pca_choices <- grep("^PC[0-9]+$", colnames(bundle$pca_scores), 
            value = TRUE)
        group_default <- if ("Condition" %in% sample_metadata_choices) 
            "Condition"
        else first_or(sample_metadata_choices, "sample_id")
        pca_x_default <- first_or(pca_choices, "PC1")
        pca_y_default <- if (length(pca_choices) >= 2L) 
            pca_choices[[2L]]
        else pca_x_default
        pca_size_default <- "None"
        shiny::updateSelectInput(session, "pca_x_pc", choices = pca_choices,
            selected = pca_x_default)
        shiny::updateSelectInput(session, "pca_y_pc", choices = pca_choices, 
            selected = pca_y_default)
        shiny::updateSelectInput(session, "pca_color", choices = c("None", 
            sample_metadata_choices), selected = if ("Condition" %in% 
            sample_metadata_choices) 
            "Condition"
        else "None")
        shiny::updateSelectInput(session, "pca_shape", choices = c("None", 
            sample_metadata_choices), selected = if ("Cell_Line" %in% 
            sample_metadata_choices) 
            "Cell_Line"
        else "None")
        shiny::updateSelectInput(session, "pca_size", choices = c("None", 
            sample_metadata_choices), selected = pca_size_default)
        shiny::updateSelectInput(session, "contrast_key", choices = stats::setNames(bundle$contrast_catalog$contrast_key, 
            bundle$contrast_catalog$contrast), selected = first_or(bundle$contrast_catalog$contrast_key))
        shiny::updateSelectInput(session, "gene_box_group", choices = sample_metadata_choices, 
            selected = group_default)
        shiny::updateSelectizeInput(session, "geneset_name",
            choices = c(stats::setNames("", ""), bundle$msigdb_names),
            server = TRUE, selected = "")
        state <- shiny::reactiveValues(active_gene_key = NULL,
            hover_gene_key = NULL,
            ot_cache = new.env(parent = emptyenv()))
        current_contrast_row <- shiny::reactive({
            key <- input$contrast_key
            shiny::req(!is.null(key), nzchar(key))
            bundle$contrast_catalog[bundle$contrast_catalog$contrast_key %in% 
                key, , drop = FALSE]
        })
        current_de_table <- shiny::reactive({
            shiny::req(input$main_tabs == "Explore DEG results")
            key <- input$contrast_key
            shiny::req(!is.null(key), nzchar(key))
            de_df <- bundle$de_tables[[key]]
            if (is.null(de_df)) {
                return(data.frame())
            }
            de_df
        })
        manual_gene_keys <- shiny::reactive({
            resolve_gene_keys(identifiers = parse_multiline_gene_input(input$manual_genes), 
                gene_lookup = bundle$gene_lookup)
        })
        selected_geneset_keys <- shiny::reactive({
            geneset_name <- input$geneset_name
            if (is.null(geneset_name) || !nzchar(geneset_name)) {
                return(character())
            }
            resolve_gene_keys(bundle$msigdb[[geneset_name]],
                bundle$gene_lookup)
        })
        shiny::observeEvent(input$geneset_name, {
            geneset_name <- input$geneset_name
            if (is.null(geneset_name) || !nzchar(geneset_name)) {
                shiny::updateTextAreaInput(session, "manual_genes", value = "")
                return(invisible(NULL))
            }
            symbols <- bundle$msigdb[[geneset_name]]
            if (length(symbols)) {
                shiny::updateTextAreaInput(session, "manual_genes",
                    value = paste(symbols, collapse = "\n"))
            }
        }, ignoreNULL = FALSE)
        fgsea_hit <- shiny::reactive({
            current_contrast <- current_contrast_row()
            if (!nrow(current_contrast)) {
                return(data.frame())
            }
            resolve_fgsea_hit(bundle = bundle, contrast_name = first_or(current_contrast$contrast), 
                geneset_name = input$geneset_name)
        })
        leading_edge_keys <- shiny::reactive({
            hit <- fgsea_hit()
            if (!nrow(hit)) {
                return(character())
            }
            resolve_gene_keys(hit$leading_edge_genes[[1L]], bundle$gene_lookup)
        })
        filtered_table <- shiny::reactive({
            de_df <- current_de_table()
            if (is.null(de_df) || !nrow(de_df)) {
                return(data.frame())
            }
            geneset_keys <- selected_geneset_keys()
            if (!length(geneset_keys)) {
                return(de_df)
            }
            de_df[de_df$gene_key %in% geneset_keys, , drop = FALSE]
        })
        shiny::observeEvent(filtered_table(), {
            de_df <- filtered_table()
            if (!nrow(de_df)) {
                state$active_gene_key <- NULL
                state$hover_gene_key <- NULL
                return(invisible(NULL))
            }
            if (is.null(state$active_gene_key) || !state$active_gene_key %in% 
                de_df$gene_key) {
                state$active_gene_key <- de_df$gene_key[[1L]]
            }
            invisible(NULL)
        }, ignoreNULL = FALSE)
        shiny::observeEvent(input$de_table_rows_selected, {
            de_df <- filtered_table()
            selected_rows <- input$de_table_rows_selected
            if (!length(selected_rows)) {
                return(invisible(NULL))
            }
            selected_index <- selected_rows[[1L]]
            if (length(selected_index) && selected_index <= nrow(de_df)) {
                state$active_gene_key <- de_df$gene_key[[selected_index]]
                state$hover_gene_key <- state$active_gene_key
            }
            invisible(NULL)
        }, ignoreNULL = TRUE)
        shiny::observeEvent(input$de_table_hover, {
            if (input$de_table_hover %in% bundle$gene_df$gene_key) {
                state$hover_gene_key <- input$de_table_hover
            }
        }, ignoreNULL = TRUE)
        shiny::observeEvent(safe_event_data("plotly_click", source = "ma_plot"), 
            {
                clicked <- get_plotly_customdata(safe_event_data("plotly_click", 
                  source = "ma_plot"))
                if (length(clicked)) {
                  state$active_gene_key <- clicked[[1L]]
                  state$hover_gene_key <- clicked[[1L]]
                }
        }, ignoreNULL = TRUE)
        shiny::observeEvent(safe_event_data("plotly_click", source = "volcano_plot"), 
            {
                clicked <- get_plotly_customdata(safe_event_data("plotly_click", 
                  source = "volcano_plot"))
                if (length(clicked)) {
                  state$active_gene_key <- clicked[[1L]]
                  state$hover_gene_key <- clicked[[1L]]
                }
        }, ignoreNULL = TRUE)
        shiny::observeEvent(safe_event_data("plotly_hover", source = "ma_plot"), 
            {
                hovered <- get_plotly_customdata(safe_event_data("plotly_hover", 
                  source = "ma_plot"))
            if (length(hovered)) {
                state$hover_gene_key <- hovered[[1L]]
            }
        }, ignoreNULL = TRUE)
        shiny::observeEvent(safe_event_data("plotly_hover", source = "volcano_plot"), 
            {
                hovered <- get_plotly_customdata(safe_event_data("plotly_hover", 
                  source = "volcano_plot"))
            if (length(hovered)) {
                state$hover_gene_key <- hovered[[1L]]
            }
        }, ignoreNULL = TRUE)
        shiny::observeEvent(input$reset_manual_genes, {
            shiny::updateSelectizeInput(session, "geneset_name",
                selected = "")
            shiny::updateTextAreaInput(session, "manual_genes",
                value = "")
        })
        shiny::observeEvent(input$contrast_key, {
            state$hover_gene_key <- NULL
        })
        output$pca_plot <- plotly::renderPlotly({
            build_pca_plot(bundle = bundle, x_pc = input$pca_x_pc,
                y_pc = input$pca_y_pc, color_col = selected_choice(input$pca_color),
                shape_col = selected_choice(input$pca_shape),
                size_col = selected_choice(input$pca_size))
        }) |> shiny::bindCache(input$pca_x_pc, input$pca_y_pc,
            input$pca_color, input$pca_shape, input$pca_size)
        output$pca_variance_plot <- plotly::renderPlotly({
            build_pca_variance_plot(bundle)
        })
        output$sample_table <- DT::renderDT({
            DT::datatable(bundle$sample_df, filter = "top", selection = "single", 
                options = list(pageLength = 10, scrollX = TRUE))
        })
        output$download_counts <- shiny::downloadHandler(filename = function() paste0("deexplorer_raw_counts_", 
            Sys.Date(), ".csv"), content = function(path) download_matrix_with_annotations(bundle$counts, 
            bundle$gene_df, path))
        output$download_cpm <- shiny::downloadHandler(filename = function() paste0("deexplorer_cpm_", 
            Sys.Date(), ".csv"), content = function(path) download_matrix_with_annotations(bundle$cpm, 
            bundle$gene_df, path))
        output$download_lcpm <- shiny::downloadHandler(filename = function() paste0("deexplorer_log_cpm_", 
            Sys.Date(), ".csv"), content = function(path) download_matrix_with_annotations(bundle$lcpm, 
            bundle$gene_df, path))
        output$selection_summary <- shiny::renderUI({
            contrast_row <- current_contrast_row()
            if (!nrow(contrast_row)) {
                return(shiny::tagList(shiny::h4("Current selection"), 
                  shiny::div(class = "deexplorer-note", "Select a contrast to populate DEG exploration panels.")))
            }
            build_selection_summary(contrast_label = first_or(contrast_row$contrast,
                "Unknown contrast"),
                manual_gene_count = length(manual_gene_keys()),
                geneset_name = if (nzchar(input$geneset_name %||% 
                  "")) 
                  input$geneset_name
                else NULL, filtered_row_count = nrow(filtered_table()))
        })
        output$fgsea_summary <- shiny::renderUI({
            if (!nzchar(input$geneset_name %||% "")) {
                return(shiny::tagList(shiny::h4("Precomputed fGSEA"), 
                  shiny::div(class = "deexplorer-note", "Choose an MSigDB gene set to inspect precomputed enrichment for the active contrast.")))
            }
            build_fgsea_summary(fgsea_hit())
        })
        shiny::observeEvent(current_contrast_row(), {
            contrast_row <- current_contrast_row()
            if (!nrow(contrast_row)) return(invisible(NULL))
            contrast_name <- first_or(contrast_row$contrast, "")
            fgsea_df <- bundle$fgsea[[contrast_name]]
            collections <- if (!is.null(fgsea_df) && nrow(fgsea_df)) {
                sort(unique(fgsea_df$collection))
            } else {
                character()
            }
            labeled <- stats::setNames(collections,
                vapply(collections, format_collection_label, character(1)))
            shiny::updateSelectInput(session, "fgsea_collection",
                choices = labeled,
                selected = if (length(labeled)) labeled[[1L]] else "")
        })
        output$fgsea_barplot <- plotly::renderPlotly({
            contrast_row <- current_contrast_row()
            shiny::validate(shiny::need(nrow(contrast_row) > 0L,
                "Select a contrast to show fGSEA results."))
            contrast_name <- first_or(contrast_row$contrast, "")
            fgsea_df <- bundle$fgsea[[contrast_name]]
            collection <- input$fgsea_collection
            if (!is.null(collection) && nzchar(collection) &&
                !is.null(fgsea_df) && "collection" %in% colnames(fgsea_df)) {
                fgsea_df <- fgsea_df[fgsea_df$collection == collection, , drop = FALSE]
            }
            build_fgsea_barplot(fgsea_df)
        })
        output$de_table <- DT::renderDT({
            de_df <- filtered_table()
            shiny::validate(shiny::need(nrow(de_df) > 0L, "No genes match the current table filter."))
            display_df <- de_df[, c("gene_key", "display_symbol", 
                "gene_id", "logFC", "AveExpr", "t", "P.Value", 
                "adj.P.Val", "B", "feature_name")]
            colnames(display_df) <- c("gene_key", "Symbol", "Gene ID", 
                "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", 
                "B", "Feature")
            DT::formatSignif(DT::formatRound(DT::datatable(display_df, 
                filter = "top", selection = "single", rownames = FALSE, 
                options = list(pageLength = 15, scrollX = TRUE, 
                  columnDefs = list(list(targets = 0, visible = FALSE))), 
                callback = DT::JS(sprintf("table.on('mouseenter', 'tbody tr', function() {\n               var data = table.row(this).data();\n               if (data) { Shiny.setInputValue('%s', data[0], {priority: 'event'}); }\n             });", 
                  "de_table_hover"))), columns = c("logFC", "AveExpr", 
                "t", "B"), digits = 3), columns = c("P.Value", 
                "adj.P.Val"), digits = 3)
        })
        output$ma_plot <- plotly::renderPlotly({
            de_df <- current_de_table()
            shiny::validate(shiny::need(nrow(de_df) > 0L, "No differential expression rows are available for this contrast."))
            contrast_row <- current_contrast_row()
            shiny::validate(shiny::need(nrow(contrast_row) > 0L,
                "Select a valid contrast to draw the MA plot."))
            plotly::layout(build_de_plot(de_df = de_df, x_col = "AveExpr",
                y_col = "logFC", title = paste0("MA plot: ",
                  first_or(contrast_row$contrast, "Unknown contrast")), source_id = "ma_plot",
                highlighted_gene_keys = selected_geneset_keys(),
                manual_gene_keys = manual_gene_keys(), leading_edge_keys = leading_edge_keys(),
                active_gene_key = character(), x_title = "AveExpr",
                y_title = "logFC"), shapes = list(list(type = "line",
                x0 = min(de_df$AveExpr, na.rm = TRUE), x1 = max(de_df$AveExpr,
                  na.rm = TRUE), y0 = 0, y1 = 0, line = list(color = "#555555",
                  dash = "dash"))))
        }) |> shiny::bindCache(input$contrast_key, input$geneset_name,
            input$manual_genes)
        output$volcano_plot <- plotly::renderPlotly({
            de_df <- current_de_table()
            shiny::validate(shiny::need(nrow(de_df) > 0L, "No differential expression rows are available for this contrast."))
            contrast_row <- current_contrast_row()
            shiny::validate(shiny::need(nrow(contrast_row) > 0L,
                "Select a valid contrast to draw the volcano plot."))
            plotly::layout(build_de_plot(de_df = de_df, x_col = "logFC",
                y_col = "neg_log10_p", title = paste0("Volcano plot: ",
                  first_or(contrast_row$contrast, "Unknown contrast")), source_id = "volcano_plot",
                highlighted_gene_keys = selected_geneset_keys(),
                manual_gene_keys = manual_gene_keys(), leading_edge_keys = leading_edge_keys(),
                active_gene_key = character(), x_title = "logFC",
                y_title = "-log10(P.Value)"), shapes = list(list(type = "line",
                x0 = 0, x1 = 0, y0 = 0, y1 = max(de_df$neg_log10_p,
                  na.rm = TRUE), line = list(color = "#555555",
                  dash = "dash"))))
        }) |> shiny::bindCache(input$contrast_key, input$geneset_name,
            input$manual_genes)
        shiny::observe({
            de_df <- current_de_table()
            key <- state$active_gene_key
            ma_proxy <- plotly::plotlyProxy("ma_plot", session)
            vol_proxy <- plotly::plotlyProxy("volcano_plot", session)
            if (is.null(key) || !length(key) || !nrow(de_df) ||
                !key %in% de_df$gene_key) {
                plotly::plotlyProxyInvoke(ma_proxy, "deleteTraces", 4L)
                plotly::plotlyProxyInvoke(vol_proxy, "deleteTraces", 4L)
                return(invisible(NULL))
            }
            row <- de_df[de_df$gene_key == key, , drop = FALSE]
            active_trace <- list(list(
                x = list(row$AveExpr),
                y = list(row$logFC),
                text = list(row$hover_text),
                customdata = list(list(row$gene_key)),
                hoverinfo = "text",
                type = "scatter", mode = "markers",
                name = "Active gene",
                marker = list(color = "#0b6e4f", size = 13, symbol = "star"),
                showlegend = TRUE
            ))
            plotly::plotlyProxyInvoke(ma_proxy, "deleteTraces", 4L)
            plotly::plotlyProxyInvoke(ma_proxy, "addTraces", active_trace)
            active_trace_vol <- active_trace
            active_trace_vol[[1L]]$x <- list(row$logFC)
            active_trace_vol[[1L]]$y <- list(row$neg_log10_p)
            plotly::plotlyProxyInvoke(vol_proxy, "deleteTraces", 4L)
            plotly::plotlyProxyInvoke(vol_proxy, "addTraces", active_trace_vol)
        })
        active_gene_row <- shiny::reactive({
            de_df <- current_de_table()
            if (is.null(state$active_gene_key) || !length(state$active_gene_key) || 
                !nrow(de_df)) {
                return(data.frame())
            }
            de_df[de_df$gene_key == state$active_gene_key, , 
                drop = FALSE]
        })
        focused_gene_row <- shiny::reactive({
            de_df <- current_de_table()
            if (!nrow(de_df)) {
                return(data.frame())
            }
            focus_key <- state$hover_gene_key
            if (is.null(focus_key) || !length(focus_key) || !focus_key %in% 
                de_df$gene_key) {
                focus_key <- state$active_gene_key
            }
            if (is.null(focus_key) || !length(focus_key) || !focus_key %in% 
                de_df$gene_key) {
                return(data.frame())
            }
            de_df[de_df$gene_key == focus_key, , drop = FALSE]
        })
        output$selected_gene_summary <- shiny::renderUI({
            build_selected_gene_summary(focused_gene_row())
        })
        output$gene_boxplot <- plotly::renderPlotly({
            gene_row <- focused_gene_row()
            shiny::validate(shiny::need(nrow(gene_row) == 1L,
                "Hover or select a gene to show the per-sample log-CPM distribution."))
            gene_expr_df <- prepare_gene_expression_df(bundle,
                gene_row$gene_key[[1L]], input$gene_box_group)
            build_gene_boxplot(gene_expr_df, gene_row$display_symbol[[1L]],
                input$gene_box_group)
        }) |> shiny::bindCache(state$hover_gene_key, state$active_gene_key,
            input$gene_box_group)
        active_ensembl_id <- shiny::reactive({
            gene_row <- active_gene_row()
            if (!nrow(gene_row)) return(NULL)
            gene_id <- gene_row$gene_id[[1L]]
            if (is.null(gene_id) || !grepl("^ENSG\\d+$", gene_id)) return(NULL)
            gene_id
        })
        active_gene_info <- shiny::reactive({
            eid <- active_ensembl_id()
            if (is.null(eid)) return(NULL)
            cache_key <- paste0("info_", eid)
            if (exists(cache_key, envir = state$ot_cache, inherits = FALSE)) {
                return(get(cache_key, envir = state$ot_cache, inherits = FALSE))
            }
            result <- fetch_opentargets_gene_info(eid)
            assign(cache_key, result, envir = state$ot_cache)
            result
        })
        active_gene_assoc <- shiny::reactive({
            eid <- active_ensembl_id()
            if (is.null(eid)) return(NULL)
            cache_key <- paste0("assoc_", eid)
            if (exists(cache_key, envir = state$ot_cache, inherits = FALSE)) {
                return(get(cache_key, envir = state$ot_cache, inherits = FALSE))
            }
            result <- fetch_opentargets_associations(eid)
            assign(cache_key, result, envir = state$ot_cache)
            result
        })
        output$gene_function_text <- shiny::renderUI({
            info <- active_gene_info()
            gene_row <- active_gene_row()
            if (is.null(info) || !nrow(gene_row)) {
                return(shiny::div(class = "deexplorer-note",
                    "Click a gene to load external annotations from Open Targets."))
            }
            shiny::tagList(
                shiny::h4(paste0(info$symbol, " \u2014 ", info$name)),
                shiny::p(info$function_text,
                    style = "font-size: 0.92rem; color: #4a5568; max-height: 120px; overflow-y: auto;"),
                shiny::div(class = "deexplorer-note",
                    shiny::tags$a(
                        href = paste0("https://platform.opentargets.org/target/", active_ensembl_id()),
                        target = "_blank",
                        paste0("View ", info$symbol, " on Open Targets")
                    )
                )
            )
        })
        output$tissue_expression_plot <- plotly::renderPlotly({
            info <- active_gene_info()
            shiny::validate(shiny::need(
                !is.null(info) && nrow(info$expressions_df) > 0L,
                "Click a gene to show tissue expression."
            ))
            build_tissue_expression_barplot(info$expressions_df)
        })
        output$disease_associations_table <- DT::renderDT({
            assoc_df <- active_gene_assoc()
            shiny::validate(shiny::need(
                !is.null(assoc_df) && nrow(assoc_df) > 0L,
                "Click a gene to show disease associations."
            ))
            col_meta <- ot_datasource_labels()
            display_names <- vapply(colnames(assoc_df), function(cn) {
                entry <- col_meta[[cn]]
                if (!is.null(entry)) entry$label else cn
            }, character(1), USE.NAMES = FALSE)
            header_tips <- vapply(colnames(assoc_df), function(cn) {
                entry <- col_meta[[cn]]
                if (!is.null(entry)) entry$tip else cn
            }, character(1), USE.NAMES = FALSE)
            colnames(assoc_df) <- display_names
            numeric_cols <- setdiff(display_names, "Disease")
            header_cb <- DT::JS(sprintf(
                "function(thead, data, start, end, display) {
                    var tips = %s;
                    $(thead).find('th').each(function(i) {
                        $(this).attr('title', tips[i]);
                    });
                }",
                jsonlite::toJSON(unname(header_tips), auto_unbox = FALSE)
            ))
            DT::formatRound(
                DT::datatable(assoc_df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 10, scrollX = TRUE,
                        headerCallback = header_cb,
                        columnDefs = list(list(
                            targets = seq_along(numeric_cols),
                            render = DT::JS(
                                "function(data, type, row, meta) {",
                                "  if (type === 'display') {",
                                "    if (data === 0 || data === null) return '';",
                                "    return parseFloat(data).toFixed(2);",
                                "  }",
                                "  return data;",
                                "}"
                            )
                        ))
                    )
                ),
                columns = numeric_cols, digits = 2
            )
        })
        shiny::updateSelectInput(session, "hm_contrast",
            choices = stats::setNames(bundle$contrast_catalog$contrast_key,
                bundle$contrast_catalog$contrast),
            selected = first_or(bundle$contrast_catalog$contrast_key))
        shiny::updateCheckboxGroupInput(session, "hm_samples",
            choices = bundle$sample_df$sample_id,
            selected = bundle$sample_df$sample_id)
        hm_col_ann_choices <- setdiff(colnames(bundle$sample_df),
            c("sample_key", "sample_id", "hover_text"))
        hm_col_ann_choices <- hm_col_ann_choices[
            vapply(hm_col_ann_choices, function(col) {
                !is.numeric(bundle$sample_df[[col]])
            }, logical(1))]
        shiny::updateCheckboxGroupInput(session, "hm_col_ann",
            choices = hm_col_ann_choices)
        hm_selected_genes <- shiny::reactive({
            shiny::req(input$main_tabs == "Heatmap maker")
            manual_ids <- parse_multiline_gene_input(input$hm_manual_genes)
            if (length(manual_ids)) {
                keys <- resolve_gene_keys(manual_ids, bundle$gene_lookup)
                keys <- intersect(keys, rownames(bundle$lcpm))
                shiny::validate(shiny::need(length(keys) >= 2L,
                    "None of the entered genes were found. Try gene symbols (e.g. TP53) or Ensembl IDs."))
                return(keys)
            }
            contrast_key <- input$hm_contrast
            shiny::req(contrast_key, nzchar(contrast_key))
            de_df <- bundle$de_tables[[contrast_key]]
            shiny::req(!is.null(de_df) && nrow(de_df) > 0L)
            top_n <- input$hm_top_n
            sorted <- de_df[order(de_df$logFC, decreasing = TRUE), , drop = FALSE]
            top_up <- utils::head(sorted, top_n)
            top_down <- utils::tail(sorted, top_n)
            unique(c(top_up$gene_key, top_down$gene_key))
        })
        shiny::observe({
            selected_genes <- hm_selected_genes()
            shiny::req(length(selected_genes) > 0L)
            gene_symbols <- toupper(bundle$gene_df[selected_genes, "symbol"])
            gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]
            if (length(bundle$msigdb) > 0L && length(gene_symbols) > 0L) {
                matching <- vapply(bundle$msigdb, function(gs) {
                    any(toupper(gs) %in% gene_symbols)
                }, logical(1))
                filtered_names <- bundle$msigdb_names[matching]
            } else {
                filtered_names <- character(0)
            }
            shiny::updateSelectizeInput(session, "hm_row_genesets",
                choices = c(stats::setNames("", ""), filtered_names),
                server = TRUE)
        })
        hm_data <- shiny::reactive({
            selected_genes <- hm_selected_genes()
            shiny::validate(shiny::need(length(selected_genes) >= 2L,
                "No DE results for selected contrast."))
            selected_samples <- input$hm_samples
            if (!length(selected_samples)) {
                selected_samples <- bundle$sample_df$sample_id
            }
            sample_keys <- bundle$sample_df$sample_key[
                bundle$sample_df$sample_id %in% selected_samples]
            sample_keys <- intersect(sample_keys, colnames(bundle$lcpm))
            shiny::validate(shiny::need(length(sample_keys) >= 2L,
                "Select at least 2 samples."))
            lcpm_sub <- bundle$lcpm[selected_genes, sample_keys, drop = FALSE]
            col_annotations <- list()
            for (ann in input$hm_col_ann) {
                if (ann %in% colnames(bundle$sample_df)) {
                    col_annotations[[ann]] <- as.character(
                        bundle$sample_df[sample_keys, ann])
                }
            }
            row_geneset_flags <- list()
            for (gs in input$hm_row_genesets) {
                if (!nzchar(gs)) next
                gs_genes <- bundle$msigdb[[gs]]
                if (length(gs_genes)) {
                    gene_symbols <- bundle$gene_df[selected_genes, "symbol"]
                    flag <- ifelse(toupper(gene_symbols) %in% toupper(gs_genes),
                        "In", "Out")
                    short_name <- sub("^HALLMARK_", "", gs)
                    short_name <- gsub("_", " ", short_name)
                    row_geneset_flags[[short_name]] <- flag
                }
            }
            list(lcpm = lcpm_sub, sample_keys = sample_keys,
                selected_genes = selected_genes,
                col_annotations = col_annotations,
                row_geneset_flags = row_geneset_flags)
        })
        output$iheatmap_plot <- iheatmapr::renderIheatmap({
            d <- hm_data()
            shiny::req(d)
            cluster_rows <- if ("rows" %in% input$hm_cluster) "hclust" else "none"
            cluster_cols <- if ("columns" %in% input$hm_cluster) "hclust" else "none"
            build_iheatmap(
                lcpm_matrix = d$lcpm,
                gene_df = bundle$gene_df,
                sample_df = bundle$sample_df,
                col_annotations = d$col_annotations,
                row_geneset_flags = d$row_geneset_flags,
                cluster_rows = cluster_rows,
                cluster_cols = cluster_cols,
                show_row_labels = input$hm_show_row_labels,
                show_col_labels = input$hm_show_col_labels)
        }) |> shiny::bindCache(input$hm_contrast, input$hm_top_n,
            input$hm_manual_genes, input$hm_samples, input$hm_col_ann,
            input$hm_row_genesets, input$hm_cluster,
            input$hm_show_row_labels, input$hm_show_col_labels)
    }
}
deexplorer_styles <-
function () 
{
    shiny::tags$style(shiny::HTML(paste(".deexplorer-title { letter-spacing: 0.02em; }", 
        ".deexplorer-sidebar { background: rgba(255, 255, 255, 0.78); border: 1px solid rgba(31, 41, 51, 0.10); border-radius: 12px; padding: 16px; }", 
        ".deexplorer-card { background: rgba(255, 255, 255, 0.92); border: 1px solid rgba(31, 41, 51, 0.10); border-radius: 12px; padding: 12px; margin-bottom: 16px; box-shadow: 0 12px 28px rgba(31, 41, 51, 0.06); }", 
        ".deexplorer-pca-sidebar .form-group { margin-bottom: 8px; }", 
        ".deexplorer-pca-sidebar .control-label { margin-bottom: 2px; font-size: 0.92rem; }", 
        ".deexplorer-pca-sidebar .form-control { height: 32px; padding: 4px 10px; font-size: 0.92rem; }", 
        ".deexplorer-pca-sidebar hr { margin-top: 10px; margin-bottom: 10px; }", 
        ".deexplorer-download-stack .btn { width: 100%; margin-bottom: 8px; }", 
        ".deexplorer-download-stack .btn:last-child { margin-bottom: 0; }", 
        ".deexplorer-note { font-size: 0.92rem; color: #4a5568; }", 
        ".deexplorer-metric { margin-bottom: 10px; }", ".tab-content { padding-top: 14px; }", 
        ".dataTables_wrapper .dataTables_filter input { width: 240px; }", 
        sep = "\n")))
}
deexplorer_theme <-
function () 
{
    bslib::bs_theme(version = 5, bg = "#f4f1ea", fg = "#1f2933", 
        primary = "#0b6e4f", secondary = "#7a4f01", base_font = bslib::font_google("Source Sans 3"), 
        heading_font = bslib::font_google("Merriweather Sans"))
}
deexplorer_ui <-
function (title = "DEExplorer") 
{
    shiny::fluidPage(theme = deexplorer_theme(), deexplorer_styles(), 
        shiny::div(class = "deexplorer-title", shiny::titlePanel(title)), 
        shiny::div(class = "deexplorer-note", "This app assumes the supplied DGEList has already been filtered to the analysis universe and that limma fits were computed on the same genes."), 
        shiny::tabsetPanel(id = "main_tabs", pca_tab_ui(), deg_tab_ui(), heatmap_tab_ui()))
}
deg_tab_ui <-
function ()
{
    shiny::tabPanel(title = "Explore DEG results", shiny::fluidRow(shiny::column(width = 3,
        shiny::div(class = "deexplorer-sidebar", shiny::h4("Contrast and Highlighting"),
            shiny::selectInput("contrast_key", "Contrast", choices = character()),
            shiny::selectInput("gene_box_group", "Group selected-gene boxplot by",
                choices = character()), shiny::selectizeInput("geneset_name",
                "MSigDB gene set", choices = NULL, selected = "",
                options = list(
                    placeholder = "Type a Hallmark, GO or pathway",
                    respect_word_boundaries = FALSE)),
            shiny::textAreaInput("manual_genes", "Highlight genes (one per line)",
                rows = 8, placeholder = "STAT1\nCXCL8\nTFRC"),
            shiny::actionButton("reset_manual_genes", "Clear gene list"),
            shiny::hr(), shiny::uiOutput("selection_summary"),
            shiny::uiOutput("fgsea_summary"),
            shiny::hr(),
            shiny::selectInput("fgsea_collection", "fGSEA collection",
                choices = character()))), shiny::column(width = 9,
        shiny::fluidRow(shiny::column(width = 6, shiny::div(class = "deexplorer-card",
            plotly::plotlyOutput("ma_plot", height = "410px"))),
            shiny::column(width = 6, shiny::div(class = "deexplorer-card",
                plotly::plotlyOutput("volcano_plot", height = "410px")))),
        shiny::fluidRow(shiny::column(width = 7, shiny::div(class = "deexplorer-card",
            shiny::h4("Ranked Differential Expression Table"),
            DT::DTOutput("de_table"))), shiny::column(width = 5,
            shiny::div(class = "deexplorer-card", shiny::uiOutput("selected_gene_summary")),
            shiny::div(class = "deexplorer-card", plotly::plotlyOutput("gene_boxplot",
                height = "310px")))),
        shiny::fluidRow(shiny::column(width = 12,
            shiny::div(class = "deexplorer-card",
                plotly::plotlyOutput("fgsea_barplot", height = "380px")))),
        shiny::fluidRow(shiny::column(width = 12,
            shiny::div(class = "deexplorer-card", shiny::uiOutput("gene_function_text")))),
        shiny::fluidRow(shiny::column(width = 4,
            shiny::div(class = "deexplorer-card",
                plotly::plotlyOutput("tissue_expression_plot", height = "420px"))),
            shiny::column(width = 8,
                shiny::div(class = "deexplorer-card",
                    shiny::h4("Disease Associations (Open Targets)"),
                    DT::DTOutput("disease_associations_table")))))))
}
download_matrix_with_annotations <-
function (matrix_data, gene_df, path) 
{
    download_df <- cbind(gene_df[, setdiff(colnames(gene_df), 
        c("symbol_upper", "gene_id_upper", "feature_upper")), 
        drop = FALSE], as.data.frame(matrix_data[gene_df$gene_key, 
        , drop = FALSE], check.names = FALSE))
    utils::write.csv(download_df, path, row.names = FALSE)
}
format_collection_label <-
function (collection)
{
    known <- list(
        hallmark = "Hallmarks",
        gobp = "GO:BP",
        gocc = "GO:CC",
        gomf = "GO:MF",
        kegg = "KEGG",
        reactome = "Reactome",
        biocarta = "BioCarta",
        wp = "WikiPathways",
        pid = "PID",
        tft = "TFT",
        mir = "miRNA Targets",
        hpo = "HPO",
        positional_x = "Positional (X)",
        positional = "Positional"
    )
    label <- known[[tolower(collection)]]
    if (!is.null(label)) return(label)
    if (grepl("^positional_", collection, ignore.case = TRUE)) {
        chr <- sub("^positional_", "", collection, ignore.case = TRUE)
        return(paste0("Positional (", toupper(chr), ")"))
    }
    toupper(collection)
}
format_number <-
function (x, digits = 3)
{
    formatC(x, digits = digits, format = "fg", flag = "#")
}
format_pvalue <-
function (x)
{
    format(x, digits = 3, scientific = TRUE)
}
fetch_opentargets_associations <-
function (ensembl_id, size = 20L)
{
    if (is.null(ensembl_id) || !grepl("^ENSG\\d+$", ensembl_id)) {
        return(NULL)
    }
    query <- '
    query TargetAssociations($id: String!) {
      target(ensemblId: $id) {
        associatedDiseases(
          page: { index: 0, size: 20 }
          orderByScore: "score"
        ) {
          count
          rows {
            disease { id name }
            score
            datasourceScores { componentId: id score }
          }
        }
      }
    }'
    tryCatch({
        resp <- httr2::req_perform(
            httr2::req_body_json(
                httr2::req_url_path_append(
                    httr2::request("https://api.platform.opentargets.org"),
                    "api", "v4", "graphql"
                ),
                data = list(query = query, variables = list(id = ensembl_id))
            )
        )
        payload <- httr2::resp_body_json(resp, simplifyVector = FALSE)
        rows <- payload$data$target$associatedDiseases$rows
        if (!length(rows)) return(NULL)
        all_source_ids <- unique(unlist(lapply(rows, function(r) {
            vapply(r$datasourceScores, function(s) s$componentId, character(1))
        })))
        parsed <- lapply(rows, function(r) {
            base <- data.frame(
                Disease = r$disease$name %||% "",
                Score = r$score %||% 0,
                stringsAsFactors = FALSE
            )
            scores <- stats::setNames(
                vapply(r$datasourceScores, function(s) s$score %||% 0, numeric(1)),
                vapply(r$datasourceScores, function(s) s$componentId, character(1))
            )
            for (src in all_source_ids) {
                base[[src]] <- unname(scores[src]) %||% NA_real_
            }
            base
        })
        result <- do.call(rbind, parsed)
        result[is.na(result)] <- 0
        result
    }, error = function(e) NULL)
}
fetch_opentargets_gene_info <-
function (ensembl_id)
{
    if (is.null(ensembl_id) || !grepl("^ENSG\\d+$", ensembl_id)) {
        return(NULL)
    }
    query <- '
    query TargetInfo($id: String!) {
      target(ensemblId: $id) {
        approvedSymbol
        approvedName
        functionDescriptions
        expressions {
          tissue { label organs }
          rna { level value unit }
        }
      }
    }'
    tryCatch({
        resp <- httr2::req_perform(
            httr2::req_body_json(
                httr2::req_url_path_append(
                    httr2::request("https://api.platform.opentargets.org"),
                    "api", "v4", "graphql"
                ),
                data = list(query = query, variables = list(id = ensembl_id))
            )
        )
        payload <- httr2::resp_body_json(resp, simplifyVector = FALSE)
        target <- payload$data$target
        if (is.null(target)) return(NULL)
        func_texts <- target$functionDescriptions
        func_text <- if (length(func_texts)) {
            sub("\\s*\\{ECO:.*$", "", func_texts[[1L]])
        } else ""
        expr_rows <- lapply(target$expressions, function(e) {
            organs <- e$tissue$organs
            if (!length(organs)) organs <- list(e$tissue$label)
            data.frame(
                tissue = e$tissue$label %||% "",
                organ = organs[[1L]] %||% "",
                rna_value = e$rna$value %||% 0,
                rna_level = e$rna$level %||% -1L,
                stringsAsFactors = FALSE
            )
        })
        expr_df <- do.call(rbind, expr_rows)
        if (!is.null(expr_df) && nrow(expr_df)) {
            expr_df <- expr_df[expr_df$rna_value > 0, , drop = FALSE]
            organ_agg <- stats::aggregate(
                rna_value ~ organ, data = expr_df, FUN = max
            )
            organ_agg <- organ_agg[order(-organ_agg$rna_value), , drop = FALSE]
        } else {
            organ_agg <- data.frame(organ = character(), rna_value = numeric())
        }
        list(
            symbol = target$approvedSymbol %||% "",
            name = target$approvedName %||% "",
            function_text = func_text,
            expressions_df = organ_agg
        )
    }, error = function(e) NULL)
}
get_plotly_customdata <-
function (event) 
{
    if (is.null(event) || is.null(event$customdata)) {
        return(character())
    }
    customdata <- event$customdata
    if (is.matrix(customdata) || is.data.frame(customdata)) {
        if (!nrow(customdata) || !ncol(customdata)) {
            return(character())
        }
        values <- customdata[, 1L]
    }
    else if (is.list(customdata)) {
        values <- vapply(customdata, function(item) {
            if (length(item) > 1L) {
                as.character(item[[1L]])
            }
            else {
                as.character(item)
            }
        }, character(1))
    }
    else {
        values <- as.character(customdata)
    }
    unique(values[nzchar(values)])
}
safe_event_data <-
function (event, source) 
{
    tryCatch(suppressWarnings(plotly::event_data(event, source = source)),
        error = function(...) NULL)
}
heatmap_tab_ui <-
function ()
{
    shiny::tabPanel(title = "Heatmap maker", shiny::fluidRow(shiny::column(width = 3,
        shiny::div(class = "deexplorer-sidebar", shiny::h4("Heatmap Settings"),
            shiny::selectInput("hm_contrast", "Contrast", choices = character()),
            shiny::sliderInput("hm_top_n", "Top up/down genes",
                min = 10, max = 100, value = 20, step = 5),
            shiny::textAreaInput("hm_manual_genes", "Custom gene list (one per line)",
                value = "", rows = 4,
                placeholder = "Leave empty to use top up/down genes"),
            shiny::checkboxGroupInput("hm_cluster", "Clustering",
                choices = c("Rows" = "rows", "Columns" = "columns"),
                inline = TRUE),
            shiny::checkboxInput("hm_show_row_labels", "Show row labels", TRUE),
            shiny::checkboxInput("hm_show_col_labels", "Show column labels", TRUE),
            shiny::hr(),
            shiny::h4("Columns"),
            shiny::tags$details(
                shiny::tags$summary(shiny::tags$strong("Samples to show")),
                shiny::checkboxGroupInput("hm_samples", NULL,
                    choices = character())),
            shiny::hr(),
            shiny::h4("Annotations"),
            shiny::tags$details(
                shiny::tags$summary(shiny::tags$strong("Column annotations (metadata)")),
                shiny::checkboxGroupInput("hm_col_ann", NULL,
                    choices = character())),
            shiny::selectizeInput("hm_row_genesets",
                "Row annotations (MSigDB gene sets)",
                choices = NULL, multiple = TRUE,
                options = list(
                    placeholder = "Type a gene set name",
                    respect_word_boundaries = FALSE)))),
        shiny::column(width = 9,
            shiny::div(class = "deexplorer-card",
                iheatmapr::iheatmaprOutput("iheatmap_plot",
                    height = "700px")))))
}
infer_fit_method <-
function (fit) 
{
    if (!is.null(fit$treat.lfc)) 
        "treat"
    else "ebayes"
}
is_missing_string <-
function (x) 
{
    is.na(x) | !nzchar(trimws(as.character(x)))
}
make_gene_hover_text <-
function (de_df) 
{
    paste0("<b>", de_df$display_symbol, "</b><br>", "Gene ID: ", 
        de_df$gene_id, "<br>", "Contrast: ", de_df$contrast, 
        "<br>", "Method: ", de_df$method_label, "<br>", "logFC: ", 
        formatC(de_df$logFC, digits = 3, format = "fg"), "<br>", 
        "AveExpr: ", formatC(de_df$AveExpr, digits = 3, format = "fg"), 
        "<br>", "P.Value: ", format(de_df$P.Value, scientific = TRUE, 
            digits = 3), "<br>", "adj.P.Val: ", format(de_df$adj.P.Val, 
            scientific = TRUE, digits = 3), "<br>", "t: ", formatC(de_df$t, 
            digits = 3, format = "fg"))
}
make_sample_hover_text <-
function (sample_df) 
{
    apply(sample_df, 1L, function(row_values) {
        labels <- paste(names(row_values), row_values, sep = ": ")
        paste0("<b>", row_values[["sample_id"]], "</b><br>", 
            paste(labels, collapse = "<br>"))
    })
}
ot_datasource_labels <-
function ()
{
    list(
        Disease              = list(label = "Disease",              tip = "Disease or phenotype name"),
        Score                = list(label = "Score",                tip = "Overall association score (0\u20131)"),
        ot_genetics_portal   = list(label = "OT Genetics",         tip = "GWAS locus-to-gene pipeline (Open Targets Genetics)"),
        gwas_credible_sets   = list(label = "GWAS Credible Sets",  tip = "Fine-mapped credible sets from curated GWAS studies"),
        eva                  = list(label = "EVA",                 tip = "ClinVar germline variant\u2013disease associations (EBI EVA)"),
        gene_burden          = list(label = "Gene Burden",         tip = "Rare-variant gene burden tests from exome/genome studies"),
        genomics_england     = list(label = "Genomics England",    tip = "Genomics England PanelApp gene\u2013disease panels"),
        gene2phenotype       = list(label = "Gene2Phenotype",      tip = "Curated developmental disorder gene\u2013phenotype pairs"),
        uniprot_literature   = list(label = "UniProt Literature",  tip = "UniProt literature-curated variant annotations"),
        uniprot_variants     = list(label = "UniProt Variants",    tip = "UniProt reviewed variant\u2013disease mappings"),
        orphanet             = list(label = "Orphanet",            tip = "Orphanet rare disease\u2013gene associations"),
        clingen              = list(label = "ClinGen",             tip = "ClinGen gene\u2013disease clinical validity curations"),
        cancer_gene_census   = list(label = "Cancer Gene Census",  tip = "COSMIC Cancer Gene Census somatic driver genes"),
        intogen              = list(label = "IntOGen",             tip = "IntOGen cancer driver gene identification"),
        eva_somatic          = list(label = "EVA Somatic",         tip = "ClinVar somatic variant\u2013disease associations"),
        cancer_biomarkers    = list(label = "Cancer Biomarkers",   tip = "Cancer Genome Interpreter biomarker\u2013drug\u2013disease links"),
        chembl               = list(label = "ChEMBL",             tip = "ChEMBL drug mechanism of action evidence"),
        crispr_screen        = list(label = "CRISPR Screen",       tip = "Project Score genome-wide CRISPR fitness screens"),
        crispr               = list(label = "CRISPR",              tip = "CRISPR functional genomics screens"),
        slapenrich           = list(label = "SLAPenrich",          tip = "SLAPenrich pathway-level somatic mutation enrichment"),
        progeny              = list(label = "PROGENy",             tip = "PROGENy pathway activity from perturbation signatures"),
        reactome             = list(label = "Reactome",            tip = "Reactome curated pathway associations"),
        sysbio               = list(label = "SysBio",              tip = "Systems biology literature-curated functional evidence"),
        europepmc            = list(label = "Europe PMC",          tip = "Europe PMC text-mined co-occurrences"),
        expression_atlas     = list(label = "Expression Atlas",    tip = "Expression Atlas differential expression evidence"),
        impc                 = list(label = "IMPC",                tip = "IMPC mouse model phenotype\u2013human disease associations"),
        ot_crispr_validation = list(label = "OT CRISPR Validation", tip = "Open Targets CRISPR validation experiments")
    )
}
map_table_gene_keys <-
function (table_df, gene_lookup) 
{
    mapped <- rep(NA, nrow(table_df))
    fill_from_lookup <- function(values, lookup_map) {
        if (is.null(values)) {
            return(invisible(NULL))
        }
        values <- toupper(as.character(values))
        keep <- is.na(mapped) & !is_missing_string(values)
        if (!any(keep)) {
            return(invisible(NULL))
        }
        matched_values <- unname(lookup_map[values[keep]])
        mapped[keep] <<- matched_values
        invisible(NULL)
    }
    fill_from_lookup(rownames(table_df), gene_lookup$by_key)
    fill_from_lookup(table_df[["gene_id"]], gene_lookup$first_by_gene_id)
    fill_from_lookup(table_df[["ENSEMBL"]], gene_lookup$first_by_gene_id)
    fill_from_lookup(table_df[["symbol"]], gene_lookup$first_by_symbol)
    fill_from_lookup(table_df[["SYMBOL"]], gene_lookup$first_by_symbol)
    fill_from_lookup(table_df[["feature_name"]], gene_lookup$first_by_feature)
    fill_from_lookup(rownames(table_df), gene_lookup$first_by_feature)
    fill_from_lookup(rownames(table_df), gene_lookup$first_by_symbol)
    if (anyNA(mapped)) {
        stop(sprintf("Failed to match %d rows from a limma topTable/topTreat result back to `dge` gene annotations.", 
            sum(is.na(mapped))), call. = FALSE)
    }
    mapped
}
match_fgsea_contrast <-
function (file_stub, contrast_names) 
{
    matches <- contrast_names[endsWith(file_stub, paste0("_", 
        contrast_names))]
    if (!length(matches)) {
        return(NA)
    }
    matches[[which.max(nchar(matches))]]
}
parse_fgsea_leading_edge <-
function (x) 
{
    pieces <- unique(trimws(unlist(strsplit(as.character(x), 
        "\\|", perl = TRUE))))
    pieces[nzchar(pieces)]
}
parse_multiline_gene_input <-
function (text) 
{
    if (is.null(text) || !nzchar(text)) {
        return(character())
    }
    genes <- unique(trimws(unlist(strsplit(text, "\\r?\\n", perl = TRUE))))
    genes[nzchar(genes)]
}
pca_tab_ui <-
function () 
{
    shiny::tabPanel(title = "PCA", shiny::fluidRow(shiny::column(width = 3, 
        shiny::div(class = "deexplorer-sidebar deexplorer-pca-sidebar", 
            shiny::h4("Sample Embedding"), shiny::selectInput("pca_x_pc", 
                "X-axis PC", choices = character()), shiny::selectInput("pca_y_pc", 
                "Y-axis PC", choices = character()), shiny::selectInput("pca_color", 
                "Color by metadata", choices = "None"), shiny::selectInput("pca_shape", 
                "Shape by metadata", choices = "None"), shiny::selectInput("pca_size", 
                "Size by metadata", choices = "None"), shiny::hr(), 
            shiny::div(class = "deexplorer-download-stack", shiny::downloadButton("download_counts", 
                "Download raw counts"), shiny::downloadButton("download_cpm", 
                "Download CPM"), shiny::downloadButton("download_lcpm", 
                "Download log-CPM")), shiny::hr(), shiny::div(class = "deexplorer-note", 
                "Hover tooltips include all columns available in `dge$samples`."))), 
        shiny::column(width = 9, shiny::fluidRow(shiny::column(width = 8, 
            shiny::div(class = "deexplorer-card", plotly::plotlyOutput("pca_plot", 
                height = "430px"))), shiny::column(width = 4, 
            shiny::div(class = "deexplorer-card", plotly::plotlyOutput("pca_variance_plot", 
                height = "430px")))), shiny::fluidRow(shiny::column(width = 12, 
            shiny::div(class = "deexplorer-card", shiny::h4("Sample Metadata"), 
                DT::DTOutput("sample_table")))))))
}
prepare_gene_expression_df <-
function (bundle, gene_key, group_col) 
{
    if (is.null(gene_key) || !length(gene_key) || !gene_key %in% 
        rownames(bundle$lcpm)) {
        return(data.frame())
    }
    sample_df <- bundle$sample_df
    group_values <- if (!is.null(group_col) && group_col %in% 
        colnames(sample_df)) {
        sample_df[colnames(bundle$lcpm), group_col]
    }
    else {
        sample_df[colnames(bundle$lcpm), "sample_id"]
    }
    expr_df <- data.frame(sample_key = colnames(bundle$lcpm),
        sample_id = sample_df[colnames(bundle$lcpm), "sample_id"],
        group_label = as.character(group_values),
        log_cpm = as.numeric(bundle$lcpm[gene_key, ]),
        stringsAsFactors = FALSE)
    expr_df$hover_text <- paste0("<b>", expr_df$sample_id, "</b><br>", 
        group_col %||% "Sample", ": ", expr_df$group_label, "<br>", 
        "log-CPM: ", format_number(expr_df$log_cpm))
    expr_df
}
read_fgsea_results <-
function (fgsea_files, contrast_names)
{
    if (is.null(fgsea_files) || !length(fgsea_files)) {
        return(list())
    }
    contrast_names <- unique(as.character(contrast_names))
    output <- list()
    for (file_path in fgsea_files) {
        file_stub <- sub("\\.tsv$", "", basename(file_path))
        contrast_name <- match_fgsea_contrast(file_stub, contrast_names)
        if (is.na(contrast_name)) {
            next
        }
        trimmed_stub <- sub("^fGSEA_", "", file_stub)
        collection_name <- substr(trimmed_stub, start = 1L, stop = nchar(trimmed_stub) - 
            nchar(contrast_name) - 1L)
        fgsea_df <- utils::read.delim(file_path, stringsAsFactors = FALSE, 
            check.names = FALSE)
        fgsea_df$collection <- collection_name
        fgsea_df$contrast <- contrast_name
        fgsea_df$leading_edge_genes <- lapply(fgsea_df$leadingEdge, 
            parse_fgsea_leading_edge)
        if (is.null(output[[contrast_name]])) {
            output[[contrast_name]] <- fgsea_df
        }
        else {
            output[[contrast_name]] <- rbind(output[[contrast_name]], 
                fgsea_df)
        }
    }
    output
}
read_msigdb_genesets <-
function (msigdb_genesets, msigdb_path) 
{
    if (is.null(msigdb_genesets) && !is.null(msigdb_path)) {
        msigdb_genesets <- readRDS(msigdb_path)
    }
    standardize_msigdb_genesets(msigdb_genesets)
}
resolve_fgsea_hit <-
function (bundle, contrast_name, geneset_name) 
{
    if (is.null(contrast_name) || is.null(geneset_name) || !nzchar(geneset_name)) {
        return(data.frame())
    }
    fgsea_df <- bundle$fgsea[[contrast_name]]
    if (is.null(fgsea_df) || !nrow(fgsea_df)) {
        return(data.frame())
    }
    fgsea_df[fgsea_df$pathway == geneset_name, , drop = FALSE]
}
resolve_gene_keys <-
function (identifiers, gene_lookup) 
{
    identifiers <- unique(toupper(trimws(as.character(identifiers))))
    identifiers <- identifiers[nzchar(identifiers)]
    if (!length(identifiers)) {
        return(character())
    }
    resolved <- character()
    for (identifier in identifiers) {
        key_hit <- unname(gene_lookup$by_key[identifier])
        if (is.na(key_hit)) key_hit <- character()
        resolved <- c(resolved, key_hit,
            gene_lookup$by_gene_id[[identifier]] %||% character(),
            gene_lookup$by_symbol[[identifier]] %||% character(),
            gene_lookup$by_feature[[identifier]] %||% character())
    }
    unique(resolved[nzchar(resolved)])
}
run_deexplorer_app <-
function (dge, ..., title = "DEExplorer", msigdb_genesets = NULL,
    msigdb_path = NULL, fgsea_files = NULL, sample_id_col = "Sample",
    gene_id_col = "ENSEMBL", gene_symbol_col = "SYMBOL", fdr_cutoff = 0.05,
    prior_count = 2, launch.browser = interactive())
{
    bundle <- create_deexplorer_bundle(dge = dge, ..., msigdb_genesets = msigdb_genesets,
        msigdb_path = msigdb_path, fgsea_files = fgsea_files,
        sample_id_col = sample_id_col, gene_id_col = gene_id_col,
        gene_symbol_col = gene_symbol_col, fdr_cutoff = fdr_cutoff,
        prior_count = prior_count)
    shiny::shinyApp(ui = deexplorer_ui(title = title), server = deexplorer_server(bundle), 
        options = list(launch.browser = launch.browser))
}
scale_marker_color <-
function (values) 
{
    if (is.null(values) || !length(values)) {
        return(rep("#0b6e4f", length(values)))
    }
    if (is.numeric(values)) {
        palette <- (grDevices::colorRampPalette(c("#e8f3f1", 
            "#0b6e4f")))(100)
        finite_values <- values
        finite_values[!is.finite(finite_values)] <- stats::median(values[is.finite(values)], 
            na.rm = TRUE)
        scaled <- round(1 + 99 * (finite_values - min(finite_values, 
            na.rm = TRUE))/max(diff(range(finite_values, na.rm = TRUE)), 
            1e-12))
        scaled[!is.finite(scaled)] <- 50
        return(palette[pmax(1, pmin(100, scaled))])
    }
    values <- as.factor(values)
    palette <- c("#0b6e4f", "#b56576", "#355070", "#7a4f01", 
        "#4d908e", "#d1495b", "#7b2cbf", "#4361ee", "#ff7f11", 
        "#2a9d8f", "#6a994e", "#c1121f")
    mapped <- palette[((as.integer(values) - 1L)%%length(palette)) + 
        1L]
    mapped[is.na(mapped)] <- "#0b6e4f"
    mapped
}
scale_marker_size <-
function (values) 
{
    if (is.null(values)) {
        return(rep(14, 0L))
    }
    if (!length(values)) {
        return(numeric())
    }
    if (is.numeric(values)) {
        rng <- range(values, na.rm = TRUE)
        if (!all(is.finite(rng)) || identical(rng[[1]], rng[[2]])) {
            return(rep(14, length(values)))
        }
        return(10 + (values - rng[[1]])/diff(rng) * 12)
    }
    values <- as.factor(values)
    levels_out <- seq(10, 22, length.out = max(2L, nlevels(values)))
    sizes <- levels_out[as.integer(values)]
    sizes[is.na(sizes)] <- 14
    sizes
}
scale_marker_symbol <-
function (values) 
{
    symbol_pool <- c("circle", "square", "diamond", "triangle-up", 
        "triangle-down", "cross", "x", "star", "hexagram", "triangle-left", 
        "triangle-right")
    if (is.null(values) || !length(values)) {
        return(character())
    }
    values <- as.factor(values)
    mapped <- symbol_pool[((as.integer(values) - 1L)%%length(symbol_pool)) + 
        1L]
    mapped[is.na(mapped)] <- "circle"
    mapped
}
selected_choice <-
function (x) 
{
    if (is.null(x) || !nzchar(x) || identical(x, "None")) 
        NULL
    else x
}
standardize_msigdb_genesets <-
function (msigdb_genesets) 
{
    if (is.null(msigdb_genesets)) {
        return(list())
    }
    if (is.data.frame(msigdb_genesets)) {
        possible_name_cols <- c("gs_name", "geneset", "pathway", 
            "name")
        possible_gene_cols <- c("gene_symbol", "gene", "symbol")
        name_col <- possible_name_cols[possible_name_cols %in% 
            colnames(msigdb_genesets)][1]
        gene_col <- possible_gene_cols[possible_gene_cols %in% 
            colnames(msigdb_genesets)][1]
        if (is.na(name_col) || is.na(gene_col)) {
            stop("When `msigdb_genesets` is a data frame it must contain a gene-set name column and a gene symbol column.", 
                call. = FALSE)
        }
        msigdb_genesets <- split(msigdb_genesets[[gene_col]], 
            msigdb_genesets[[name_col]])
    }
    if (!is.list(msigdb_genesets) || is.null(names(msigdb_genesets))) {
        stop("`msigdb_genesets` must be a named list or data frame.", 
            call. = FALSE)
    }
    genesets <- lapply(msigdb_genesets, function(genes) {
        genes <- unique(toupper(trimws(as.character(genes))))
        genes[nzchar(genes)]
    })
    genesets[order(names(genesets))]
}
string_literal <-
function (x) 
{
    encodeString(x, quote = "\"")
}
subset_plot_layers <-
function (de_df, highlighted_gene_keys, manual_gene_keys, leading_edge_keys, 
    active_gene_key) 
{
    highlighted_gene_keys <- intersect(highlighted_gene_keys, 
        de_df$gene_key)
    manual_gene_keys <- intersect(manual_gene_keys, de_df$gene_key)
    leading_edge_keys <- intersect(leading_edge_keys, de_df$gene_key)
    geneset_keys <- setdiff(highlighted_gene_keys, leading_edge_keys)
    manual_only_keys <- setdiff(manual_gene_keys, c(leading_edge_keys, highlighted_gene_keys))
    list(nonsig = de_df[!de_df$significant, , drop = FALSE],
        sig = de_df[de_df$significant, , drop = FALSE], geneset = de_df[de_df$gene_key %in%
            geneset_keys, , drop = FALSE], leading = de_df[de_df$gene_key %in%
            leading_edge_keys, , drop = FALSE], manual = de_df[de_df$gene_key %in%
            manual_only_keys, , drop = FALSE],
        active = de_df[de_df$gene_key %in% active_gene_key, ,
            drop = FALSE])
}
validate_dge_input <-
function (dge, sample_id_col) 
{
    if (!inherits(dge, "DGEList")) {
        stop("`dge` must inherit from edgeR::DGEList.", call. = FALSE)
    }
    if (is.null(dge$counts) || !is.matrix(dge$counts)) {
        stop("`dge$counts` must be a matrix.", call. = FALSE)
    }
    if (is.null(dge$samples) || !nrow(dge$samples)) {
        stop("`dge$samples` must contain sample metadata.", call. = FALSE)
    }
    if (ncol(dge$counts) != nrow(dge$samples)) {
        stop("`dge$counts` columns and `dge$samples` rows must align.", 
            call. = FALSE)
    }
    if (!is.null(sample_id_col) && !sample_id_col %in% colnames(dge$samples)) {
        warning(sprintf("Sample metadata column `%s` was not found; sample IDs will fall back to count matrix column names.", 
            sample_id_col), call. = FALSE)
    }
    invisible(TRUE)
}
validate_fit_input <-
function (fit, fit_name, dge, gene_id_col, gene_symbol_col) 
{
    if (!inherits(fit, "MArrayLM")) {
        stop(sprintf("`%s` must inherit from limma::MArrayLM.", 
            fit_name), call. = FALSE)
    }
    if (is.null(fit$coefficients) || !ncol(fit$coefficients)) {
        stop(sprintf("`%s` does not contain any contrasts in `coefficients`.", 
            fit_name), call. = FALSE)
    }
    if (nrow(fit$coefficients) != nrow(dge$counts)) {
        stop(sprintf("`%s` must contain the same number of genes as `dge`.", 
            fit_name), call. = FALSE)
    }
    if (!is.null(fit$genes) && nrow(fit$genes) != nrow(dge$counts)) {
        stop(sprintf("`%s$genes` must contain the same number of rows as `dge$counts`.", 
            fit_name), call. = FALSE)
    }
    has_gene_id <- !is.null(fit$genes) && gene_id_col %in% colnames(fit$genes)
    has_symbol <- !is.null(fit$genes) && gene_symbol_col %in% 
        colnames(fit$genes)
    if (!has_gene_id && !has_symbol && is.null(rownames(fit$coefficients))) {
        stop(sprintf("`%s` must expose gene identifiers through `%s`, `%s`, or row names.", 
            fit_name, gene_id_col, gene_symbol_col), call. = FALSE)
    }
    invisible(TRUE)
}
write_app_wrapper_files <-
function (app_dir, title) 
{
    ui_lines <- c(sprintf("DEExplorer::deexplorer_ui(title = %s)", 
        string_literal(title)))
    server_lines <- c("bundle <- readRDS(\"app-data.rds\")", "", 
        "DEExplorer::deexplorer_server(bundle)")
    writeLines(ui_lines, con = file.path(app_dir, "ui.R"))
    writeLines(server_lines, con = file.path(app_dir, "server.R"))
}
write_deexplorer_app <-
function (dge, ..., app_dir = "inst/shiny/deexplorer", title = "DEExplorer",
    overwrite = FALSE, launch = FALSE, msigdb_genesets = NULL,
    msigdb_path = NULL, fgsea_files = NULL, sample_id_col = "Sample",
    gene_id_col = "ENSEMBL", gene_symbol_col = "SYMBOL", fdr_cutoff = 0.05,
    prior_count = 2)
{
    app_dir <- normalizePath(app_dir, winslash = "/", mustWork = FALSE)
    dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)
    target_files <- file.path(app_dir, c("ui.R", "server.R",
        "app-data.rds"))
    existing_files <- file.exists(target_files)
    if (any(existing_files) && !isTRUE(overwrite)) {
        stop(sprintf("Refusing to overwrite existing files in `%s`. Re-run with `overwrite = TRUE` if replacement is intended.",
            app_dir), call. = FALSE)
    }
    bundle <- create_deexplorer_bundle(dge = dge, ..., msigdb_genesets = msigdb_genesets,
        msigdb_path = msigdb_path, fgsea_files = fgsea_files,
        sample_id_col = sample_id_col, gene_id_col = gene_id_col,
        gene_symbol_col = gene_symbol_col, fdr_cutoff = fdr_cutoff,
        prior_count = prior_count)
    saveRDS(bundle, file = file.path(app_dir, "app-data.rds"))
    write_app_wrapper_files(app_dir = app_dir, title = title)
    if (isTRUE(launch)) {
        shiny::runApp(app_dir)
    }
    invisible(list(app_dir = app_dir, ui = file.path(app_dir, 
        "ui.R"), server = file.path(app_dir, "server.R"), data = file.path(app_dir, 
        "app-data.rds")))
}
