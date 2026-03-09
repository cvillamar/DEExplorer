# DEExplorer

DEExplorer is an R package for building interactive Shiny apps to explore bulk RNA-seq differential expression results from:

- one `edgeR::DGEList` object
- one or more `limma::MArrayLM` fits (`eBayes()` and/or `treat()`)
- optional MSigDB gene sets (`.rds`) and precomputed fGSEA results (`.tsv`)

The main entrypoint is `deexplorer()`, which generates `ui.R`, `server.R`, and `app-data.rds` in an output folder.

## Installation

```r
# install.packages("remotes")
remotes::install_github("cvillamar/DEExplorer")
```

For SSH/private access:

```r
remotes::install_git("git@github.com:cvillamar/DEExplorer.git")
```

## Quick Start

```r
library(DEExplorer)

dge  <- readRDS("example_input/rds/dge.Rds")
efit <- readRDS("example_input/rds/efit.Rds")
efit2 <- readRDS("example_input/rds/efit2.Rds")

fgsea_files <- list.files(
  path = "example_input",
  pattern = "fGSEA.*",
  full.names = TRUE
)

deexplorer(
  input_DGEList  = dge,
  input_MArrayLM = list(efit = efit, efit2 = efit2),
  input_gsea_files = fgsea_files,
  msigdb_path = "example_input/msigdbr_hallmark_go_pathway_genesets_hs.rds",
  out = "deexplorer_app",
  title = "This is the first deexplorer app",
  overwrite = TRUE
)

shiny::runApp("deexplorer_app")
```

## App Features

The generated Shiny app has two tabs:

### Tab 1: PCA

- Interactive sample PCA scatter plot (Plotly)
- PC x/y axis selection
- Metadata mapping to color, shape, and size
- Variance explained bar plot
- Sortable sample metadata table
- Downloads for raw counts, CPM, and log-CPM matrices

### Tab 2: Explore DEG Results

**Contrast selection and DE table**
- Contrast selector across all supplied fits (eBayes and treat)
- Ranked DE table (`topTable()` / `topTreat()`, sorted by test statistic)
- Interactive MA and volcano plots with FDR coloring (`adj.P.Val < 0.05`)
- Hover-linked gene summary and log-CPM boxplot grouped by sample metadata

**Gene highlighting**
- *Highlight genes* text box: paste gene symbols or Ensembl IDs (one per line) to mark them with distinct markers on MA and volcano plots
- *MSigDB gene set* autocomplete: start typing any word (e.g. "interferon", "apoptosis") to search loaded gene sets. Selecting a gene set filters the DE table to its member genes and highlights them on plots with orange ring markers. If precomputed fGSEA results exist for that gene set and contrast, the leading edge genes are additionally marked with diamond markers.

**fGSEA enrichment**
- Barplot of top 10 enriched pathways by |NES|, colored by FDR significance
- Collection selector dropdown (Hallmarks, GO:BP, KEGG, etc.)
- Per-gene-set fGSEA summary (NES, p-value, leading edge) shown in the sidebar

**Gene annotation (Open Targets)**
- Click any gene in the DE table, MA plot, or volcano plot to load external annotations
- Protein/gene function description (UniProt via Open Targets API)
- Tissue expression barplot (GTEx RNA TPM, aggregated by organ)
- Disease associations table with per-datasource scores and column tooltips

## Parameters

| Parameter | Description |
|---|---|
| `input_DGEList` | `edgeR::DGEList` with aligned `counts` and `samples` |
| `input_MArrayLM` | One `MArrayLM` or a named list of fits from the same gene universe |
| `input_gsea_files` | Character vector of fGSEA `.tsv` file paths (e.g. from `list.files(..., pattern = "fGSEA.*", full.names = TRUE)`) |
| `msigdb_path` | Path to an `.rds` file containing a named list of gene sets (names = gene-set names, values = character vectors of gene symbols). Populates the MSigDB autocomplete selector. |
| `out` | Output directory for generated app files |
| `title` | App title |
| `overwrite` | Overwrite existing files (default `TRUE`) |
| `sample_id_col` | Column name in `dge$samples` for sample IDs (default `"Sample"`) |
| `gene_id_col` | Column name in `dge$genes` for gene IDs (default `"ENSEMBL"`) |
| `gene_symbol_col` | Column name in `dge$genes` for gene symbols (default `"SYMBOL"`) |
| `fdr_cutoff` | FDR threshold for significance coloring (default `0.05`) |
| `prior_count` | Prior count for log-CPM transformation (default `2`) |

## Generated Files

`deexplorer(..., out = "deexplorer_app")` writes:

- `deexplorer_app/ui.R`
- `deexplorer_app/server.R`
- `deexplorer_app/app-data.rds`

Run with `shiny::runApp("deexplorer_app")`.
