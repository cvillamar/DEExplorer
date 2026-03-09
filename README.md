# DEExplorer

DEExplorer is an R package component for building a bulk RNA-seq differential expression exploration app from:

- one `edgeR::DGEList` object
- one or more `limma::MArrayLM` fits (`eBayes()` and/or `treat()`)
- optional MSigDB genesets and precomputed fGSEA results

The main entrypoint is `deexplorer()`, which generates `ui.R`, `server.R`, and `app-data.rds` in an output folder.

## Basic Functionality

The generated Shiny app has two tabs:

1. `PCA`
- Interactive sample PCA (Plotly)
- PC x/y selection
- Metadata mapping to color, shape, and size
- Variance explained bar plot
- Sortable sample metadata table
- Downloads for raw counts, CPM, and log-CPM

2. `Explore DEG results`
- Contrast selector across all supplied fits
- Ranked DE table (`topTable()` or `topTreat()`, sorted by test statistic)
- Interactive MA and volcano plots with FDR coloring (`adj.P.Val < 0.05`)
- Hover-linked gene summary and log-CPM boxplot
- Manual gene highlighting (one gene per line; accepts gene symbols or Ensembl IDs)
- MSigDB geneset highlighting/filtering + precomputed fGSEA summary when available
- fGSEA enrichment barplot (top 10 pathways by |NES|, colored by FDR) with collection selector
- Gene annotation panel (Open Targets): protein/gene function description, tissue expression barplot, and disease associations table

## Installation

```r
# install.packages("remotes")
remotes::install_github("cvillamar/DEExplorer")
```

For SSH/private access:

```r
remotes::install_git("git@github.com:cvillamar/DEExplorer.git")
```

Then load the package:

```r
library(DEExplorer)
```

## Quick Start

```r
library(DEExplorer)

dge <- readRDS("example_input/rds/dge.Rds")
efit <- readRDS("example_input/rds/efit.Rds")
efit2 <- readRDS("example_input/rds/efit2.Rds")

fgsea_files <- list.files(
  path = "example_input",
  pattern = "fGSEA.*",
  full.names = TRUE
)

deexplorer(
  input_DGEList = dge,
  input_MArrayLM = list(efit = efit, efit2 = efit2),
  input_gsea_files = fgsea_files,
  out = "deexplorer_app",
  title = "This is the first deexplorer app",
  overwrite = TRUE
)
```

Run the generated app:

```r
shiny::runApp("deexplorer_app")
```

## Notes on Inputs

- `input_DGEList` must be an `edgeR::DGEList` with aligned `counts` and `samples`.
- `input_MArrayLM` should be a list of `MArrayLM` fits from the same `dge` gene universe.
- `input_gsea_files` is optional and can be:
  - a character vector of fGSEA TSV file paths (e.g. from `list.files(..., pattern = "fGSEA.*", full.names = TRUE)`)
  - a named list with `fgsea_files` (character vector) and/or `msigdb_genesets` (named list of gene sets)

## Generated Files

`deexplorer(..., out = "deexplorer_app")` writes:

- `deexplorer_app/ui.R`
- `deexplorer_app/server.R`
- `deexplorer_app/app-data.rds`
