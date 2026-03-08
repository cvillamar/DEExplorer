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
- Lasso selection + Enrichr support (minimum 30 genes)
- Manual gene highlighting (one gene per line)
- MSigDB geneset highlighting/filtering + precomputed fGSEA summary when available

## Installation

```r
# install.packages("remotes")
remotes::install_github("user/repo")
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

deexplorer(
  in_dgelist = dge,
  in_marralm = list(efit = efit, efit2 = efit2),
  in_gsea = list(
    fgsea_dir = "example_input",
    msigdb_path = "example_input/msigdbr_hallmark_go_pathway_genesets_hs.rds"
  ),
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

- `in_dgelist` must be an `edgeR::DGEList` with aligned `counts` and `samples`.
- `in_marralm` should be a list of `MArrayLM` fits from the same `dge` gene universe.
- `in_gsea` is optional and can be:
  - a folder with `fGSEA_*.tsv` files
  - an `.rds` geneset file
  - a named geneset list/data.frame
  - a config list containing `fgsea_dir`, `msigdb_path`, and/or `msigdb_genesets`

## Generated Files

`deexplorer(..., out = "deexplorer_app")` writes:

- `deexplorer_app/ui.R`
- `deexplorer_app/server.R`
- `deexplorer_app/app-data.rds`
