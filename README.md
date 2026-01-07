# parseMetab

parseMetab: Metabolic activity inference and analysis for omics data

Utilities to parse proteomics / transcriptomics data, infer sample-wise metabolic
activities using CellFie or pathway-based GSVA, run limma differential testing
on inferred activities, and produce summary tables and visualizations tailored
for multi-cohort analyses (e.g. TCGA).

## Highlights

- Convert CellFie outputs and pathway mappings into tidy tables for downstream analysis.
- Run GSVA or CellFie-based scoring and perform limma-based differential testing.
- Plotting helpers: dotplots, volcano plots, heatmaps, circular/sample-size summaries, and more.
- Convenience wrappers to run multi-cohort analyses and produce publication-ready figures.

## Installation

This package is under development. To install locally from the package directory:

```r
# from R, run in the repo root (one level above the parseMetab folder)
install.packages("remotes")
remotes::install_github("https://github.com/Maj18/parseMetab", dependencies = TRUE)
```

## Quick start (usage sketch)

Below is a short example illustrating a typical workflow. Replace placeholder objects with your data.

```r
library(parseMetab)

# 1) Produce or load per-sample pathway/task scores (CellFie or GSVA)
#    Example: `cf` is a data.frame/tibble with columns Depth1, Depth2, Depth3, Sample, TaskScore

# 2) Run per-cohort differential testing (limma wrappers)
# rlst_list <- runLimmaCellFie(cf_list, meta = metadata, ...)

# 3) Create a dotplot summarizing frequently dysregulated pathways across cohorts
# makeLimmaDotplot_TCGA(rlst_list, adj.P.Val.cutoff = 0.05, OUTDIRV = "./plots/")

# 4) Other visual summaries
# doVolcano_CellFie(...)
# makeTaskScoreHeatmap_CellFie(...)
```

See the function help pages for argument details and examples: e.g. `?makeLimmaDotplot_TCGA`, `?runLimmaCellFie`.

## Development & documentation

While developing, keep docs and NAMESPACE in sync using roxygen2:

```r
# from R in the package root
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::document(pkg = "parseMetab")  # generate Rd and NAMESPACE from roxygen

# Run package checks
devtools::check(pkg = "parseMetab")
# or from a shell:
# R CMD check parseMetab
```

If `R CMD check` reports missing `Imports`, add them to `parseMetab/DESCRIPTION` under `Imports:`. If examples fail because they rely on heavy data or unavailable packages, use `\dontrun{}` in roxygen examples.

## Contributing

- Please open issues for bugs or feature requests.
- When editing documentation, prefer roxygen comments in `R/*.R` and then run `devtools::document()` to regenerate `man/` pages and `NAMESPACE`.

## License

This project uses the GPL-3 license (see `parseMetab/DESCRIPTION`).

## Contact

Author/Maintainer: YUAN LI <yuan.li@nbis.se>

---



