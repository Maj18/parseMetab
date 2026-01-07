\name{parseMetab-package}
\alias{parseMetab-package}
\alias{parseMetab}
\docType{package}
\title{parseMetab: Metabolic activity inference and analysis for omics data}
\description{
Utilities to parse proteomics/transcriptomics data, infer sample-wise
metabolic activities using CellFie or pathway-based GSVA, run limma
differential testing on inferred activities, and produce summary tables
and visualizations tailored for multi-cohort (e.g. TCGA) analyses.
}
\details{
This package provides a collection of convenience wrappers and plotting
helpers commonly used when performing pathway- or task-level metabolic
activity analysis from high-throughput omics data. Key capabilities
include:

\itemize{
	\item Importing and reshaping CellFie outputs and pathway mappings.
	\item Running GSVA for KEGG or curated gene-set collections and
		summarizing results.
	\item Running limma (optionally paired) across cohorts and extracting
		top-tables per cohort (helpers: \code{GSVAlimmaTest},
		\code{limmaTest_CellFie}, \code{runLimmaCellFie}).
	\item Compact visual summaries: dotplots of frequently dysregulated
		pathways (\code{makeLimmaDotplot_TCGA}, \code{makeLimmaDotplot_CellFie}),
		volcano plots (\code{doVolcano_CellFie}), heatmaps
		(\code{makeTaskScoreHeatmap_CellFie}), and circular/sample-size
		summaries (\code{SampleSizeCircularBarplot}).
	\item Additional helpers for aggregating and plotting task scores,
		including \code{sumUpTaskScores}, \code{plotTaskSizes},
		\code{SigNrBarplot}, and the \code{diffEffectBoxplot_*} series.
}

The package is intended to speed up exploratory analyses and figure
generation for multi-cohort metabolic activity studies. Many functions
assume tidy input (tibbles/data.frames) and integration with the
\pkg{dplyr} / \pkg{ggplot2} ecosystem.
}
\author{YUAN LI <yuan.li@nbis.se>}
\references{
If you use this package in published work, please cite the methods and
databases used (e.g. limma, GSVA, KEGG). Suggested citations:
\itemize{
	\item Smyth, G. K. (2004). Linear models and empirical bayes methods
		for assessing differential expression in microarray experiments.
		Statistical Applications in Genetics and Molecular Biology.
	\item Hanzelmann, S., Castelo, R. and Guinney, J. (2013). GSVA: gene set
		variation analysis for microarray and RNA-seq data. BMC Bioinformatics.
}
}
\keyword{package}

\section{Usage}{
Load the package and call the helpers. Most plotting functions accept a
named list of per-cohort limma results (as returned by
\code{runLimmaCellFie} / \code{GSVAlimmaTest}) and an output directory
where figures will be saved.
}
\seealso{\code{\link{GSVAlimmaTest}}, \code{\link{makeLimmaDotplot_TCGA}}, \code{\link{makeTaskScoreHeatmap_CellFie}}}
\examples{
\dontrun{
library(parseMetab)

# Typical workflow outline (pseudocode):
# 1. Prepare expression matrix and metadata
# 2. Run GSVA or CellFie to obtain per-sample pathway/task scores
# 3. Run limma across cohorts: rlst_list <- runLimmaCellFie(...)
# 4. Summarize and plot: makeLimmaDotplot_TCGA(rlst_list, adj.P.Val.cutoff = 0.05, OUTDIRV = "./plots")
}
}
