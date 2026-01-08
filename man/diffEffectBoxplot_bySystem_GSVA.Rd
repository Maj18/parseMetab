\name{diffEffectBoxplot_bySystem_GSVA}
\alias{diffEffectBoxplot_bySystem_GSVA}
\title{Differential effect boxplot by metabolic class}
\description{
Generate a boxplot summarizing differential effect logFC (weighted by gene expression level) for
KEGG metabolic pathways aggregated by metabolic class. The function collects
limma test results across input groups, maps pathways to classes via the
provided hierarchy table, and writes a PDF boxplot to \code{OUT_DIR}.
}
\usage{
diffEffectBoxplot_bySystem_GSVA(rlst_list, taskHierarchy, OUT_DIR, w = 2, h = 3.5, title = "Zhou2020")
}
\arguments{
  \item{rlst_list}{A named list of results returned by \code{GSVAlimmaTest()} (one
    element per group). Each element must contain a component \code{rslt_of_interest}
    (a data frame produced by \code{limma::topTable}) with at least the columns
    \code{t}, \code{logFC}, \code{AveExpr} and \code{adj.P.Val}.}
  \item{taskHierarchy}{A data.frame or tibble with at least two columns:
    \code{gs_name} (pathway name or ID) and \code{class} (the higher-level
    metabolic class). This table maps pathways to classes for grouping and plotting.}
  \item{OUT_DIR}{Character scalar. Directory where the output PDF
    (\file{diffEffectBoxplot_system.pdf}) will be written.}
  \item{w}{Numeric. Output figure width (in inches).}
  \item{h}{Numeric. Output figure height (in inches).}
  \item{title}{Character. Plot title.}
}
\details{
The function aggregates limma t-statistics across input results, assigns each
pathway to a metabolic class using \code{taskHierarchy}, computes the median
effect per class to determine plotting order, and renders a boxplot of the
effects per class. The plot is saved as \file{diffEffectBoxplot_system.pdf}
in \code{OUT_DIR}.
}
\value{
Invisibly returns the path to the saved PDF file; the primary side-effect is
that \file{diffEffectBoxplot_system.pdf} is created in \code{OUT_DIR}.
}
\author{
  Yuan Li /email{yuan.li@nbis.se}
}
\examples{
## Example (assuming results and hierarchy are available):
## diffEffectBoxplot_bySystem_GSVA(rlst_list, taskHierarchy, OUT_DIR = "/path/to/out")
}
% For non-standard keywords, use \concept instead of \keyword:
