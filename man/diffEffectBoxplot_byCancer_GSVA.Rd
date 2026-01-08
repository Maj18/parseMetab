\name{diffEffectBoxplot_byCancer_GSVA}
\alias{diffEffectBoxplot_byCancer_GSVA}
\title{Differential effect boxplot by group (e.g. cancer type)}
\description{
Generate a boxplot summarizing differential effect logFC (weighted by gene expression level) for
KEGG metabolic pathways across groups (for example, cancer types). The
function collects limma test results from each supplied group, computes a
group-level ordering by median effect, and writes a PDF boxplot to the
specified output directory.
}
\usage{
diffEffectBoxplot_byCancer_GSVA(rlst_list, OUT_DIR, w = 2.75, h = 3, title = "Zhou2020")
}
\arguments{
  \item{rlst_list}{Named list. Results as returned by \code{GSVAlimmaTest()}
    (one element per group). Each element must contain a component
    \code{rslt_of_interest}, the data frame produced by \code{limma::topTable},
    and that data frame must include at minimum the columns \code{t}, \code{AveExpr},
    \code{logFC} and \code{adj.P.Val}.}
  \item{OUT_DIR}{Character scalar. Directory where the output PDF
    (\file{diffEffectBoxplot_cancer.pdf}) will be written. The directory
    must be writable by the caller.}
  \item{w}{Numeric. Figure width in inches (default: 2.75).}
  \item{h}{Numeric. Figure height in inches (default: 3.0).}
  \item{title}{Character. Plot title (default: "Zhou2020").}
}
\details{
The function binds per-group limma results into a single table, assigns each
pathway to its group, computes the median t-statistic per group to determine
the plotting order, and renders a boxplot of t-statistics for each group.
The resulting plot is saved as \file{diffEffectBoxplot_cancer.pdf} in
\code{OUT_DIR}.
}
\value{
Invisibly returns the path to the saved PDF file. The primary side-effect
is creation of \file{diffEffectBoxplot_cancer.pdf} in \code{OUT_DIR}.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
## Example (assuming results are available in `rlst_list`):
## diffEffectBoxplot_byCancer_GSVA(rlst_list, OUT_DIR = "/path/to/out")
}
% Keywords
\keyword{hplot}
