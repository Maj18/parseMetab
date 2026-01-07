\name{doVolcano_CellFie}
\alias{doVolcano_CellFie}
\title{Volcano plot for limma differential results}
\description{
Create a volcano plot showing log2 fold-change vs -log10(adjusted
p-value) for limma differential testing results. Significant features
are colored by direction and the largest fold-changes are labeled.
}
\usage{
doVolcano_CellFie(resi, NAME, adj.P.Val.cutoff = 0.05, logFC.cutoff = 0.3, lfcShrink = F)
}
\arguments{
  \item{resi}{A data.frame-like limma top-table containing at least the columns \code{Task}, \code{logFC}, \code{P.Value}, \code{adj.P.Val}, and \code{Depth}.}
  \item{NAME}{Character. Short name used for the plot title (and useful when saving or annotating figures).}
  \item{adj.P.Val.cutoff}{Numeric. Threshold on adjusted p-value used to call significance (default 0.05).}
  \item{logFC.cutoff}{Numeric. Log2 fold-change threshold used to classify up/down regulation (default 0.3).}
  \item{lfcShrink}{Logical. If \code{TRUE} the plot title indicates that log-fold-change shrinkage was applied.}
}
\details{
Rows with missing adjusted p-values are removed before plotting. Features
are classified into \dQuote{up}, \dQuote{down}, or \dQuote{nonSig}
based on \code{adj.P.Val.cutoff} and \code{logFC.cutoff}. The function
scales point sizes according to the \code{Depth} column (e.g. Depth3 =
smallest). Labels are placed on extreme fold-changes (top ~2% by
absolute \code{logFC}) using \pkg{ggrepel} to reduce overlap.
}
\value{
A \code{ggplot} object is returned. The caller can print or save the
plot as needed.
}
\examples{
\dontrun{
# p <- doVolcano_CellFie(resi, NAME = "MyContrast")
# print(p)
}
}
\author{Yuan Li {yuan.li@nbis.se}}
\keyword{visualization}
%% ~~objects to See Also as \code{\link{help}}, ~~~
