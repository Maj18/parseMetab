\name{visualizeSigNr}
\alias{visualizeSigNr}
\title{Number of significant pathways per group}
\description{
Generate a bar chart summarizing the number of significantly different
pathways per group (up- and down-regulated) using limma results produced by
\code{GSVAlimmaTest}.
}
\usage{
visualizeSigNr(OUTDIR, GSVA_limma_rslt)
}
\arguments{
  \item{OUTDIR}{Character scalar. Directory where the output PDF
    (\file{Sig_pathway_nr.pdf}) will be written.}
  \item{GSVA_limma_rslt}{A named list of results returned by \code{GSVAlimmaTest()}
    (one element per group). Each element must contain a component
    \code{rslt_of_interest} (a data frame, as produced by \code{limma::topTable})
    with at least the columns \code{logFC} and \code{adj.P.Val}.}
}
\details{
The function counts pathways with adjusted p-value < 0.05 per group, splits
them by direction (up/down) and renders a stacked bar chart. The plot is
written to \file{Sig_pathway_nr.pdf} in \code{OUTDIR}.
}
\value{
Invisibly returns the path to the saved PDF file; the primary side-effect is
creating the PDF in \code{OUTDIR}.
}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\examples{
## Example (assuming results are available):
## visualizeSigNr(OUTDIR = "/path/to/out", GSVA_limma_rslt = my_results)
}

