\name{getCirPackingPlot_GSVA}
\alias{getCirPackingPlot_GSVA}
\title{Circular packing plot of KEGG metabolic pathway hierarchy}
\description{
Create a circular packing plot that summarizes KEGG metabolic pathways and
their classification. Nodes represent pathway classes and individual pathways;
node size indicates the number of member pathways.
}
\usage{
getCirPackingPlot_GSVA(kegg_metab_db_table, outdir)
}
\arguments{
    \item{kegg_metab_db_table}{A data frame or tibble with at least three columns:
        \code{gs_name} (pathway name or ID), \code{gene_symbol} (member gene symbol),
        and \code{class} (pathway class).}
    \item{outdir}{A character scalar: directory where the plot PDF will be written.
        The output file will be saved as \file{KEGG_metabolism/TaskSummary.pdf} inside
        \code{outdir}.}
}
\details{
The function builds a hierarchical representation (Metabolism -> Class -> Pathway),
creates a graph suitable for circle-packing, and renders a labeled plot. Labels are
shortened to improve readability. The plot is saved as a PDF to
\file{KEGG_metabolism/TaskSummary.pdf} under the supplied \code{outdir}.
}
\value{
Invisibly returns the generated \code{ggplot2} object. The primary side-effect is
that a PDF file is written to the specified \code{outdir}.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
## Example (assuming the KEGG table and outdir exist):
## getCirPackingPlot_GSVA(kegg_metab_db_table, outdir = "/path/to/output")
}

