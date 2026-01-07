\name{getCirPackingPlot_CellFie}
\alias{getCirPackingPlot_CellFie}
\name{getCirPackingPlot_CellFie}
\alias{getCirPackingPlot_CellFie}
\title{Circular packing plot of CellFie task hierarchy}
\description{
Create a circular packing plot that summarizes CellFie task hierarchy and
their grouping by higher-level systems. Labels are shortened for display and
the plot is written as a PDF to the supplied output directory.
}
\usage{
getCirPackingPlot_CellFie(cfs, out_dir)
}
\arguments{
    \item{cfs}{Data frame or tibble. CellFie output (as returned by
        \code{processCellFieOutput()}) containing at minimum the hierarchical depth
        columns (e.g. \code{Depth1}, \code{Depth2}, \code{Depth3}) and task
        identifiers; rows represent task membership or task->gene mappings.}
    \item{out_dir}{Character(1). Directory where the output PDF
        (\file{TaskSummary.pdf}) will be written. The directory must be writable.}
}
\details{
The function assembles a three-level hierarchy from \code{Depth1/Depth2/Depth3}
and plots it as a circle packing diagram using \pkg{ggraph}. Shortened
labels are used for readability. The plot file is saved as
\file{TaskSummary.pdf} in \code{out_dir}.
}
\value{
Invisibly returns the \pkg{ggplot2} object used to draw the plot. The primary
side-effect is creation of the PDF file in \code{out_dir}.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
## Example (assuming a CellFie summary `cfs` is available):
## getCirPackingPlot_CellFie(cfs, out_dir = "./out")
}
\keyword{hplot}
