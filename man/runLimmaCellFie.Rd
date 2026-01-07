\name{runLimmaCellFie}
\alias{runLimmaCellFie}
\title{Run limma analysis on CellFie task scores and save summaries}
\description{
Prepare task-score matrices at a specified hierarchical depth from CellFie
output, run per-group limma analysis, generate heatmaps, and save a
workbook summarizing per-group top-table results.
}
\usage{
runLimmaCellFie(depth, cfs, meta, OUTDIRV, w = 9.5, height = 23)
}
\arguments{
  \item{depth}{Character. Name of the depth column in \code{cfs} to aggregate tasks by (e.g. "Depth3").}
  \item{cfs}{Data frame. CellFie output containing at minimum the columns \code{Sample}, \code{TaskScore}, and the hierarchical depth columns (e.g. \code{Depth1}, \code{Depth2}, \code{Depth3}).}
  \item{meta}{Data frame. Sample metadata that must include at least \code{Sample}, \code{Cancer_type}, \code{Condition}, \code{Batch}, \code{Paired} (values 'Paired'/'No') and \code{Patient}.}
  \item{OUTDIRV}{Character. Directory where output files (heatmaps and Excel summary) will be written. The directory will be created if it does not exist.}
  \item{w}{Numeric. Base width used to calculate heatmap widths (default 9.5).}
  \item{height}{Numeric. Height used for heatmap figures (default 23).}
}
\details{
The function reshapes \code{cfs} into a matrix of task scores per sample for
the requested depth, subsets samples by cancer type, decides whether a
paired design is appropriate (based on the \code{Paired} column and the
number of available pairs), and calls \code{limmaTest_CellFie} per cancer
type. It produces two heatmaps per cancer type (with and without sample
names) and writes a consolidated Excel workbook containing the top-table
results for all cancer types analyzed.
}
\value{
A list with components:
\item{cf5_Depth}{The reshaped task-by-sample matrix used for analysis.}
\item{rlst_list}{A named list of results returned by \code{limmaTest_CellFie} for each cancer type analyzed.}
}
\examples{
\dontrun{
# Example (illustrative only):
# res <- runLimmaCellFie("Depth3", cfs, meta, OUTDIRV = "./results/limma", w = 9.5, height = 23)
# str(res)
}
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\keyword{utilities}

