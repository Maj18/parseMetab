\name{makeLimmaDotplot_CellFie}
\alias{makeLimmaDotplot_CellFie}
\title{Dot-plot summary of limma results across cancer types}
\description{
Summarize per-cancer limma results in a dot-plot where each point
represents a task (or subsystem/system) in a given cancer type. Point
shape indicates direction (Up/Down/NonSig), color maps the test
statistic (t), and point size reflects adjusted p-value.
}
\usage{
makeLimmaDotplot_CellFie(rlst_list, adj.P.Val.cutoff, w = NULL, h = NULL, OUTDIRV, Depth = "Depth3", lowfigure = F, title = "Pan-cancer tasks", suffix = "Hu2025", taskHierarchy = cfs, nested = T, filter = FALSE, probs_cutoff = 0.35)
}
\arguments{
  \item{rlst_list}{Named list. Output returned by \code{runLimmaCellFie} or entries produced by \code{limmaTest_CellFie}; each element must contain a \code{rslt_of_interest} data.frame (topTable).}
  \item{adj.P.Val.cutoff}{Numeric. Adjusted p-value cutoff used to select significant tasks for plotting.}
  \item{w}{Numeric or \code{NULL}. Plot width in inches; sensible default width is computed when \code{NULL}.}
  \item{h}{Numeric or \code{NULL}. Plot height in inches; sensible default height is computed when \code{NULL}.}
  \item{OUTDIRV}{Character. Directory where the generated PDF dotplot will be saved.}
  \item{Depth}{Character. Which hierarchical depth column from \code{taskHierarchy} to use as the task label (e.g. \code{"Depth3"}).}
  \item{lowfigure}{Logical. If \code{TRUE} apply compact theming suitable for small/embedded figures.}
  \item{title}{Character. Plot title.}
  \item{suffix}{Character. Filename suffix appended to the saved PDF.}
  \item{taskHierarchy}{Data.frame. A mapping table that contains columns \code{Depth1}, \code{Depth2}, \code{Depth3}, and a task identifier column used to attach hierarchy labels.}
  \item{nested}{Logical. If \code{TRUE}, facets are nested by \code{Depth1} (requires \code{ggh4x}).}
  \item{filter}{Logical. If \code{TRUE}, restrict tasks to the most frequent ones according to \code{probs_cutoff}.}
  \item{probs_cutoff}{Numeric between 0 and 1. Quantile threshold used when filtering tasks (default 0.35).}
}
\details{
The function collapses per-cancer topTable results into a long table,
joins hierarchical labels from \code{taskHierarchy}, and selects tasks
that meet the adjusted p-value threshold. The dot-plot encodes
direction, significance, and test statistic; optional filtering keeps
only the most recurrently significant tasks across cancers.
}
\value{
A \code{ggplot} object is returned (invisibly) and a PDF file is written
to \code{OUTDIRV}. The returned object can be printed or saved by the
user.
}
\examples{
\dontrun{
# p <- makeLimmaDotplot_CellFie(rlst_list, adj.P.Val.cutoff = 0.05, OUTDIRV = "./results", Depth = "Depth3", taskHierarchy = cfs)
# print(p)
}
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\keyword{visualization}
