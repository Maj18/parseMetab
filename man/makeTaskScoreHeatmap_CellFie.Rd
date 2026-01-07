\name{makeTaskScoreHeatmap_CellFie}
\alias{makeTaskScoreHeatmap_CellFie}
\title{Generate a task-score heatmap across all samples and groups}
\description{
Create a combined heatmap of tasks that are frequently significant
across cancer types using per-group limma results. The function selects
significant tasks (up or down) from \code{rlst_list}, orders them by
frequency, and plots them across all samples in \code{cf5_Depth} with
column annotations derived from \code{meta}.
}
\usage{
makeTaskScoreHeatmap_CellFie(rlst_list, cf5_Depth, OUTDIRV, w = 10, h = 23, adj.P.Val.cutoff = 0.05, meta)
}
\arguments{
  \item{rlst_list}{Named list. Output from \code{runLimmaCellFie} or \code{limmaTest_CellFie}; each element must contain a \code{rslt_of_interest} data.frame (topTable).}
  \item{cf5_Depth}{Matrix or data.frame. Task-by-sample matrix (rows = tasks at the chosen depth, columns = samples).}
  \item{OUTDIRV}{Character. Directory where generated heatmap files will be saved.}
  \item{w}{Numeric. Base width used to compute heatmap figure widths (default 10).}
  \item{h}{Numeric. Height used for the heatmap output (default 23).}
  \item{adj.P.Val.cutoff}{Numeric. Adjusted p-value threshold used to select significant tasks (default 0.05).}
  \item{meta}{Data frame. Sample metadata containing at least columns \code{Sample}, \code{Cancer_type}, and \code{Condition}, used to construct column annotations.}
}
\details{
The function collects significant tasks across all cancer types in
\code{rlst_list} (splitting into up- and down-regulated sets), ranks
them by how often they occur as significant, and then subsets
\code{cf5_Depth} to those tasks for plotting. Two heatmaps are
produced: one with sample names and one without. The column annotation
frame contains cancer type and group (Condition).
}
\value{
Invisibly returns the plotting object (if available) or \code{NULL}. The
primary effect is side-effect: PDF heatmap files are written to
\code{OUTDIRV}.
}
\examples{
\dontrun{
# Example (illustrative only):
# makeTaskScoreHeatmap_CellFie(rlst_list, cf5_Depth, OUTDIRV = "./results", w = 10, h = 20, adj.P.Val.cutoff = 0.05, meta)
}
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\keyword{visualization}
%%  \item{comp2 }{Description of 'comp2'}
