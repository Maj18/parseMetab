\name{diffEffectBoxplot_bySystem}
\alias{diffEffectBoxplot_bySystem}
\title{Boxplot of differential effect logFC (weighted by gene expression level) by system}
\description{
Create a boxplot summarizing the differential effect logFC (weighted by gene expression level) for tasks
grouped at the system level (Depth1). The function maps task-level
statistics (Depth3) to systems using the provided task hierarchy in
\code{cf} and writes a PDF plot to the specified output directory.
}
\usage{
diffEffectBoxplot_bySystem(rlst_list, cf, OUT_DIR, w = 2, h = 3.5, title = "Pan-cancer tasks")
}
\arguments{
  \item{rlst_list}{A named list. Each element should contain a
    data.frame \code{rslt_of_interest} with at least the limma columns
    \code{t}, \code{adj.P.Val}, \code{AveExpr} and \code{logFC}. Typical input is the
    output from \code{runLimmaCellFie} or \code{limmaTest_CellFie}.}
  \item{cf}{A data.frame. CellFie output containing at least
    \code{Depth1} and \code{Depth3} used to map tasks to systems.}
  \item{OUT_DIR}{Character. Directory where the PDF file
    (\code{diffEffectBoxplot_system.pdf}) will be written.}
  \item{w}{Numeric. Width of the output PDF in inches (default 2).}
  \item{h}{Numeric. Height of the output PDF in inches (default 3.5).}
  \item{title}{Character. Plot title (default: "Pan-cancer tasks").}
}
\details{
The function concatenates limma results from \code{rlst_list}, annotates
each row with its originating cohort, merges task-level rows with the
task hierarchy in \code{cf}, and computes median t-statistics per
system for ordering. The resulting boxplot shows the distribution of
t-statistics for each system across all tasks.
}
\value{Invisibly returns the ggplot object used to render the figure; the main
side-effect is the written PDF file.}
\author{Yuan Li {yuan.li@nbis.se}}
\examples{
\dontrun{
## Example usage (assuming 'rlst_list' and 'cf' are prepared):
diffEffectBoxplot_bySystem(rlst_list, cf, OUT_DIR = "./plots/", w = 3, h = 4,
                          title = "Study Systems")
}
}
\keyword{hplot}
