\name{ActivityBoxplot_byCancer}
\alias{ActivityBoxplot_byCancer}
\title{Boxplot of mean task scores by cancer type}
\description{
Generate a boxplot showing mean CellFie task scores (aggregated at
the task level) across cancer types. The plot is annotated with total
sample counts per cancer and saved as a PDF.
}
\usage{
ActivityBoxplot_byCancer(meta, cf, OUT_DIR, w = 2.75, h = 3, title = "Pan-cancer tasks")
}
\arguments{
  \item{meta}{A data.frame with at least the columns \code{Sample}, \code{Condition}, and \code{Cancer_type}.}
  \item{cf}{A data.frame with CellFie output containing at minimum \code{Sample}, \code{TaskScore}, and hierarchical depth columns such as \code{Depth1}, \code{Depth2}, and \code{Depth3}.}
  \item{OUT_DIR}{Character. Directory where the output PDF \code{meanTaskScores_cancers.pdf} will be written.}
  \item{w}{Numeric. Figure width in inches (default 2.75).}
  \item{h}{Numeric. Figure height in inches (default 3).}
  \item{title}{Character. Plot title (default \code{"Pan-cancer tasks"}).}
}
\details{
The function joins \code{cf} with \code{meta}, computes the mean task
score for each task and cancer type (using \code{Depth3}), orders
cancer types by their average activity, and renders a boxplot. The
resulting plot is written to \code{OUT_DIR}/meanTaskScores_cancers.pdf.
}
\value{
Invisibly returns the ggplot object after writing the PDF to
\code{OUT_DIR}. The side-effect (saved PDF) is the primary output.
}
\examples{
\dontrun{
# ActivityBoxplot_byCancer(meta, cf, OUT_DIR = "./plots/", w = 3, h = 4, title = "My Study")
}
}
\author{Yuan Li {yuan.li@nbis.se}}
\keyword{visualization}

