\name{ActivityBoxplot_bySystem}
\alias{ActivityBoxplot_bySystem}
\title{Boxplot of mean CellFie task scores by system (Depth1)}
\description{
Create a boxplot that summarizes mean CellFie task scores at the
system level (Depth1). The plot shows the distribution of mean task
scores (aggregated over Depth3 tasks) for each system and is written
to a PDF file in the specified output directory.
}
\usage{
ActivityBoxplot_bySystem(cf, OUT_DIR, w = 2, h = 3.5, title = "Pan-cancer tasks")
}
\arguments{
  \item{cf}{A data.frame. CellFie output with at least the columns
    \code{Sample}, \code{TaskScore}, \code{Depth1} and \code{Depth3}.}
  \item{OUT_DIR}{Character. Directory where the PDF plot will be saved.}
  \item{w}{Numeric. Width of the output PDF in inches (default 2).}
  \item{h}{Numeric. Height of the output PDF in inches (default 3.5).}
  \item{title}{Character. Title to display on the plot (default: "Pan-cancer tasks").}
}
\details{
The function computes the mean task score for each combination of
\code{Depth1} (system) and \code{Depth3} (task), then orders systems by
their median/mean activity before plotting. The resulting PDF file is
named \code{meanTaskScores_system.pdf} and saved inside \code{OUT_DIR}.
}
\value{Invisibly returns the ggplot object used to render the figure. The main
side effect is the creation of the PDF file.}
\examples{
\dontrun{
## Example (assuming 'cf' is a CellFie-style data.frame):
ActivityBoxplot_bySystem(cf, OUT_DIR = "./plots/", w = 2, h = 3.5,
                         title = "Study Systems")
}
}
/author{Yuan Li {yuan.li@nbis.se}}
\keyword{hplot}
