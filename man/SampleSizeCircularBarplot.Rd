\name{SampleSizeCircularBarplot}
\alias{SampleSizeCircularBarplot}
\title{Circular barplot summarizing sample counts per cancer type}
\description{
Create a circular stacked barplot showing the number of samples per
cancer type, split by condition (for example Normal vs Tumor). Each
cancer label is annotated with the total sample count to aid
interpretation.
}
\usage{
SampleSizeCircularBarplot(meta, title = "Pan-cancer tasks", OUT_DIR, h = 4, w = 4.5)
}
\arguments{
  \item{meta}{A data.frame with at least the columns \code{Sample}, \code{Condition} and \code{Cancer_type}.}
  \item{title}{Character. Plot title (default \code{"Pan-cancer tasks"}).}
  \item{OUT_DIR}{Character. Path or directory prefix where the resulting PDF \code{SampleSizeCircularBarplot.pdf} will be written.}
  \item{h}{Numeric. Figure height in inches (default 4).}
  \item{w}{Numeric. Figure width in inches (default 4.5).}
}
\details{
The function aggregates sample counts by cancer type and condition,
reformats the table for plotting, and renders a stacked barplot in
polar coordinates. Cancer labels include the total sample count between
parentheses for easy reference.
}
\value{
Invisibly returns the ggplot object after writing a PDF to \code{OUT_DIR}.
The primary side-effect is the saved PDF file.
}
\examples{
\dontrun{
# SampleSizeCircularBarplot(meta, title = "MyStudy", OUT_DIR = "./plots/", h = 4, w = 4.5)
}
}
\author{Yuan Li {yuan.li@nbis.se}}
\keyword{visualization}
