\name{plotTaskSizes}
\alias{plotTaskSizes}
\name{plotTaskSizes}
\alias{plotTaskSizes}
\title{Plot distribution of gene counts per task (by system)}
\description{
Create a boxplot summarizing the number of genes associated with each task,
grouped by higher-level system (e.g. metabolic system). Useful to compare
task sizes across systems and samples.
}
\usage{
plotTaskSizes(taskInfo, OUTDIR)
}
\arguments{
  \item{taskInfo}{Data frame or tibble. Must contain at least the columns
    \code{Sample}, \code{Depth3}, \code{Depth1} and
    \code{GeneAssociatedToEssentialRxnsTask}. Rows represent task->gene
    associations.}
  \item{OUTDIR}{Character(1). Directory where the output PDF
    (\file{Task_size_system.pdf}) will be written. The directory must be
    writable.}
}
\details{
The function collapses task->gene associations to compute the number of genes
per task per sample, merges each task to its higher-level system using
\code{Depth1}, and renders a boxplot showing the distribution of task sizes
for each system. The plot is written as \file{Task_size_system.pdf} to
\code{OUTDIR}.
}
\value{
A data.frame with one row per task per sample summarising gene counts
(`\code{geneCount}`) is returned invisibly. The primary side-effect is the
creation of the PDF file in \code{OUTDIR}.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
## Example (assuming `taskInfo` is available):
## plotTaskSizes(taskInfo, OUTDIR = "./out")
}
\keyword{hplot}
