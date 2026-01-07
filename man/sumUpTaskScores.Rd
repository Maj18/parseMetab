\name{sumUpTaskScores}
\alias{sumUpTaskScores}
\name{sumUpTaskScores}
\alias{sumUpTaskScores}
\title{Aggregate CellFie task scores to a higher-level category}
\description{
Summarize CellFie task-level scores at a specified hierarchical depth (for
example \code{Depth1} or \code{Depth2}) and produce a per-sample summary
table plus a PDF visualization showing the fraction of active tasks per sample.
}
\usage{
sumUpTaskScores(cfs, Depth, h = 12, w = 30, OUTDIR)
}
\arguments{
  \item{cfs}{Data frame or tibble. CellFie output (as produced by
    \code{processCellFieOutput()}) containing at least the columns
    \code{Sample}, the chosen \code{Depth} (e.g. \code{Depth1}),
    \code{TaskScore}, \code{BinaryTaskScore} and \code{ActiveTaskScore}.}
  \item{Depth}{Character(1). Name of the hierarchical depth to aggregate
    (e.g. \code{"Depth1"} or \code{"Depth2"}).}
  \item{h}{Numeric. Height of the output PDF in inches (default: 12).}
  \item{w}{Numeric. Width of the output PDF in inches (default: 30).}
  \item{OUTDIR}{Character(1). Directory where the output PDF
    (\file{BinaryTaskScores_<Depth>.pdf}) will be written. The directory must
    be writable.}
}
\details{
The function groups task-level scores by the supplied \code{Depth} and
\code{Sample}, computes summary statistics (means, sums, fraction and
count of active tasks), and writes a PDF named
\file{BinaryTaskScores_<Depth>.pdf} to \code{OUTDIR} showing per-sample
activity. The returned data frame contains the computed summary columns.
}
\value{
A data.frame with one row per sample/category, containing columns:
\item{MeanTaskScore}{Mean of raw task scores for the group and sample.}
\item{MeanActiveTaskScore}{Mean of active task scores (TaskScore * Binary flag).}
\item{SumTaskScore}{Sum of raw task scores.}
\item{SumActiveTaskScore}{Sum of active task scores.}
\item{FractionOfActiveTasks}{Fraction of tasks marked active for the sample.}
\item{NrOfActiveTasks}{Number of active tasks for the sample.}
Invisibly returns the same data.frame; the primary side-effect is the PDF
written to \code{OUTDIR}.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
## Example (assuming `cfs` is available):
## sumUpTaskScores(cfs, Depth = "Depth1", h = 12, w = 30, OUTDIR = "./out")
}
\keyword{hplot}
