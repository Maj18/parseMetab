\name{processCellFieOutput}
\alias{processCellFieOutput}
\title{Process CellFie output and produce summaries and figures}
\description{
Import per-group CellFie output, reshape and summarise task-level scores,
produce aggregated tables at multiple hierarchical depths, and render
diagnostic plots (boxplots and a circle-packing diagram). Results and plots
are written to a Summary subdirectory under the provided `outdir`.
}
\usage{
processCellFieOutput(combined = FALSE, outdir = paste0(OUTDIR, "/CellFieOut/"), SampleNames, meta, samples2keep)
}
\arguments{
  \item{combined}{Logical(1). If TRUE, treat all samples as coming from a single combined run.}
  \item{outdir}{Character(1). Path to the directory that contains per-group CellFie output folders.}
  \item{SampleNames}{Character vector. Vector of sample IDs in the same order used when running CellFie.}
  \item{meta}{Data frame. Sample metadata containing at least the columns \code{Sample} and \code{Cancer_type}; additional columns used include \code{Condition} and \code{Cancer_Condition}.}
  \item{samples2keep}{Character vector. Sample IDs to keep after import; used to filter the merged results.}
}
\details{
The function performs the following steps: (1) reads per-group CellFie CSVs;
(2) cleans and reshapes task-level results so tasks are rows and samples
are columns; (3) merges sample-level metadata; (4) computes aggregated
task-score summaries at Depth1/Depth2/Depth3; and (5) writes Excel tables,
RDS objects and PDF figures into \file{<outdir>/../Summary/}.
}
\value{
Returns invisibly a single data frame containing the concatenated, post-processed CellFie results across groups; the function's primary effect is writing summary files and figures to disk.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
## processCellFieOutput(combined = FALSE, outdir = "./CellFieOut/", SampleNames = samples, meta = meta, samples2keep = samples)
}
\keyword{utilities}

