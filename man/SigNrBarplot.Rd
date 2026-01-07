\name{SigNrBarplot}
\alias{SigNrBarplot}
\title{Number of significantly changed tasks per cohort}
\description{
Count and visualize the number of tasks that are significantly
up- or down-regulated per cancer cohort based on limma results and
plot the counts as a bar chart. Up-regulated counts are plotted as
positive bars and down-regulated counts as negative bars.
}
\usage{
SigNrBarplot(rlst_list, sig.cutoff = 0.05, w = 2.75, h = 3, OUT_DIR)
}
\arguments{
  \item{rlst_list}{A named list. Each element should contain a
    data.frame \code{rslt_of_interest} with at least the columns
    \code{adj.P.Val} (adjusted p-value) and \code{logFC}. Typical
    input is the output from \code{runLimmaCellFie} or
    \code{limmaTest_CellFie}.}
  \item{sig.cutoff}{Numeric. Adjusted p-value threshold used to call a
    task significant (default 0.05).}
  \item{w}{Numeric. Width of the output PDF in inches (default 2.75).}
  \item{h}{Numeric. Height of the output PDF in inches (default 3).}
  \item{OUT_DIR}{Character. Directory where the plot PDF
    (\code{SigDiffNr.pdf}) will be written.}
}
\details{
For each list element (named by cohort), the function counts rows in
\code{rslt_of_interest} with \code{adj.P.Val < sig.cutoff} separately
for positive and negative log fold changes. Counts for down-regulated
tasks are negated so the plot shows up/down directionality from a
single baseline (y = 0).
}
\value{Invisibly returns the ggplot object used to render the plot; the
main side-effect is the written PDF file.}
\author{Yuan Li {yuan.li@nbis.se}}
\examples{
\dontrun{
## Example usage (assuming 'rlst_list' is prepared):
SigNrBarplot(rlst_list, sig.cutoff = 0.05, w = 3, h = 4, OUT_DIR = "./plots/")
}
}
\keyword{hplot}
