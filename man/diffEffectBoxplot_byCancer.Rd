\name{diffEffectBoxplot_byCancer}
\alias{diffEffectBoxplot_byCancer}
	itle{Boxplot of limma differential effect (t statistic) by cancer}
\description{
Create a boxplot summarizing limma test statistics (t) for each cancer
cohort. The function concatenates task-level limma results from
\code{rlst_list}, orders cancer types by median t-statistic and writes
the resulting plot as a PDF to the specified output directory.
}
\usage{
diffEffectBoxplot_byCancer(rlst_list, OUT_DIR, w = 2.75, h = 3, title = "Pan-cancer tasks")
}
\arguments{
  \item{rlst_list}{A named list. Each element should contain a data.frame
    \code{rslt_of_interest} with at least the columns \code{t},
    \code{adj.P.Val} and \code{logFC}.}
  \item{OUT_DIR}{Character. Directory where the PDF file
    (\code{diffEffectBoxplot_cancer.pdf}) will be written.}
  \item{w}{Numeric. Width of the output PDF in inches (default 2.75).}
  \item{h}{Numeric. Height of the output PDF in inches (default 3).}
  \item{title}{Character. Plot title (default: "Pan-cancer tasks").}
}
\details{
The function labels each limma result row with its originating cohort,
optionally removes duplicate rows, computes ordering of cancer types by
median t-statistic and shows distributions of t-statistics per cancer
in a boxplot. This helps compare effect-size distributions across
cohorts.
}
\value{Invisibly returns the ggplot object used to render the plot; the
main side-effect is the written PDF file.}
\author{Yuan Li {yuan.li@nbis.se}}
\examples{
\dontrun{
## Example usage (assuming 'rlst_list' is prepared):
diffEffectBoxplot_byCancer(rlst_list, OUT_DIR = "./plots/", w = 3, h = 4,
                           title = "Study Cancer Types")
}
}
\keyword{hplot}
