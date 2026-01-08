\name{limmaTest_CellFie}
\alias{limmaTest_CellFie}
\title{Limma differential testing for CellFie-inferred metabolic activities}
\description{
Performs differential analysis on CellFie task scores using limma. The
function supports paired or unpaired designs and can include an optional
covariate alongside the primary `Condition` factor.
}
\usage{
limmaTest_CellFie(dat, paired = TRUE, currentCovariate = NULL, checkCovariate = FALSE, meta = meta, paired_variable = "Patient")
}
\arguments{
  \item{dat}{A numeric matrix or data.frame of CellFie scores (log2 transformed) with features as rownames and samples as column names. Rows are tested features (tasks), columns are samples.}
  \item{paired}{Logical scalar. If TRUE a paired analysis is fitted using the blocking column specified by `paired_variable`. Default: \code{TRUE}.}
  \item{currentCovariate}{Character scalar or \code{NULL}. Optional covariate name (a column in \code{meta}) to include in the model in addition to \code{Condition}. If \code{NULL} only \code{Condition} is used.}
  \item{checkCovariate}{Logical scalar. If \code{TRUE} perform basic checks on \code{currentCovariate} (existence in \code{meta} and not constant). Default: \code{FALSE}. if checkCovariate=TRUE, currentCovariate must be one a single covariate, this is for testing each covariate individually to decide whether to keep them in the model or not.}
  \item{meta}{A data.frame containing at least the columns \code{Sample} and \code{Condition}. \code{Sample} must match the column names of \code{dat}. Currently \code{Condition} is expected to have two categories.}
  \item{paired_variable}{Character scalar giving the column name in \code{meta} that links paired samples (for example a patient identifier). Default: \code{"Patient"}.}
}
\details{
The function filters features that do not have sufficient non-missing
observations across groups. When \code{paired = TRUE} a block model is
fitted using \code{limma::duplicateCorrelation} to estimate the within-block
correlation; the consensus correlation is then used in \code{limma::lmFit}.
}
\value{
A list with the following components:
\item{rslt_of_interest}{A data.frame (topTable) with statistics for the primary contrast of interest (the second level of \code{meta$Condition}).}
\item{rslt_co}{A data.frame topTable results for additional covariate coefficients (when \code{checkCovariate = TRUE}) or \code{NA} otherwise.}
\item{data}{The filtered input data matrix actually used for fitting (rows/features kept after filtering).}
\item{fitted.model}{The fitted limma model object (result of \code{limma::lmFit}).}
}

\examples{
## Example (illustrative only; user should supply real data):
\dontrun{
# dat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
# rownames(dat) <- paste0("Task", 1:100)
# meta <- data.frame(Sample = paste0("S", 1:10), Condition = rep(c("Normal", "Tumor"), each = 5), Patient = rep(1:5, 2))
# res <- limmaTest_CellFie(dat, paired = TRUE, meta = meta, paired_variable = "Patient")
# str(res)
}
}

\author{Yuan Li \email{yuan.li@nbis.se}}

\keyword{models}

