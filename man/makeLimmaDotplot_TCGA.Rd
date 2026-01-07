\name{makeLimmaDotplot_TCGA}
\alias{makeLimmaDotplot_TCGA}
\title{Dot-plot summary of pathway-level limma results across TCGA cohorts}
\description{
Create a dot-plot summarizing pathway-level limma results across multiple
cohorts. Each dot represents a pathway–cohort pair sized by adjusted
p-value and colored by the limma t-statistic. Pathways can be grouped
by metabolic class using the provided \code{taskHierarchy} table and the
plot is saved as a PDF to \code{OUTDIRV}.
}
\usage{
makeLimmaDotplot_TCGA(rlst_list, adj.P.Val.cutoff, w = NULL, h = NULL, OUTDIRV, Depth = "gs_name", lowfigure = F, title = "Dysregulated pathways_TCGA", suffix = "TCGA", taskHierarchy = cfs, nested = T, top = FALSE, color_low = "skyblue", color_high = "purple")
}
\arguments{
  \item{rlst_list}{A list. Each element should be a data.frame (limma
    output) for a cohort and include at least the columns
    \code{logFC}, \code{adj.P.Val} and \code{t}.}
  \item{adj.P.Val.cutoff}{Numeric. Adjusted p-value cutoff used to select
    significant pathways for plotting.}
  \item{w}{Numeric or \code{NULL}. Plot width in inches; computed
    automatically when \code{NULL}.}
  \item{h}{Numeric or \code{NULL}. Plot height in inches; computed
    automatically when \code{NULL}.}
  \item{OUTDIRV}{Character. Directory where the generated PDF dotplot
    will be saved.}
  \item{Depth}{Character. Column name in \code{taskHierarchy} used as
    the pathway label (default: \code{"gs_name"}).}
  \item{lowfigure}{Logical. If \code{TRUE} use a compact layout with
    legends at the bottom.}
  \item{title}{Character. Plot title (default: \code{"Dysregulated pathways_TCGA"}).}
  \item{suffix}{Character. Filename suffix appended to the saved PDF
    (default: \code{"TCGA"}).}
  \item{taskHierarchy}{Data.frame. Mapping table that contains pathway
    metadata including the columns named by \code{Depth} and \code{class}.}
  \item{nested}{Logical. If \code{TRUE} add metabolic class information as
    nested facets (default: \code{TRUE}).}
  \item{top}{Logical. If \code{TRUE} show only the most frequently
    significant pathways (default: \code{FALSE}).}
  \item{color_low}{Character. Color used for negative t-statistics.
    Default \code{"skyblue"}.}
  \item{color_high}{Character. Color used for positive t-statistics.
    Default \code{"purple"}.}
}
\details{
The function concatenates limma results from \code{rlst_list},
filters rows by \code{adj.P.Val.cutoff}, annotates directionality
(Up/Down) using \code{logFC}, and joins pathway metadata from
\code{taskHierarchy}. The plot provides a compact overview of which
pathways are most frequently dysregulated across cohorts.
}
\value{Invisibly returns the ggplot object used to render the plot; the
main side-effect is the written PDF file.}
\author{Yuan Li {yuan.li@nbis.se}}
\examples{
\dontrun{
## Example (assumes 'rlst_list' and 'taskHierarchy' are prepared):
makeLimmaDotplot_TCGA(rlst_list, adj.P.Val.cutoff = 0.05,
                      OUTDIRV = "./plots/", suffix = "TCGA")
}
}
\keyword{hplot}
