\name{visualizeGSVSscoresClass}
\alias{visualizeGSVSscoresClass}
\title{Visualize pathway-level GSVA means grouped by metabolic class}
\description{
Create boxplots that summarize mean GSVA pathway activities by metabolic
class. Two panels are produced: an overall cohort view and a normals-only view
computed from matched normal samples when available in \code{paired_data}.
}
\usage{
visualizeGSVSscoresClass(GSVA_limma_rslt, kegg_metab_db_table, paired_data, OUTDIR)
}
\arguments{
  \item{GSVA_limma_rslt}{A named list of GSVA result objects or matrices.
    Each element should correspond to a group (list names indicate groups,
    e.g. cancer types) and contain GSVA scores with pathways as row names
    and samples as column names (or provide component \code{gsva}).}
  \item{kegg_metab_db_table}{A data frame or tibble with at least two columns:
    \code{gs_name} (pathway name or ID) and \code{class} (the higher-level
    metabolic class). This table maps pathways to classes for grouping and plotting.}
  \item{paired_data}{A named list of metadata data frames (one per group) that
    correspond to \code{GSVA_limma_rslt}. Each metadata table must contain at
    least the columns \code{sample_type} (e.g. "Solid Tissue Normal" or
    "Primary tumor") and \code{sample.submitter_id} used to select normal samples.}
  \item{OUTDIR}{Character scalar or \code{NULL}. Directory where output plots and
    tables will be saved. If \code{NULL}, plots are returned instead of being written.}
}
\details{
The function computes mean GSVA scores per pathway for each group and maps
each pathway to its metabolic class via \code{kegg_metab_db_table}. It then
creates two boxplot panels (overall and normals-only) and writes a PDF to
\file{GSVA_pathway_class.pdf} in \code{OUTDIR} when provided.
}
\value{
A list with two \code{ggplot2} objects (overall and normals-only). Invisibly
returns the path to the saved PDF file when \code{OUTDIR} is supplied.
}
\examples{
## Example (assuming GSVA results and metadata are available):
## visualizeGSVSscoresClass(GSVA_limma_rslt, kegg_metab_db_table, paired_data, OUTDIR)
}

\author{
Yuan Li /email{yuan.li@nbis.se}
}


