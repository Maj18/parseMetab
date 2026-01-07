\name{visualizeGSsizeClass}
\alias{visualizeGSsizeClass}
\title{Visualize KEGG metabolic pathway sizes by class}
\description{
Summarize and plot the distribution of KEGG metabolic pathway sizes
  (number of member genes) grouped by pathway class. This function
  produces a pathway-size summary useful for comparing pathway gene
  counts across metabolic classes and datasets.
}
\usage{
visualizeGSsizeClass(GSVA_limma_rslt_gsva, genes4GSVA,
  kegg_metab_db, kegg_metab_db_table, OUTDIR)
}
\arguments{
  \item{GSVA_limma_rslt_gsva}{A named list of GSVA result matrices or data frames.
    Each element should correspond to a group (list names indicate groups,
    e.g. cancer types) and contain KEGG pathway identifiers as row names
    and sample submitter IDs as column names.}
  \item{genes4GSVA}{A character vector of gene identifiers used to compute GSVA
    scores (typically the genes present in the normalized expression matrix).}
  \item{kegg_metab_db}{A pathway-to-genes list mapping KEGG pathway names (or IDs)
    to character vectors of their member genes.}
  \item{kegg_metab_db_table}{A data frame or tibble with at least two columns:
    \code{gs_name} (pathway name or ID) and \code{class} (the higher-level metabolic class).
    This table is used to map pathways to their classes for grouping and plotting.}
  \item{OUTDIR}{Directory where output plots and summary tables will be written.
    If \code{NULL}, plots may be returned invisibly instead of saved.}
}
\value{
The function saves a PDF plot, Pathway_size_class.pdf, to OUTDIR.
}
\author{
Yuan Li {yuan.li@nbis.se}
}
\examples{
## Example usage (assuming appropriate objects are available):
## visualizeGSsizeClass(GSVA_limma_rslt_gsva, genes4GSVA, kegg_metab_db, kegg_metab_db_table, OUTDIR)
}

