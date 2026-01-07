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

## The function is currently defined as
function (GSVA_limma_rslt_gsva, genes4GSVA, kegg_metab_db, kegg_metab_db_table, 
    OUTDIR) 
{
    pathway_class_mapping = setNames(kegg_metab_db_table$class, 
        kegg_metab_db_table$gs_name)
    P = lapply(GSVA_limma_rslt_gsva, function(g) {
        pathways = rownames(g)
        pathway_sizes = sapply(kegg_metab_db[pathways], function(d) {
            intersect(d, genes4GSVA) \%>\% length()
        })
        data.frame(Pathway = names(pathway_sizes), Psize = as.numeric(pathway_sizes), 
            Class = pathway_class_mapping[names(pathway_sizes)])
    }) \%>\% Reduce(rbind, .) \%>\% mutate(Class = gsub("etabolism", 
        ".", Class)) \%>\% mutate(Class = gsub("thesis", ".", Class)) \%>\% 
        mutate(Class = gsub("biodegradation", "biodegrad.", Class)) \%>\% 
        mutate(Class = gsub(" m.", "", Class)) \%>\% mutate(Class = gsub(" biosyn. and m.", 
        "", Class)) \%>\% mutate(Class = gsub("M. of ", "", Class)) \%>\% 
        mutate(Class = gsub(" and m.", "", Class)) \%>\% mutate(Class = gsub(" biosyn. and", 
        "", Class)) \%>\% mutate(Class = gsub(" biodegrad. and", 
        "", Class)) \%>\% mutate(Class = gsub("Biosyn. of ", "", 
        Class)) \%>\% ggplot2::ggplot(., ggplot2::aes(x = Class, 
        y = Psize)) + ggplot2::geom_boxplot(fill = "green4") + 
        ggplot2::labs(x = "", y = "Pathway size") + ggplot2::theme_minimal() + 
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, 
            hjust = 1, vjust = 0))
    pdf(paste0(OUTDIR, "/Pathway_size_class.pdf"), w = 2.5, h = 3.5)
    print(P)
    dev.off()
}    
