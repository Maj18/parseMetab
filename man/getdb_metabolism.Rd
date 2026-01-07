\name{getdb_metabolism}
\alias{getdb_metabolism}
\title{
Extract KEGG Metabolic Pathway Database
}
\description{
Extract all the metabolic pathways from the KEGG database
}
\usage{
getdb_metabolism(databaseDIR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{databaseDIR}{a diretory to save the KEGG metabolic pathway database
}
}
\details{
Extract all the metabolic pathways from the KEGG database
}
\value{
\item{kegg_metab_db }{A list of metabolic pathways showing their gene members (symbols) separated by ';'}
\item{kegg_metab_db_table }{A table with three columns: gs_name (pathway), gene_symbol, class (pathway class)}
}
\references{

}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (databaseDIR) 
{
    query = KEGGREST::keggGet(c("hsa01100"))
    metab_pathways = names(query[[1]]$REL_PATHWAY)
    get_pathway_info = function(path_id) {
        info = KEGGREST::keggGet(path_id)[[1]]
        gene_ids = if (!is.null(info$GENE)) 
            info$GENE[seq(1, length(info$GENE), 2)]
        else character()
        gene_symbols = if (!is.null(info$GENE)) 
            info$GENE[seq(2, length(info$GENE), 2)] \%>\% strsplit(., 
                ";") \%>\% sapply(., function(g) g[1])
        else character()
        gene_infos = if (!is.null(info$GENE)) 
            info$GENE[seq(2, length(info$GENE), 2)] \%>\% strsplit(., 
                ";") \%>\% sapply(., function(g) g[2])
        else character()
        compound_ids = if (!is.null(info$COMPOUND)) 
            names(info$COMPOUND)
        else character()
        compound_Names = if (!is.null(info$COMPOUND)) 
            as.character(info$COMPOUND)
        else character()
        class = if (!is.null(info$CLASS)) 
            as.character(info$CLASS)
        else character()
        list(ID = as.character(info$ENTRY), pathway_name = info$NAME, 
            genes = gene_ids, gene_symbols = gene_symbols, gene_infos = gene_infos, 
            compounds = compound_ids, compound_Names = compound_Names, 
            class = class)
    }
    metabolism_Pathways = lapply(metab_pathways, function(path_id) get_pathway_info(path_id))
    metabolism_Pathways \%>\% lapply(., function(p) {
        if (length(p$genes) != 0) {
            data.frame(ID = p$ID, pathway_name = p$pathway_name \%>\% 
                gsub(" - Homo sapiens \\(human\\)", "", .), genes = paste(p$genes, 
                collapse = ";"), gene_symbols = paste(p$gene_symbols, 
                collapse = ";"), gene_infos = paste(p$gene_infos, 
                collapse = ";"), compounds = paste(p$compounds, 
                collapse = ";"), compound_Names = paste(p$compound_Names, 
                collapse = ";"), class = p$class)
        }
    }) \%>\% Filter(Negate(is.null), .) \%>\% Reduce(rbind, .) \%>\% 
        write.table(paste0(databaseDIR, "/KEGG_metabolism.csv"), 
            quote = F, sep = "\t", row.names = F)
    kegg_metab = read.delim(paste0(outdir, "/KEGG_metabolism/KEGG_metabolism.csv"), 
        header = TRUE, stringsAsFactors = FALSE) \%>\% mutate(class = gsub("Metabolism; ", 
        "", class))
    kegg_metab_db = stringr::str_split(kegg_metab[, "gene_symbols", 
        drop = T], ";") \%>\% setNames(., kegg_metab$pathway_name)
    path_class = setNames(kegg_metab$class, kegg_metab$pathway_name)
    kegg_metab_db_table = lapply(seq_along(kegg_metab_db), function(i) {
        p = kegg_metab_db[[i]]
        data.frame(gs_name = rep(names(kegg_metab_db)[i], length(p)), 
            gene_symbol = p, class = path_class[names(kegg_metab_db)[i]])
    }) \%>\% Reduce(rbind, .) \%>\% tibble::as_tibble()
    saveRDS(kegg_metab_db, paste0(databaseDIR, "/KEGG_metabolism_db.RDS"))
    saveRDS(kegg_metab_db_table, paste0(databaseDIR, "/KEGG_metabolism_db_table.RDS"))
    return(list(kegg_metab_db = kegg_metab_db, kegg_metab_db_table = kegg_metab_db_table))
  }
}

