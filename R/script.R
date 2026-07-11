 # Functions for inferring metabolic acitivities using GSVA:
 

 #' Limma differential analysis of GSVA-inferred KEGG metabolic activities
 #'
 #' Normalize raw counts, compute GSVA pathway enrichment scores for KEGG
 #' metabolic pathways, and perform differential testing using limma. The
 #' function supports optional feature filtering and paired study designs.
 #'
 #' @param dat A numeric matrix or data frame of raw transcriptomics counts with
 #'   features (genes) as row names and samples as column names.
 #' @param paired Logical scalar. Set to \code{TRUE} to fit a paired analysis
 #'   when paired samples are available (at least two pairs). Paired samples
 #'   must be linked by a shared identifier provided in \code{paired_variable}.
 #' @param Feature.filtering Logical scalar. If \code{TRUE}, features are
 #'   filtered to retain only those with counts > 5 in at least 15% of the
 #'   samples in each group (default: \code{TRUE}).
 #' @param currentCovariate Optional character. Name of an additional covariate
 #'   (column in \code{meta}) to include in the model in addition to \code{Group}.
 #' @param meta A data frame of sample metadata that contains at minimum the
 #'   columns \code{Sample} (sample IDs matching \code{colnames(dat)}) and
 #'   \code{Group} (a factor with two levels for the comparison). Additional
 #'   covariates may be provided and referenced via \code{currentCovariate}.
 #' @param paired_variable Character scalar. Name of the metadata column used to
 #'   identify paired samples (for example, patient or subject ID). This is used
 #'   when \code{paired = TRUE}.
 #' @param OUTDIR Character scalar. Directory where outputs (plots and tables)
 #'   will be written (e.g. normalized counts, GSVA scores, MDS plot).
 #' @param kegg_metab_db A list of KEGG metabolic pathways mapping pathway
 #'   names/IDs to vectors of member gene symbols (one of the outputs of
 #'   \code{getdb_metabolism}).
 #' @param min.count Numeric scalar. Minimum count threshold for feature filtering.
 #' 
 #' @importFrom dplyr %>%
 #'
 #' @return A list with components:
 #' \item{rslt_of_interest}{A data frame of limma topTable results for the
 #'   comparison of interest (including statistics and adjusted p-values).}
 #' \item{gsva}{The matrix of GSVA scores (pathways x samples).}
 #' \item{fitted.model}{The limma fit object (returned invisibly).}
 #' 
 #' @export
 #'
 GSVAlimmaTest = function(dat, paired = TRUE, Feature.filtering = TRUE,
        currentCovariate = NULL, meta = meta, paired_variable = "Patient",
        OUTDIR, kegg_metab_db, min.count=5) {
    meta = filter(meta, Sample%in%colnames(dat))
    # Make sure the sample are in the same order in both meta and dat:
    dat = dat[, meta$Sample]
    meta = meta %>% dplyr::mutate(across(where(is.factor), droplevels))
    mapping = setNames(meta$Group, meta$Sample)
    min.nr = round(max((ncol(dat)/2)*0.15, 2), 0)
    print(paste0("Only keep features that have at least ", min.nr, 
                        " values that are >", min.count, " in each group..."))
    if (Feature.filtering) {
        B = mapping[colnames(dat)] %>% as.character()
        featureskeep = apply(dat>min.count, 1, function(row) {
            A = row
            sum(aggregate(A~B, data=data.frame(A=A,B=B), sum, na.rm=T)[,"A"]>=min.nr)>1
        })
        dat = dat[featureskeep, ]
    }

    if (nrow(dat)>0) {
        # conda install conda-forge::r-locfit
        # conda install bioconda::bioconductor-edger
        # BiocManager::install("edgeR")
        Group = factor(meta$Group)
        # Create a DEGList from the edgeR package
        dge = edgeR::DGEList(counts=dat)
        # Normalize the data
        dge = edgeR::calcNormFactors(dge)
        # Note: calcNormFactors does not normalize the data, it just calculates normalization factors for use downstream.
        pdf(paste0(OUTDIR, "/MDS_plot.pdf"))
            print(limma::plotMDS(dge, col=as.numeric(Group)))
        dev.off()

        #if (!is.null(currentCovariate)) {currentCovariate=NULL}
        baseFormula = ~ Group
        if (is.null(currentCovariate)) {
            currentFormula = baseFormula
        } else {
            currentFormula = as.formula(paste0("~ ", currentCovariate, " + Group"))}
        
        design = model.matrix(currentFormula, data=meta)
        colnames(design) = gsub("Group", "", colnames(design))
        pdf("Temp.pdf")
        # options(bitmapType = "png")
        dat2 = limma::voom(dge, design, plot=T)
        dev.off()
        # What is voom doing
        # Counts are transformed to log2 counts per million reads (CPM), 
        # where per million reads is defined based on the normalization factors we calculated earlier
        # A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
        # A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
        # The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs
        # Save normalized count
        normalizedCount = dat2$E
        write.table(
            normalizedCount %>% as.data.frame() %>%
                    tibble::rownames_to_column(var="Features"),
            paste0(OUTDIR, "/NormalizedCount.csv"),
            quote = F, row.names = F, col.names = T, sep="\t")
        
        # Calcualte GSVA per sample
        # BiocManager::install("GSVA")
        normalizedCount = as.matrix(normalizedCount)
        normalizedCount = matrix(as.numeric(normalizedCount),
            nrow = nrow(normalizedCount),
            ncol = ncol(normalizedCount), dimnames = dimnames(normalizedCount))

        param = GSVA::gsvaParam(normalizedCount, kegg_metab_db, minSize=5, maxSize=500)
        gsva_scores = tryCatch(GSVA::gsva(param, verbose=TRUE) %>% 
                                as.data.frame(), error = function(e) NA)
        if (class(gsva_scores)=="data.frame"&length(unique(meta$Patient))>1) {
            write.table(
                gsva_scores %>% as.data.frame() %>%
                        tibble::rownames_to_column(var="Features"),
                paste0(OUTDIR, "/GSVAscores.csv"),
                quote = F, row.names = F, col.names = T, sep="\t")
            if (paired) {
                #Fit with correlated arrays
                dupcor = limma::duplicateCorrelation(gsva_scores, design,
                    block=meta[,paired_variable,drop=T])
                fit = limma::lmFit(gsva_scores, design, block=meta[,paired_variable,drop=T],
                    correlation=dupcor$consensus)
            } else { fit = limma::lmFit(gsva_scores, design) }

            # Limma trend to refine gene-level variance estimates
            # gssizes = sapply(kegg_metab_db[rownames(fit$coefficients)], function(gs) {
            #     length(gs)
            # })
            # fit.con = limma::eBayes(fit, robust = TRUE, trend=gssizes)
            fit.con = limma::eBayes(fit, robust = TRUE, trend=TRUE)
            rlst_interest =
                limma::topTable(fit.con, n=Inf, coef=levels(meta$Group)[2]) %>%
                arrange(-t)
            write.table(
                rlst_interest %>% as.data.frame() %>%
                        tibble::rownames_to_column(var="Features"),
                paste0(OUTDIR, "/limmaResult.csv"),
                quote = F, row.names = F, col.names = T, sep="\t")
        } else {
            fit = NULL
            rlst_interest = NULL
        }
    } else {
        fit = NULL
        rlst_interest = NULL
    }
    return(list(rslt_of_interest=rlst_interest, gsva=gsva_scores, fitted.model=fit))
}



#' @title Get KEGG Metabolic Pathway Database
#' @description
#' Retrieves human metabolic pathways from the KEGG database, extracts pathway
#' information including genes and compounds, and creates a structured database
#' for metabolic analysis.
#'
#' @param databaseDIR Character string specifying the directory path where the
#'   metabolic database files will be saved (KEGG_metabolism.csv,
#'   KEGG_metabolism_db.RDS, KEGG_metabolism_db_table.RDS).
#'
#' @return A list containing two elements:
#'   \item{kegg_metab_db}{Named list of gene symbols split by pathway name}
#'   \item{kegg_metab_db_table}{Tibble with columns: gs_name, gene_symbol, class}
#'
#' @details
#' This function queries the KEGG REST API to obtain human metabolic pathways,
#' extracts genes and compounds for each pathway, and creates both wide (list)
#' and long (table) format databases. Results are saved as TSV and RDS files.
#'
#' @importFrom KEGGREST keggGet
#' @importFrom dplyr mutate
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#'   result <- getdb_metabolism("/path/to/database")
#'   head(result$kegg_metab_db_table)
#' }
#'
#' @export
getdb_metabolism = function(databaseDIR) {
    # Get all human metabolic pathways
    query = KEGGREST::keggGet(c("hsa01100"))
    metab_pathways = names(query[[1]]$REL_PATHWAY)

    get_pathway_info = function(path_id) {
      info = KEGGREST::keggGet(path_id)[[1]]
      gene_ids = if (!is.null(info$GENE)) info$GENE[seq(1, length(info$GENE), 2)] else character()
      gene_symbols = if (!is.null(info$GENE)) info$GENE[seq(2, length(info$GENE), 2)]%>%
            strsplit(., ";") %>% sapply(., function(g) g[1]) else character()
      gene_infos = if (!is.null(info$GENE)) info$GENE[seq(2, length(info$GENE), 2)]%>%
            strsplit(., ";") %>% sapply(., function(g) g[2]) else character()
      compound_ids = if (!is.null(info$COMPOUND)) names(info$COMPOUND) else character()
      compound_Names = if (!is.null(info$COMPOUND)) as.character(info$COMPOUND) else character()
      class = if (!is.null(info$CLASS)) as.character(info$CLASS) else character()
      list(
        ID = as.character(info$ENTRY),
        pathway_name = info$NAME,
        genes = gene_ids,
        gene_symbols = gene_symbols,
        gene_infos = gene_infos,
        compounds = compound_ids,
        compound_Names = compound_Names,
        class = class
      )
    }

    metabolism_Pathways =
        lapply(metab_pathways, function(path_id) get_pathway_info(path_id))

    # Remove empty pathways
    metabolism_Pathways %>% lapply(., function(p) {
        if(length(p$genes)!=0) {
            data.frame(
            ID = p$ID,
            pathway_name = p$pathway_name %>% gsub(" - Homo sapiens \\(human\\)", "", .),
            genes = paste(p$genes, collapse=";"),
            gene_symbols = paste(p$gene_symbols, collapse=";"),
            gene_infos = paste(p$gene_infos, collapse=";"),
            compounds = paste(p$compounds, collapse=";"),
            compound_Names = paste(p$compound_Names, collapse=";"),
            class = p$class)
    }}) %>% Filter(Negate(is.null), .) %>%
        Reduce(rbind, .) %>%
        write.table(paste0(databaseDIR, "/KEGG_metabolism.csv"),
            quote=F, sep="\t", row.names=F)

    # Create KEGG metabolic database
    kegg_metab = read.delim(
        paste0(databaseDIR, "/KEGG_metabolism.csv"),
        header = TRUE, stringsAsFactors = FALSE) %>%
        dplyr::mutate(class = gsub("Metabolism; ", "", class))
    kegg_metab_db =
        stringr::str_split(kegg_metab[,"gene_symbols",drop=T], ";") %>%
        setNames(., kegg_metab$pathway_name)

    path_class = setNames(kegg_metab$class, kegg_metab$pathway_name)
    kegg_metab_db_table = lapply(seq_along(kegg_metab_db), function(i) {
        p = kegg_metab_db[[i]]
        data.frame(gs_name=rep(names(kegg_metab_db)[i],length(p)),
            gene_symbol = p, class=path_class[names(kegg_metab_db)[i]])
        }) %>% Reduce(rbind, .) %>% tibble::as_tibble()
    saveRDS(kegg_metab_db, paste0(databaseDIR, "/KEGG_metabolism_db.RDS"))
    saveRDS(kegg_metab_db_table, paste0(databaseDIR, "/KEGG_metabolism_db_table.RDS"))

    return(list(kegg_metab_db=kegg_metab_db, kegg_metab_db_table=kegg_metab_db_table))

}


 #' Circular packing plot of KEGG metabolic pathway hierarchy
 #'
 #' Create a circular packing plot that summarizes KEGG metabolic pathways and
 #' their grouping by metabolic class. Labels are shortened for display. The plot is saved as a PDF in
 #' the provided output directory.
 #'
 #' @param kegg_metab_db_table A data.frame or tibble with at least three
 #'   columns: \code{gs_name} (pathway name or ID), \code{gene_symbol}
 #'   (member gene symbol), and \code{class} (pathway class).
 #' @param outdir Character scalar. Directory where the plot PDF will be written.
 #'   The file is saved to \code{file.path(outdir, "KEGG_metabolism", "TaskSummary.pdf")}.
 #' 
 #' @importFrom dplyr %>%
 #'
 #' @return Invisibly returns the generated \code{ggplot2} object. The primary
 #'   side-effect is writing the PDF file to \code{outdir} when provided.
 #'
 #' @export
 #' 
 getCirPackingPlot_GSVA = function(kegg_metab_db_table, outdir) {
    pathway_hierachy =  kegg_metab_db_table %>% mutate(Depth0="Metabolism") %>%
        dplyr::mutate(Depth2=gs_name) %>% dplyr::mutate(Depth1=class) %>%
        dplyr::select(Depth0, Depth1, Depth2) %>% unique()
    edges = pathway_hierachy %>% .[,c("Depth0","Depth1")] %>%
        unique()%>%setNames(c("from", "to")) %>%
        rbind(pathway_hierachy[, c("Depth1","Depth2")] %>%
            unique() %>% setNames(c("from", "to")))
    vertices = pathway_hierachy %>% .[, c("Depth0","Depth2")] %>%
        unique()%>%pull(Depth0)%>%table()%>%as.data.frame()%>%setNames(c("name", "size")) %>%
        rbind(pathway_hierachy[,c("Depth1","Depth2")]%>%unique()%>%pull(Depth1)%>%table()%>%
            as.data.frame()%>%setNames(c("name", "size"))) %>%
        rbind(pathway_hierachy %>% pull(Depth2) %>% unique() %>% table() %>%
            as.data.frame() %>% setNames(c("name","size"))) %>%
            dplyr::mutate(shortName = gsub(" and metabolism", "", paste0(name,"(",size,")"))) %>%
            dplyr::mutate(shortName = gsub(" metabolism", "", shortName)) %>%
            dplyr::mutate(shortName = gsub("Metabolism of ", "", shortName)) %>%
            dplyr::mutate(shortName = gsub(" biosynthesis and", "", shortName)) %>%
            dplyr::mutate(shortName = gsub(" biosynthesis", "", shortName)) %>%
            dplyr::mutate(shortName = gsub(" biodegradation and", "", shortName)) %>%
            dplyr::mutate(shortName = gsub("Biosynthesis of ", "", shortName)) %>%
            dplyr::mutate(shortName = stringr::str_to_title(shortName))
    vertices$shortName[!vertices$name%in%unique(pathway_hierachy$Depth1)] = NA
    mygraph = igraph::graph_from_data_frame(edges, vertices=vertices)
    p = ggraph::ggraph(mygraph, layout = 'circlepack', weight=size) +
        ggraph::geom_node_circle(aes(col=depth)) +
        ggraph::geom_node_label(aes(label=shortName),repel=T,size=4,color="darkgreen") +
        ggplot2::theme_void() +
        viridis::scale_color_viridis() +
        ggplot2::theme(legend.position="FALSE") +
        ggplot2::ggtitle("                                           TCGA")
    dir.create(paste0(outdir, "/KEGG_metabolism/"), recursive=TRUE, showWarnings=FALSE)
    pdf(paste0(outdir, "/KEGG_metabolism/TaskSummary.pdf"), h=5, w=4.85)
        print(p)
    dev.off()
}


 #' @title Visualize GSVA pathway means across groups (and normals)
 #'
 #' @description Generate side-by-side boxplots that summarize GSVA pathway activity means
 #' across provided groups and, separately, for matched normal samples.
 #' This is intended for cohort-level comparisons (for example, across cancer types)
 #' and assumes GSVA has been run with KEGG pathway gene sets.
 #'
 #' @param GSVA_limma_rslt_gsva A named list of GSVA result matrices or data frames.
 #'   Each element should correspond to a group (list names indicate the groups,
 #'   e.g. cancer types) and contain pathways (features) as row names and
 #'   sample submitter IDs as column names.
 #' @param paired_data A named list of metadata data frames (one per group) that
 #'   corresponds to `GSVA_limma_rslt_gsva`. Each metadata table must contain at
 #'   least the columns `sample_type` ("Solid Tissue Normal" or "Primary tumor")
 #'   and `sample.submitter_id` which are used to select normal samples.
 #' @param OUTDIR Character scalar. Directory where output plots (PDF) will be
 #'   written. If `NULL`, plots may be returned instead of being saved.
  #' @param Condition_column Character string indicating the column name in `paired_data`
 #'   data frames used to classify samples (e.g., "sample_type"). Default: "sample_type".
 #' @param Condition.control Character string specifying the value in `Condition_column`
 #'   that identifies control/normal samples (e.g., "Solid Tissue Normal").
 #'   Default: "Solid Tissue Normal".
 #' @param sample_column Character string indicating the column name in `paired_data`
 #'   data frames that contains sample identifiers matching the column names in
 #'   GSVA matrices (e.g., "sample.submitter_id"). Default: "sample.submitter_id".
 #' @param w Numeric scalar. Width of the output PDF in inches. Default: 5.5.
 #' @param h Numeric. Height of the output PDF in inches. Default: 3.0.
 #' @param title Character string. Title to display on the plots. Default: "TCGA".
 #'
 #' @importFrom dplyr %>% mutate group_by summarise arrange pull
 #' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal theme element_text
 #' 
 #' @return Invisibly returns a list containing the two ggplot2 objects (all
 #'   samples and normal-only samples) and the path to the saved PDF (if
 #'   `OUTDIR` is provided). A PDF named `GSVA_samples.pdf` is written to
 #'   `OUTDIR` when `OUTDIR` is non-NULL.
 #' @export
 #' 
 visualizeGSVSscoresGroup = function(
    GSVA_limma_rslt_gsva, paired_data, OUTDIR,
    Condition_column = "sample_type",
    Condition.control = "Solid Tissue Normal",
    sample_column = "sample.submitter_id",
    w=5.5, h=3.0, title="TCGA",
    plot_healthy_control=TRUE) {
    # GSVA scores for all samples
    dat = lapply(seq_along(GSVA_limma_rslt_gsva), function(i) {
        data.frame(
            means = (
                rowSums(GSVA_limma_rslt_gsva[[i]],na.rm=TRUE)/ncol(GSVA_limma_rslt_gsva[[i]])) %>%
                    as.numeric(),
            Cancer_type = names(GSVA_limma_rslt_gsva)[i])
        }) %>% Reduce(rbind, .)
    Cancer_type_order = dat %>% group_by(Cancer_type) %>%
            summarise(Means=mean(means)) %>%
            arrange(Means) %>% pull(Cancer_type) %>% unique()
    P1 = dat %>% dplyr::mutate(Cancer_type = factor(Cancer_type, levels=Cancer_type_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Cancer_type, y=means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(title = title,
                x="Cancer_type", y="Mean Metab. Act.")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))
    if (plot_healthy_control) {
        # GSVA scores for healthy control samples only
        dat = lapply(seq_along(GSVA_limma_rslt_gsva), function(j) {
            normal_samples = paired_data[[j]] # %>%
                # filter(sample_type=="Solid Tissue Normal") %>% pull(sample.submitter_id)
            normal_samples = normal_samples[normal_samples[[Condition_column]]==Condition.control, ]
            normal_samples = normal_samples[[sample_column]]
            gsva = GSVA_limma_rslt_gsva[[j]] %>% .[, normal_samples]
            data.frame(
                means = (rowSums(gsva, na.rm=TRUE)/ncol(gsva)) %>% as.numeric(),
                Cancer_type = names(GSVA_limma_rslt_gsva)[j])
            }) %>% Reduce(rbind, .)
        Cancer_type_order = dat %>% group_by(Cancer_type) %>%
                summarise(Means=mean(means)) %>%
                arrange(Means) %>% pull(Cancer_type) %>% unique()
        P2 = dat %>% dplyr::mutate(Cancer_type = factor(Cancer_type, levels=Cancer_type_order)) %>%
                ggplot2::ggplot(., ggplot2::aes(x=Cancer_type, y=means)) +
                ggplot2::geom_boxplot(fill="green4") +
                ggplot2::labs(title = title,
                    x="Cancer_type", y="Mean Norm. Metab. Act.")+
                ggplot2::theme_minimal() +
                ggplot2::theme(legend.position="none",
                    axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))
    } else {P2=P1}

    dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
    pdf(paste0(OUTDIR, "/GSVA_samples.pdf"), w=w, h=h)
        print(cowplot::plot_grid(P1, P2))
    dev.off()
}


 #' Visualize KEGG metabolic pathway sizes by class
 #'
 #' Summarize and plot the distribution of KEGG metabolic pathway sizes (number of member genes)
 #' grouped by pathway class. This function produces a pathway-size summary useful for
 #' comparing pathway gene counts across metabolic classes and datasets.
 #'
 #' @param GSVA_limma_rslt_gsva A named list of GSVA result matrices (or data frames). Each list
 #'   element should correspond to a group (the list names indicate the groups, e.g. cancer
 #'   types) and contain KEGG pathway identifiers as row names and sample submitter IDs as
 #'   column names.
 #' @param genes4GSVA A character vector of gene identifiers used to compute GSVA scores
 #'   (typically the genes present in the normalized expression matrix).
 #' @param kegg_metab_db A pathway-to-genes list mapping KEGG pathway names (or IDs) to
 #'   character vectors of their member genes.
 #' @param kegg_metab_db_table A data frame or tibble with at least two columns:
 #'   \code{gs_name} (pathway name or ID) and \code{class} (the higher-level metabolic class).
 #'   This table is used to map pathways to their classes for grouping and plotting.
 #' @param OUTDIR Directory where output plots/tables will be written (optional).
 #'
 #' @importFrom dplyr %>%
 #' 
 #' @return A pdf file, Pathway_size_class.pdf in the OUTDIR.
 #' @examples
 #' # visualizeGSsizeClass(GSVA_limma_rslt_gsva, genes4GSVA, kegg_metab_db, kegg_metab_db_table, OUTDIR)
 #' @export
 #' 
 visualizeGSsizeClass = function(
        GSVA_limma_rslt_gsva, genes4GSVA,
        kegg_metab_db, kegg_metab_db_table,
        OUTDIR) {
    pathway_class_mapping =
        setNames(kegg_metab_db_table$class, kegg_metab_db_table$gs_name)
    P = lapply(GSVA_limma_rslt_gsva, function(g) {
        # get pathways
        pathways = rownames(g)
        pathway_sizes = sapply(kegg_metab_db[pathways],function(d) {
            intersect(d, genes4GSVA) %>% length()
        })
        data.frame(Pathway = names(pathway_sizes),
                   Psize = as.numeric(pathway_sizes),
                   Class = pathway_class_mapping[names(pathway_sizes)])
    }) %>% Reduce(rbind, .) %>%
        dplyr::mutate(Class = gsub(" and metabolism", "", Class)) %>%
        dplyr::mutate(Class = gsub(" metabolism", "", Class)) %>%
        dplyr::mutate(Class = gsub("Metabolism of ", "", Class)) %>%
        dplyr::mutate(Class= gsub(" biosynthesis and", "", Class)) %>%
        dplyr::mutate(Class = gsub(" biosynthesis", "", Class)) %>%
        dplyr::mutate(Class = gsub(" biodegradation and", "", Class)) %>%
        dplyr::mutate(Class = gsub("Biosynthesis of ", "", Class)) %>%
        dplyr::mutate(Class = stringr::str_to_title(Class)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=Class, y=Psize)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(# title = "TCGA",
            x="", y="Pathway size")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))
    pdf(paste0(OUTDIR, "/Pathway_size_class.pdf"), w=2.5, h=3.5)
        print(P)
    dev.off()
}

 #' @title Visualize GSVA Metabolic Activity by Pathway Class
 #'
 #' @description
 #' Generates comparative boxplots of mean GSVA pathway activity scores grouped
 #' by metabolic class across multiple cohorts or groups. The function creates
 #' two complementary panels: one showing overall cohort means and another
 #' showing means computed from matched normal samples only. Pathway class names
 #' are automatically cleaned and standardized for improved readability.
 #'
 #' @param GSVA_limma_rslt A named list of GSVA result objects. Each list
 #'   element should correspond to a group (list names indicate groups, e.g.
 #'   cancer types) and contain at least a component `gsva` with pathways as
 #'   row names and samples as column names, or be a matrix/data.frame of
 #'   GSVA scores (pathways x samples).
 #' @param kegg_metab_db_table A data.frame or tibble with at least two columns:
 #'   \code{gs_name} (pathway name or ID) and \code{class} (the higher-level
 #'   metabolic class). This table is used to map pathways to their classes for
 #'   grouping and plotting.
 #' @param paired_data A named list of metadata data.frames (one per group) that
 #'   correspond to \code{GSVA_limma_rslt}. Each metadata table must contain at
 #'   least the columns \code{sample_type} (e.g. "Solid Tissue Normal" or
 #'   "Primary tumor") and \code{sample.submitter_id}, which are used to select
 #'   normal samples.
 #' @param OUTDIR Character scalar or \code{NULL}. Directory where output
 #'   plots/tables will be written. If \code{NULL}, plots are returned instead
 #'   of being saved.
 #' @param Condition_column Character string naming the metadata column that
 #'   contains sample classification (e.g., "sample_type"). Default: "sample_type".
 #'
 #' @param Condition.control Character string specifying the value in
 #'   `Condition_column` that identifies control/normal samples
 #'   (e.g., "Solid Tissue Normal"). Default: "Solid Tissue Normal".
 #' @param sample_column Character string naming the metadata column containing
 #'   sample identifiers that match column names in the GSVA matrices
 #'   (e.g., "sample.submitter_id"). Default: "sample.submitter_id".
 #' @param w Numeric. Width of output PDF in inches. Default: 5.5.
 #' @param h Numeric. Height of output PDF in inches. Default: 3.5.
 #'
 #' @details
 #' The function performs the following operations:
 #' \enumerate{
 #'   \item Computes mean GSVA scores per pathway across all samples in each group
 #'   \item Maps pathways to metabolic classes using `kegg_metab_db_table`
 #'   \item Cleans class names by removing common KEGG descriptors
 #'     ("metabolism", "biosynthesis", "biodegradation", etc.)
 #'   \item Creates two boxplots (all samples vs. control samples only),
 #'     sorted by median pathway activity
 #'   \item Combines both plots horizontally using cowplot for PDF output
 #' }
 #' Class name standardization helps create more readable axis labels
 #' across different metabolic categories.
 #'
 #' @importFrom dplyr %>% select group_by summarise arrange pull mutate
 #' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal theme element_text
 #' @importFrom stringr str_to_title
 #'
 #' @return Invisibly returns a list with the ggplot2 objects for the two
 #'   panels (overall and normals-only) and, when \code{OUTDIR} is provided,
 #'   the path to the saved PDF file. The primary side-effect is creation of a
 #'   PDF file in \code{OUTDIR} named \file{GSVA_pathway_class.pdf}.
 #'
 #' @export
 #' 
 visualizeGSVSscoresClass = function(
    GSVA_limma_rslt, kegg_metab_db_table, paired_data, OUTDIR,
    Condition_column = "sample_type",
    Condition.control = "Solid Tissue Normal",
    sample_column = "sample.submitter_id",
    w=5.5, h=3.5) {
    dat = lapply(seq_along(GSVA_limma_rslt), function(i) {
        s = names(GSVA_limma_rslt)[i]
        pathway_means = rowSums(GSVA_limma_rslt[[i]]$gsva,na.rm=TRUE)/ncol(GSVA_limma_rslt[[i]]$gsva)
        kegg_metab_db_table =
            kegg_metab_db_table %>% dplyr::select(gs_name, class) %>%
            unique()
        kegg_metab_db_mapping = setNames(kegg_metab_db_table$class, kegg_metab_db_table$gs_name)
        data.frame(Pathway = names(pathway_means),
            Means = as.numeric(pathway_means),
            Class = kegg_metab_db_mapping[names(pathway_means)]) %>%
            dplyr::mutate(Class = gsub(" and metabolism", "", Class)) %>%
            dplyr::mutate(Class = gsub(" metabolism", "", Class)) %>%
            dplyr::mutate(Class = gsub("Metabolism of ", "", Class)) %>%
            dplyr::mutate(Class= gsub(" biosynthesis and", "", Class)) %>%
            dplyr::mutate(Class = gsub(" biosynthesis", "", Class)) %>%
            dplyr::mutate(Class = gsub(" biodegradation and", "", Class)) %>%
            dplyr::mutate(Class = gsub("Biosynthesis of ", "", Class)) %>%
            dplyr::mutate(Class = stringr::str_to_title(Class))
        }) %>% Reduce(rbind, .)
        Class_order = dat %>% group_by(Class) %>%
            summarise(Means=mean(Means)) %>% arrange(Means) %>% pull(Class) %>% unique()
    P1 = dat %>% dplyr::mutate(Class = factor(Class, levels=Class_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Class, y=Means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(# title = "TCGA",
                x="", y="Mean Metab. Act.")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    dat = lapply(seq_along(GSVA_limma_rslt), function(j) {
        s = names(GSVA_limma_rslt)[j]
        normal_samples = paired_data[[s]] # %>%
            # filter(sample_type=="Solid Tissue Normal") %>% pull(sample.submitter_id)
        normal_samples = normal_samples[normal_samples[[Condition_column]]==Condition.control, ]
        normal_samples = normal_samples[[sample_column]]
        gsva = GSVA_limma_rslt[[s]]$gsva %>% .[, normal_samples]
        pathway_means = rowSums(gsva,na.rm=TRUE)/ncol(gsva)
        kegg_metab_db_table =
            kegg_metab_dbs$kegg_metab_db_table %>% dplyr::select(gs_name, class) %>%
            unique()
        kegg_metab_db_mapping = setNames(kegg_metab_db_table$class, kegg_metab_db_table$gs_name)
        data.frame(Pathway = names(pathway_means),
            Means = as.numeric(pathway_means),
            Class = kegg_metab_db_mapping[names(pathway_means)]) %>%
            dplyr::mutate(Class = gsub(" and metabolism", "", Class)) %>%
            dplyr::mutate(Class = gsub(" metabolism", "", Class)) %>%
            dplyr::mutate(Class = gsub("Metabolism of ", "", Class)) %>%
            dplyr::mutate(Class= gsub(" biosynthesis and", "", Class)) %>%
            dplyr::mutate(Class = gsub(" biosynthesis", "", Class)) %>%
            dplyr::mutate(Class = gsub(" biodegradation and", "", Class)) %>%
            dplyr::mutate(Class = gsub("Biosynthesis of ", "", Class)) %>%
            dplyr::mutate(Class = stringr::str_to_title(Class))
    }) %>% Reduce(rbind, .)
    Class_order = dat %>% group_by(Class) %>%
        summarise(Means=mean(Means)) %>% arrange(Means) %>% pull(Class) %>% unique()
    P2 = dat %>% dplyr::mutate(Class = factor(Class, levels=Class_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Class, y=Means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(# title = "TCGA",
                x="", y="Mean Norm. Metab.")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUTDIR, "/GSVA_pathway_class.pdf"), w=w, h=h)
        print(cowplot::plot_grid(P1, P2))
    dev.off()
}


 #' @title Visualize Significant Pathways Count by Group
 #'
 #' @description
 #' Generates a stacked bar chart showing the count of significantly dysregulated
 #' pathways (up- and down-regulated) across multiple groups (e.g., cancer types).
 #' Significance is determined by statistical threshold applied to differential
 #' expression results from limma testing on GSVA pathway activity scores.
 #' The visualization helps identify which groups have the most pathway dysregulation.
 #'
 #' @param OUTDIR Character scalar. Directory where the output PDF
 #'   (\file{Sig_pathway_nr.pdf}) will be written.
 #' @param GSVA_limma_rslt A named list of results returned by
 #'   \code{GSVAlimmaTest()} (one element per group). Each element must
 #'   contain a component \code{rslt_of_interest} (a data frame produced by
 #'   \code{limma::topTable}) with at least \code{logFC} and \code{adj.P.Val}
 #'   columns.
 #' @param significance_statistic Character string naming the column in the
 #'   result data frames to use for significance filtering
 #'   (e.g., "adj.P.Val", "p.value"). Default: "adj.P.Val".
 #' @param sig.cutoff Numeric. Significance threshold for pathway selection.
 #'   Only pathways with `significance_statistic < sig.cutoff` are counted
 #'   as significant. Default: 0.05.
 #' @param w Numeric. Width of output PDF in inches. Default: 3.
 #' @param h Numeric. Height of output PDF in inches. Default: 3.5.
 #'
 #' @details
 #' The function:
 #' \enumerate{
 #'   \item Filters pathways for each group using the specified significance threshold
 #'   \item Counts up-regulated (logFC > 0) and down-regulated (logFC < 0) pathways
 #'   \item Down-regulated counts are displayed as negative values
 #'   \item Creates a stacked bar chart with groups on x-axis and
 #'     up-regulated pathways in red (firebrick) and down-regulated in blue (steelblue)
 #'   \item Groups are ordered by total significant pathway count
 #' }
 #'
 #' @importFrom dplyr %>% filter arrange mutate
 #' @importFrom tidyr gather
 #' @importFrom tibble rownames_to_column
 #' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual geom_hline labs theme_minimal theme element_text
 #'
 #' @return Invisibly returns NULL. A PDF file named "Sig_pathway_nr.pdf" is
 #'   written to `OUTDIR` containing the stacked bar chart.
 #'
 #' @export
 #' 
 visualizeSigNr = function(
    OUTDIR, GSVA_limma_rslt, 
    significance_statistic="adj.P.Val", 
    sig.cutoff=0.05, w=3, h=3.5) {
    pdf(paste0(OUTDIR, "/Sig_pathway_nr.pdf"), w=w, h=h)
        sig_nr = lapply(names(GSVA_limma_rslt), function(cancer) {
            rslt = GSVA_limma_rslt[[cancer]]$rslt_of_interest
            rslt = rslt[rslt[[significance_statistic]]<sig.cutoff, ]
            data.frame(up = rslt %>% filter(logFC>0) %>% nrow(),
                down = -(rslt %>% filter(logFC<0) %>% nrow())) %>%
                t() %>% as.data.frame() %>% setNames(cancer)
            }) %>% Reduce(cbind, .) %>%
                tibble::rownames_to_column(var="Regulation") %>%
                tidyr::gather(key="Cancer_type", value="Count", -Regulation) %>%
                arrange(-Count) %>%
                dplyr::mutate(Cancer_type=factor(Cancer_type,levels=unique(Cancer_type)))
        print(ggplot2::ggplot(sig_nr, ggplot2::aes(x=Cancer_type, y=Count, fill=Regulation)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values=c("up"="firebrick","down"="steelblue")) +
            ggplot2::geom_hline(yintercept=0, color="black") +
            ggplot2::labs(title = "Sig. diff. systems",
                x="Cancer_type", y="Sig. Count")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()
}


#' @title Create Dot-Plot of Dysregulated Metabolic Pathways Across Cancer Types
#'
#' @description
#' Generates a comprehensive dot-plot summarizing pathway-level differential
#' expression results (from limma) across multiple cancer cohorts in TCGA.
#' Each dot represents a pathway-cohort pair, positioned by cancer type and pathway,
#' with size representing significance level and color representing the direction
#' and magnitude of dysregulation (t-statistic). Pathways can be hierarchically
#' organized by metabolic class with nested faceting.
#'
#' @param sig_statistic Character string naming the significance column in
#'   `rlst_list` data frames. Default: "adj.P.Val".
#'
#' @param sig.cutoff Numeric. Significance threshold for selecting significant
#'   pathways. Only pathways with `sig_statistic < sig.cutoff` are displayed.
#'
#' @param adj.P.Val.cutoff Numeric or NULL. Deprecated alias for `sig.cutoff`.
#'   If provided, overrides `sig.cutoff`. Default: NULL.
#'
#' @param w Numeric or NULL. Plot width in inches. When NULL, width is
#'   automatically computed based on number of cohorts. Default: NULL.
#'
#' @param h Numeric or NULL. Plot height in inches. When NULL, height is
#'   automatically computed based on number of pathways. Default: NULL.
#'
#' @param OUTDIRV Character string specifying the output directory where the
#'   PDF will be saved. Must be a valid, writable directory path.
#'
#' @param Depth Character string naming the column in `taskHierarchy` to use
#'   as the pathway label displayed on the y-axis.
#'   Default: "gs_name".
#'
#' @param lowfigure Logical. When TRUE, produces a compact layout with
#'   legends positioned below and horizontally aligned. Default: FALSE.
#'
#' @param title Character string for the plot title.
#'   Default: "Dysregulated pathways_TCGA".
#'
#' @param suffix Character string appended to the output PDF filename to
#'   distinguish between different runs or datasets. Default: "TCGA".
#'
#' @param taskHierarchy Data.frame containing pathway metadata. Must include:
#'   \itemize{
#'     \item Column named by `Depth` parameter (typically "gs_name"): pathway identifiers
#'     \item `class` column: metabolic class for each pathway
#'   }
#'   Default: `cfs` (expects global variable).
#'
#' @param nested Logical. When TRUE, use nested faceting to group pathways
#'   hierarchically by metabolic class. When FALSE, display all pathways
#'   without class-based grouping. Default: TRUE.
#'
#' @param top Logical. When TRUE, filter to show only pathways that are
#'   significantly dysregulated in more than the median number of cohorts.
#'   Useful for focusing on frequently dysregulated pathways. Default: FALSE.
#'
#' @param color_low Character string specifying the color for negative
#'   t-statistics (down-regulated pathways). Default: "skyblue".
#'
#' @param color_high Character string specifying the color for positive
#'   t-statistics (up-regulated pathways). Default: "purple".
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Filters pathways for significance using the specified threshold
#'   \item Classifies pathways as Up, Down, or NonSignificant based on logFC and threshold
#'   \item Maps pathway names using `taskHierarchy`
#'   \item Clips t-statistics to ±5 range to improve visualization
#'   \item Abbreviates metabolic class names for clarity (e.g., "CAR" for Carbohydrate,
#'     "LIP" for Lipid, "ENE" for Energy)
#'   \item Creates a point plot with:
#'     \itemize{
#'       \item X-axis: cancer cohorts
#'       \item Y-axis: pathways (wrapped text)
#'       \item Point size: inverse of significance (smaller p-value = larger point)
#'       \item Point color: t-statistic (blue for negative, purple for positive)
#'       \item Point shape: regulation direction (diamond for Up, circle for Down)
#'     }
#'   \item Optional nested faceting by metabolic class with free y-axis scales
#' }
#'
#' @importFrom dplyr %>% left_join mutate select filter pull
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes aes_string geom_point scale_shape_manual
#'   ggtitle ylab labs scale_color_gradient2 scale_size theme_bw scale_y_discrete guides
#'   guide_legend element_blank element_text unit theme
#' @importFrom stringr str_wrap str_to_title
#' @importFrom ggh4x facet_nested
#'
#' @return Invisibly returns the ggplot2 object representing the dot-plot.
#'   Primary side-effect is creation of a PDF file in `OUTDIRV` named
#'   "Dotplot_all_adj.P.ValCutoff<cutoff>_<suffix>.pdf"
#'
#' @examples
#' \dontrun{
#'   makeLimmaDotplot_TCGA(rlst_list, sig.cutoff = 0.05, OUTDIRV = "./plots/",
#'                         suffix = "TCGA")
#' }
#'
#' @export
#' 
makeLimmaDotplot_TCGA = function(rlst_list, sig_statistic="adj.P.Val", 
        sig.cutoff, adj.P.Val.cutoff=NULL, w=NULL, h=NULL, OUTDIRV, 
        Depth="gs_name", lowfigure=F, title = "Dysregulated pathways_TCGA", 
        suffix="TCGA", taskHierarchy=cfs, nested=T, top=FALSE,
        color_low="skyblue", color_high="purple"){
    if (!is.null(adj.P.Val.cutoff)) {
        sig.cutoff = adj.P.Val.cutoff
    }
    taskHierarchy$Task = taskHierarchy[[Depth]]
    rlst_long = lapply(seq_along(rlst_list), function(i) {
        rlst = rlst_list[[i]]
        rlst$Cancer = names(rlst_list)[i]
        rlst$Regulation = 
            ifelse(rlst$logFC>0&rlst[[sig_statistic]]<sig.cutoff, "Up", 
            ifelse (rlst$logFC<0&rlst[[sig_statistic]]<sig.cutoff, "Down", "NonSig"))
        rlst = rlst %>% tibble::rownames_to_column(var="Task")
        rlst[rlst[[sig_statistic]]<sig.cutoff, ]
    }) %>% Reduce(rbind, .) %>%
        dplyr::left_join(taskHierarchy) %>%
            dplyr::mutate(logFC=round(logFC, 2)) %>%
            dplyr::mutate(Class = gsub("Carbohydrate metabolism","CAR",class)) %>%
            dplyr::mutate(Class = gsub("Lipid metabolism", "LIP", Class)) %>%
            dplyr::mutate(Class = gsub("Metabolism of cofactors and vitamins", "COF", Class)) %>%
            dplyr::mutate(Class = gsub("Energy metabolism", "ENE", Class)) %>%
            dplyr::mutate(Class = gsub("Amino acid metabolism", "AMI", Class)) %>%
            dplyr::mutate(Class = gsub("Nucleotide metabolism", "NUC", Class)) %>%
            dplyr::mutate(Class = gsub("Biosynthesis of other secondary metabolites", "OSM", Class)) %>%
            dplyr::mutate(Class = gsub("Metabolism of other amino acids", "OAA", Class)) %>%
            dplyr::mutate(Class = gsub("Glycan biosynthesis and metabolism", "GLY", Class)) %>%
            dplyr::mutate(Class = gsub("Metabolism of terpenoids and polyketides", "TER", Class)) %>%
            dplyr::mutate(Class = gsub("Xenobiotics biodegradation and metabolism" , "XEN", Class))

    if (top) {
        temp = rlst_long %>% dplyr::select(Task, Cancer, Regulation) %>% 
            unique() %>% filter(Regulation %in% c("Up", "Down")) %>%
            pull(Task) %>% table()
        tops = names(temp)[temp > median(temp)]
        rlst_long = rlst_long %>% filter(Task %in% tops)
    }
    rlst_long$t[rlst_long$t>5] = 5
    rlst_long$t[rlst_long$t< -5] = -5
  
    p = ggplot2::ggplot(data=rlst_long, mapping=ggplot2::aes_string(x="Cancer", y="Task")) +
        ggplot2::geom_point(ggplot2::aes_string(size=sig_statistic, color="t", shape="Regulation"))
    if (nested) p = p + ggh4x::facet_nested(Class~., scales = "free_y", space = "free")
    p = p + ggplot2::scale_shape_manual(values = c("Up" = 18, "Down" = 20)) +
        ggplot2::ggtitle(title) +
        ggplot2::ylab("") +  ggplot2::labs(size=sig_statistic) +
        ggplot2::scale_color_gradient2(low=color_low, mid="snow2", high=color_high, midpoint=0) +
        # ggplot2::scale_color_gradient2(low="blue", mid="snow2", high="red", midpoint=0) +
        ggplot2::scale_size(trans="reverse") +
        ggplot2::theme_bw(base_size = 15) + 
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 1000)) +
        ggplot2::guides(alpha="none", size=ggplot2::guide_legend(override.aes=list(shape=21))) +
        ggplot2::theme(
            axis.title =  ggplot2::element_blank(),
            strip.text.y = ggplot2::element_text(size = 10, angle = 90, hjust=0.5),
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)) +
        ggplot2::labs(caption=paste0(sig_statistic,"Cutoff = ",sig.cutoff,", pAdjustMethod = BH"))

    if (lowfigure) {
      p = p + ggplot2::theme(
              axis.title = ggplot2::element_blank(),
              legend.position = "bottom",
              legend.box = "vertical",                  # Stack legends vertically
              legend.box.just = "left",                 # Center the whole legend box
              legend.text = ggplot2::element_text(size = 7),
              legend.title = ggplot2::element_text(size = 11),
              legend.spacing.y = ggplot2::unit(0, "cm") # Reduce spacing between legend rows
            ) +
            ggplot2::guides(
              fill  = ggplot2::guide_legend(ncol=1, title.position = "top", byrow = F),
              shape = ggplot2::guide_legend(title.position = "left"),
              size  = ggplot2::guide_legend(title.position = "left"))
    }
    if (is.null(h)|is.null(w)) {
        h = ((unique(rlst_long$gs_name) %>% length())/6) + 3
        w = (nchar(unique(rlst_long$gs_name)) %>% max())/11 + 3.5
    }
    pdf(paste0(OUTDIRV,"/Dotplot_all_", sig_statistic, "Cutoff",
        sig.cutoff,"_",suffix,".pdf"),w=w,h=h)
        show(p)
    dev.off()
}



 #' @title Visualize Differential Metabolic Effects by Pathway Class
 #'
 #' @description
 #' Generates violin plots showing the distribution of differential metabolic effects
 #' (represented by t-statistics or logFC) across KEGG metabolic pathway classes.
 #' Differential effects can optionally be weighted by pathway activity level
 #' (average expression) and winsorized to handle outliers. The function supports
 #' both single-level (cohort-based) and hierarchical (study-cohort) data structures.
 #'
 #' @param rlst_list A named list of results from \code{GSVAlimmaTest()} testing.
 #'   Supports two structures:
 #'   \itemize{
 #'     \item **Single-level**: List of cancer types, each containing
 #'       \code{$rslt_of_interest} data frame with differential expression results
 #'     \item **Hierarchical (depth=5)**: Nested list structure with studies
 #'       at top level, cancer types nested below, then \code{$rslt_of_interest}
 #'   }
 #'   Each result data frame must contain columns: \code{t} (or other
 #'   effect statistic), \code{logFC}, and significance column (default \code{adj.P.Val}).
 #'
 #' @param taskHierarchy A data.frame or tibble mapping pathways to metabolic classes,
 #'   with at least two columns:
 #'   \itemize{
 #'     \item `gs_name`: pathway identifier (must match row names in results)
 #'     \item `class`: metabolic class category for grouping
 #'   }
 #'
 #' @param OUT_DIR Character string specifying the output directory where the PDF
 #'   will be saved. Must be a valid, writable directory path.
 #'
 #' @param weight.effect.by.gene Logical. When TRUE, weights the differential effect
 #'   statistic by the inverse rank of pathway activity level (AveExpr), giving
 #'   more weight to effects in highly active pathways. Default: TRUE.
 #'
 #' @param effect.statistic Character string naming the column in result data frames
 #'   to use as the differential effect measure (e.g., "t" for t-statistic).
 #'   Default: "t".
 #'
 #' @param w Numeric. Width of output PDF in inches. Default: 2.
 #'
 #' @param h Numeric. Height of output PDF in inches. Default: 3.5.
 #'
 #' @param title Character string for the plot title (typically dataset or study name).
 #'   Default: "Zhou2020".
 #'
 #' @param winsorize Logical. When TRUE, clips extreme values to the specified
 #'   probability quantiles to reduce the influence of outliers. Default: TRUE.
 #'
 #' @param probs Numeric. Probability level for winsorization symmetrically applied
 #'   at both tails (e.g., 0.05 clips at 5th and 95th percentiles).
 #'   Ignored if \code{winsorize=FALSE}. Default: 0.05.
 #'
 #' @param suffix Character string appended to the output PDF filename for
 #'   distinguishing between different analyses. Default: "".
 #'
 #' @param sig.statistic Character string naming the significance column in result
 #'   data frames (e.g., "adj.P.Val", "p.value"). Default: "adj.P.Val".
 #'
 #' @param sig.cutoff Numeric. Significance threshold for filtering pathways.
 #'   Only pathways with `sig.statistic < sig.cutoff` are included.
 #'   Default: 0.2.
 #'
 #' @details
 #' The function performs the following steps:
 #' \enumerate{
 #'   \item Detects data structure depth to handle single-level or nested hierarchies
 #'   \item Filters pathways for significance using the specified threshold
 #'   \item Maps pathway results to metabolic classes
 #'   \item Standardizes class names by removing common KEGG descriptors
 #'   \item Optionally weights effects by inverse pathway activity rank
 #'   \item Optionally winsorizes extreme weighted effects
 #'   \item Creates violin plots with:
 #'     \itemize{
 #'       \item X-axis: metabolic classes (ordered by median effect)
 #'       \item Y-axis: weighted differential effects
 #'       \item Fill color: Source/Study (for hierarchical data)
 #'     }
 #' }
 #' Violin plots effectively show the full distribution of effects within
 #' each pathway class, revealing both typical and extreme dysregulation patterns.
 #'
 #' @importFrom dplyr %>% left_join mutate group_by summarise arrange pull
 #' @importFrom tibble rownames_to_column
 #' @importFrom ggplot2 ggplot aes geom_violin position_dodge geom_boxplot labs
 #'   theme_minimal theme element_text
 #' @importFrom stringr str_to_title
 #'
 #' @return Invisibly returns NULL. A PDF file named "diffEffectBoxplot_system<suffix>.pdf"
 #'   is written to `OUT_DIR` containing the violin plot visualization.
 #'
 #' @export
 #' 
 diffEffectBoxplot_bySystem_GSVA = function(
        rlst_list, taskHierarchy, OUT_DIR, 
        weight.effect.by.gene=F, effect.statistic="t",
        w = 2, h = 3.5, title = "Zhou2020",
        winsorize=TRUE, probs=0.05, suffix="",
        sig.statistic="adj.P.Val", sig.cutoff=0.2) {
    depth = function(x) {
        if (!is.list(x)) return(0)
        1 + max(sapply(x, depth))
        }
    if (depth(rlst_list)==5) {
        B_Z_N0 = lapply(names(rlst_list), function(study){
            lapply(names(rlst_list[[study]]), function(c) {
                temp = rlst_list[[study]][[c]]$rslt_of_interest
                temp = temp[temp[[sig.statistic]] < sig.cutoff, ]
                temp$Cancer = c
                temp = as.data.frame(temp) %>% tibble::rownames_to_column(var="gs_name")
                temp$Source = study
                temp
            }) %>% Reduce(rbind, .) %>% as.data.frame()
        }) %>% Reduce(rbind, .) %>% as.data.frame() 
    } else {
        B_Z_N0 = lapply(names(rlst_list), function(cancer){
            temp = rlst_list[[cancer]]$rslt_of_interest
            temp = temp[temp[[sig.statistic]] < sig.cutoff, ]
            if (nrow(temp)>0) {
                temp$Cancer = cancer
                return(temp %>% tibble::rownames_to_column(var="gs_name"))
            } else {
                return(NA)
            }
        })
        B_Z_N0 = B_Z_N0[!is.na(B_Z_N0)] %>% Reduce(rbind, .) %>% as.data.frame()
    }
    B_Z_N0 = B_Z_N0 %>%
        dplyr::left_join(taskHierarchy %>% unique()) %>%
        dplyr::mutate(class = gsub(" metabolism", "", class)) %>%
        dplyr::mutate(class = gsub("Metabolism of ", "", class)) %>%
        dplyr::mutate(class= gsub(" biosynthesis and", "", class)) %>%
        dplyr::mutate(class = gsub(" biosynthesis", "", class)) %>%
        dplyr::mutate(class = gsub(" biodegradation and", "", class)) %>%
        dplyr::mutate(class = gsub("Biosynthesis of ", "", class)) %>%
        dplyr::mutate(class = stringr::str_to_title(class))
    B_Z_N0$effect = B_Z_N0[[effect.statistic]]
    if (weight.effect.by.gene) {
        # B_Z_N0$AveExpr_scaled = 
        #     (B_Z_N0$AveExpr-min(B_Z_N0$AveExpr))/(max(B_Z_N0$AveExpr)-min(B_Z_N0$AveExpr))
        B_Z_N0$AveExpr_scaled = rev(rank(B_Z_N0$AveExpr, ties.method = "average"))
        B_Z_N0$effect_weighted = (B_Z_N0$effect)/(B_Z_N0$AveExpr_scaled)
    } else {
        B_Z_N0$effect_weighted = B_Z_N0$effect
    }
    if (winsorize){
        B_Z_N0$effect_weighted[B_Z_N0$effect_weighted>quantile(B_Z_N0$effect_weighted, probs=(1-probs))] = 
            quantile(B_Z_N0$effect_weighted, probs=(1-probs))
        B_Z_N0$effect_weighted[B_Z_N0$effect_weighted<quantile(B_Z_N0$effect_weighted, probs=probs)] = 
            quantile(B_Z_N0$effect_weighted, probs=probs)   
    }

    System_orderZ = B_Z_N0 %>% unique() %>% group_by(class) %>%
        # dplyr::summarise(meanEffect = median(effect_weighted, na.rm=T), .groups = "drop") %>% 
        dplyr::summarise(meanEffect = median(effect_weighted, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(class) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(class=factor(class, levels = System_orderZ))

    if (depth(rlst_list)==5) {
        B_Z_N = B_Z_N %>% 
            ggplot2::ggplot(., ggplot2::aes(x=class, y=effect_weighted, fill=Source)) +
            ggplot2::geom_violin(position = position_dodge(width = 0.85),  # reduce spacing
                width=2, trim=FALSE, alpha=0.75, col=NA)+
            scale_fill_manual(values = 
                c("#E69F00", "#56B4E9","#999999","#F0E442","#009E73","#0072B2","#D55E00","#CC79A7"))
    } else {
        B_Z_N = B_Z_N %>%  
            ggplot2::ggplot(., ggplot2::aes(x=class, y=effect_weighted, fill=class)) +
            ggplot2::geom_violin(trim=FALSE, color=NA, alpha=0.6, width=1.1) #fill="#E69F00", 
    }
        # ggplot2::ggplot(., ggplot2::aes(x=class, y=effect_weighted)) +
        # ggplot2::geom_boxplot(fill="green4") +
     B_Z_N = B_Z_N + ggplot2::labs(title = title,
            x="Class", y=effect.statistic)+
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="right",
            axis.text.x = ggplot2::element_text(angle = 135, hjust=1, vjust=1))

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_system",suffix,".pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
 }



 #' Differential effect boxplot by cancer/group
 #' @title Visualize Differential Metabolic Effects by Cancer Type
 #'
 #' @description
 #' Generates boxplots showing the distribution of differential metabolic effects
 #' (t-statistics) across different cancer types or cohorts. Effects are weighted
 #' by pathway activity level (inverse rank of average expression) and winsorized
 #' at the 5th and 95th percentiles to mitigate outlier influence. Useful for
 #' comparing metabolic dysregulation patterns across multiple cancer cohorts.
 #'
 #' @param rlst_list A named list of results from \code{GSVAlimmaTest()} 
 #'   (one element per cancer type or cohort). Each element must contain
 #'   a component \code{$rslt_of_interest}, which should be a data frame produced by
 #'   \code{limma::topTable()} with at least these columns:
 #'   \itemize{
 #'     \item `t`: t-statistic from limma testing
 #'     \item `logFC`: log2 fold-change values
 #'     \item `AveExpr`: average expression level
 #'     \item `adj.P.Val`: adjusted p-values (optional, for filtering)
 #'   }
 #'   Row names should be pathway identifiers. List names are used as cancer type labels.
 #'
 #' @param OUT_DIR Character string specifying the output directory where the PDF
 #'   will be saved. Must be a valid, writable directory path.
 #'
 #' @param w Numeric. Width of output PDF in inches. Default: 2.75.
 #'
 #' @param h Numeric. Height of output PDF in inches. Default: 3.0.
 #'
 #' @param title Character string for the plot title (typically dataset or cohort name).
 #'   Default: "Zhou2020".
 #'
 #' @details
 #' The function:
 #' \enumerate{
 #'   \item Collects differential effect results across all cancer types
 #'   \item Computes inverse rank of average expression for each pathway
 #'   \item Weights t-statistics by the inverse expression rank to emphasize
 #'     effects in highly active pathways
 #'   \item Winsorizes weighted effects at 5th and 95th percentiles to
 #'     reduce extreme value influence
 #'   \item Orders cancer types by median t-statistic
 #'   \item Creates boxplots with:
 #'     \itemize{
 #'       \item X-axis: cancer types (ordered by median effect)
 #'       \item Y-axis: weighted differential effects
 #'       \item Green colored boxes showing quartile distributions
 #'     }
 #' }
 #' This visualization allows rapid identification of cancer types with
 #' predominantly up-regulated vs. down-regulated metabolic pathways.
 #'
 #' @importFrom dplyr %>% mutate group_by summarise arrange pull
 #' @importFrom tibble rownames_to_column
 #' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal theme element_text
 #'
 #' @return Invisibly returns NULL. A PDF file named "diffEffectBoxplot_cancer.pdf"
 #'   is written to `OUT_DIR` containing the boxplot visualization.
 #'
 #' @export
 #' 
diffEffectBoxplot_byCancer_GSVA = function(
    rlst_list, OUT_DIR, w=2.75, h=3.0, title="Zhou2020") {
    B_Z_N0 = lapply(names(rlst_list), function(cancer){
        temp = rlst_list[[cancer]]$rslt_of_interest
        temp$Cancer_type = cancer
        temp %>% tibble::rownames_to_column(var="gs_name") 
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>% unique()

    # B_Z_N0$AveExpr_scaled = 
    #     (B_Z_N0$AveExpr-min(B_Z_N0$AveExpr))/(max(B_Z_N0$AveExpr)-min(B_Z_N0$AveExpr))
    B_Z_N0$AveExpr_scaled = rev(rank(B_Z_N0$AveExpr, ties.method = "average")) 
    B_Z_N0$logFC_weighted = (B_Z_N0$logFC)/(B_Z_N0$AveExpr_scaled)
    B_Z_N0$logFC_weighted[B_Z_N0$logFC_weighted>quantile(B_Z_N0$logFC_weighted, probs=0.95)] = 
        quantile(B_Z_N0$logFC_weighted, probs=0.95)
    B_Z_N0$logFC_weighted[B_Z_N0$logFC_weighted<quantile(B_Z_N0$logFC_weighted, probs=0.05)] = 
        quantile(B_Z_N0$logFC_weighted, probs=0.05)   

    Cancer_orderZ = B_Z_N0 %>% unique() %>% group_by(Cancer_type) %>%
        dplyr::summarise(meanEffect = median(t, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(Cancer_type) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(Cancer=factor(Cancer_type, levels = Cancer_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=Cancer, y=t)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="Cancer", y="Wei. Diff. Effect")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_cancer.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}







# Functions for parsing cellFie metabolic task scores

## How to run CellFie
    # -   Normalize the data to make them comparable across samples
    # -   Sumbit the normalized data file to the [genePatten webserver](https://cloud.genepattern.org/gp/pages/index.jsf), and run it there with the following settings:
        # -   Data="PsoPproteomics.csv",
        # -   Sample.Number=,
        # -   Reference.Model= Human_iHsa,
        # -   Thresholding.Approach="local",
        # -   Percentile.Or.Value="percentile",
        # -   Local.Threshold.Type="minmaxmean",
        # -   Local.Lower.Bound="25",
        # -   Local.Upper.Bound="75"
        # -   Uses gene-specific thresholds, customized per gene (based on its expression distribution)
        # -   For each gene, CellFie computes its 25th and 75th percentile and considers it ON if expression > mean of [25th, 75th].


 #' Sum up task scores to higher level categories
 #'
 #' Aggregate task-level CellFie scores to a higher hierarchical depth
 #' (for example, \code{Depth1} or \code{Depth2}) and produce a per-sample
 #' summary table plus a PDF visualization. The summary contains mean and sum
 #' statistics for raw and active task scores as well as the number and
 #' fraction of active tasks per sample.
 #'
 #' @param cfs data.frame or tibble. A CellFie output object (as returned by
 #'   \code{processCellFieOutput()}) containing at least the columns
 #'   \code{Sample}, the requested depth column (e.g. \code{Depth1} or
 #'   \code{Depth2}), \code{TaskScore}, \code{BinaryTaskScore} and
 #'   \code{ActiveTaskScore}.
 #' @param Depth character(1). Name of the hierarchical depth to aggregate
 #'   (e.g. \code{"Depth1"} or \code{"Depth2"}).
 #' @param h numeric(1). Height of the output PDF in inches (default: 12).
 #' @param w numeric(1). Width of the output PDF in inches (default: 30).
 #' @param OUTDIR character(1). Directory where the output PDF
 #'   (\file{BinaryTaskScores_<Depth>.pdf}) will be written. Directory must be
 #'   writable; caller should create it if it does not exist.
 #' 
 #' @importFrom dplyr %>%
 #'
 #' @return A data.frame with one row per sample and category containing the
 #'   aggregated statistics: \code{MeanTaskScore}, \code{MeanActiveTaskScore},
 #'   \code{SumTaskScore}, \code{SumActiveTaskScore}, \code{FractionOfActiveTasks}
 #'   and \code{NrOfActiveTasks}. Invisibly returns the same data.frame; the
 #'   primary side-effect is writing the PDF to \code{OUTDIR}.
 #'
 #' @examples
 #' # depths <- sumUpTaskScores(cfs, Depth = "Depth1", h = 12, w = 30, OUTDIR = "./out")
 #'
 #' @export
 #' 
 sumUpTaskScores = function(cfs, Depth, h=12, w=30, OUTDIR) {
    cols_to_group = c("Sample", Depth)
    Depths = cfs %>%
        group_by(across(all_of(cols_to_group))) %>%
        summarise(
            MeanTaskScore = mean(TaskScore, na.rm = TRUE),
            MeanActiveTaskScore = mean(ActiveTaskScore, na.rm = TRUE),
            SumTaskScore = sum(TaskScore, na.rm = TRUE),
            SumActiveTaskScore = sum(ActiveTaskScore, na.rm = TRUE),
            FractionOfActiveTasks = mean(BinaryTaskScore==1, na.rm = TRUE),
            NrOfActiveTasks = sum(BinaryTaskScore==1, na.rm = TRUE)
        ) %>% dplyr::left_join(
            cfs[, c("Sample", "Condition", "Cancer_Condition", "Cancer_type")] %>% unique()
        )
    L = Depths$Cancer_Condition %>% unique() %>% length()
    Depths$Condition = Depths$Condition %>% gsub("NAT","Normal",.)
    P = ggplot2::ggplot(Depths, ggplot2::aes_string(y=Depth, x="Sample")) +
        ggplot2::geom_point(shape=18, ggplot2::aes(size=abs(FractionOfActiveTasks), color=Condition)) +
        ggplot2::labs(size = "FractionOfActiveTasks") + 
        ggplot2::scale_color_manual(values=c("Tumor"="darkorange", "Normal"="darkturquoise")) +
        ggplot2::xlab("") + 
        ggplot2::ggtitle(paste0("Fraction of active tasks_", Depth)) + 
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          strip.text.x = ggplot2::element_text(size = 14),
          axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))
    pdf(paste0(OUTDIR, "/BinaryTaskScores_", Depth,".pdf"),h=h,w=w)
        show(P)
    dev.off()
    return(Depths)
}

 #' Generate a boxplot for visualizing the number of genes in each task.
 #'
 #' Create a boxplot that summarizes the distribution of gene counts per task
 #' (grouped by higher-level system). This is useful to compare task sizes
 #' across systems and samples. The function writes a PDF to the provided
 #' output directory and returns invisibly the summarized data frame.
 #'
 #' @param taskInfo data.frame or tibble. Table with at least the columns
 #'   \code{Sample}, \code{Depth3}, \code{Depth1} and
 #'   \code{GeneAssociatedToEssentialRxnsTask}. Rows should represent task->gene
 #'   assignments (one row per gene associated with a task and sample).
 #' @param OUTDIR character(1). Directory where the output PDF
 #'   (\file{Task_size_system.pdf}) will be written. The directory must be
 #'   writable; the caller should create it if necessary.
 #' 
 #' @importFrom dplyr %>%
 #'
 #' @return Invisibly returns the summarized data frame of gene counts per
 #'   task and system. The main side-effect is writing
 #'   \file{Task_size_system.pdf} to \code{OUTDIR}.
 #'
 #' @examples
 #' # plotTaskSizes(taskInfo, OUTDIR = "./out")
 #'
 #' @export
 #' 
 plotTaskSizes = function(taskInfo, OUTDIR) {
    taskInfo = taskInfo %>% unique() %>% 
        dplyr::select(Sample, Depth3, Depth1, GeneAssociatedToEssentialRxnsTask) %>%
        unique()
    geneCounts = taskInfo %>% group_by(Sample, Depth3) %>%
        summarise (geneCount = n())
    taskInfo = taskInfo %>% dplyr::select(Depth3, Depth1) %>% unique()
    geneCounts = geneCounts[geneCounts$Depth3!="", ] %>% left_join(
        taskInfo, by = "Depth3")
    rm(taskInfo)
    gc()
    P = geneCounts %>% dplyr::mutate(System = gsub("METABOLISM","M.", Depth1)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=System, y=geneCount)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(# title = "TCGA",
                x="", y="Task size")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))
    pdf(paste0(OUTDIR, "/Task_size_system.pdf"), w=2.5, h=3.5)
        print(P)
    dev.off()
}



 #' Circular packing plot of CellFie task hierarchy
 #'
 #' Build and save a circular packing plot that summarizes tasks and their
 #' grouping by higher-level systems using CellFie output. Labels are
 #' shortened for display; the plot is written as a PDF to the provided
 #' output directory.
 #'
 #' @param cfs data.frame or tibble. A CellFie output object (as returned by
 #'   \code{processCellFieOutput()}) containing at minimum the columns
 #'   \code{Sample}, the hierarchical depth columns (e.g. \code{Depth1},
 #'   \code{Depth2}, \code{Depth3}), and any task identifiers needed to
 #'   reconstruct the hierarchy. Rows should represent task->gene or task
 #'   membership entries.
 #' @param out_dir character(1). Directory where the output PDF
 #'   (\file{TaskSummary.pdf}) will be written. The directory must be
 #'   writable; caller should create it if needed.
 #' 
 #' @importFrom dplyr %>%
 #'
 #' @return Invisibly returns the ggplot2 object used to render the figure.
 #'   The primary side-effect is writing \file{TaskSummary.pdf} to
 #'   \code{out_dir}.
 #'
 #' @export
 #' 
 getCirPackingPlot_CellFie = function(cfs, out_dir) {
    edges = cfs%>%dplyr::mutate(Depth0="Metabolism")%>%.[,c("Depth0","Depth1")]%>%
        unique()%>%setNames(c("from", "to")) %>%
        rbind(cfs[, c("Depth1","Depth2")]%>%unique()%>%setNames(c("from", "to"))) %>%
        rbind(cfs[, c("Depth2","Depth3")]%>%unique()%>%setNames(c("from","to")))
    vertices = cfs%>%dplyr::mutate(Depth0="Metabolism")%>%.[,c("Depth0","Depth3")]%>%unique()%>%
        pull(Depth0)%>%table()%>%as.data.frame()%>%setNames(c("name", "size")) %>%
        rbind(cfs[,c("Depth1","Depth3")]%>%unique()%>%pull(Depth1)%>%table()%>%
        as.data.frame()%>%setNames(c("name", "size"))) %>%
        rbind(cfs[,c("Depth2","Depth3")]%>%unique()%>%pull(Depth2)%>%table()%>%
        as.data.frame()%>%setNames(c("name", "size"))) %>%
        rbind(cfs%>%pull(Depth3)%>%unique()%>%table()%>%as.data.frame()%>%
        setNames(c("name","size"))) %>%
        dplyr::mutate(size = size) %>%
        dplyr::mutate(shortName=gsub("METABOLISM","M.",paste0(name,"(",size,")")))
    vertices$shortName[!vertices$name%in%unique(cfs$Depth1)] = NA
    mygraph = igraph::graph_from_data_frame(edges, vertices=vertices)       
    p1 = ggraph::ggraph(mygraph, layout = 'circlepack', weight=size) +
        ggraph::geom_node_circle(ggplot2::aes(col=depth)) +
        ggraph::geom_node_label(ggplot2::aes(label=shortName),repel=T,size=4,color="darkgreen") +
        ggplot2::theme_void() +
        viridis::scale_color_viridis() +
        ggplot2::theme(legend.position="FALSE") +
        ggplot2::ggtitle("")

    pdf(paste0(out_dir,"/TaskSummary.pdf"), h=5, w=4.85)
        print(p1)
    dev.off()
}



 #' Process the CellFie output
 #'
 #' Import, clean and summarise CellFie output directories. This function
 #' parses per-cancer (or per-group) CellFie outputs, reshapes task-level
 #' results so tasks are rows and samples are columns, merges sample
 #' metadata, computes aggregated summaries at multiple hierarchical depths
 #' (Depth1, Depth2, Depth3), and writes summary tables and diagnostic
 #' figures (boxplots and a circular packing plot) to the provided output
 #' location.
 #'
 #' @param combined logical(1). If TRUE, treat all samples (including different cancer types) as coming from a
 #'   single combined run; otherwise, CellFie output is expected in
 #'   per-group subfolders and group names are inferred from these folder
 #'   names (default: FALSE).
 #' @param outdir character(1). Path to the directory containing CellFie
 #'   output folders. Each group's output should be in its own subdirectory
 #'   under this path. The function will write summaries under
 #'   file.path(outdir, "../Summary/").
 #' @param SampleNames character(). Vector of sample names in the same order
 #'   as used when running CellFie. 
 #' @param meta data.frame. Sample metadata that must contain at least the
 #'   columns \code{Sample} and \code{Cancer_type}. When \code{combined=FALSE}
 #'   the function uses \code{meta$Cancer_type} to match samples to group
 #'   folders. Additional columns used include \code{Condition} and
 #'   \code{Cancer_Condition}.
 #' @param samples2keep character(). Vector of sample IDs to retain; useful
 #'   to exclude unwanted samples after import (default: all samples present).
 #' 
 #' @importFrom dplyr %>%
 #'
 #' @return A single data.frame (invisibly) containing the concatenated,
 #'   post-processed CellFie results across groups. Side-effects: multiple
 #'   summary files and PDF figures are written to \code{file.path(outdir, "../Summary/")},
 #'   including \file{TaskHierarchy_summary.xlsx}, \file{cfs.RDS}, task score
 #'   tables for each depth, and several PDF plots.
 #'
 #' @examples
 #' # processCellFieOutput(combined = FALSE, outdir = "./CellFieOut/", SampleNames = samples, meta = meta, samples2keep = samples)
 #'
 #' @export
 #' 
 processCellFieOutput = function(combined=FALSE,
     outdir = paste0(OUTDIR, "/CellFieOut/"), SampleNames, meta, samples2keep) {
    # list folders 
    cancers = list.dirs(path = outdir, full.names = FALSE, recursive = FALSE) %>%
        gsub("_.*", "", .)

    # Import files
    files = list.files(outdir, 
        pattern="detailScoring.csv", full.names=TRUE, recursive=T)
    cfs = lapply(files, function(f) {
        read.csv2(f, 
            sep=",", stringsAsFactors = FALSE) 
        }) %>% setNames(cancers)
    
    # Convert complex number to only keep the real part
    cfs = lapply(cfs, function(cf) {
        complexCol = grep("TaskScore_|ExpressionScoreEssentialRxnsTask_", colnames(cf))
        for (i in complexCol) {
            cf[, i] = gsub("\\+\\-", "-", cf[, i])
            cf[, i] = Re(as.complex(cf[, i]))
        }
        cf
    })

    # Summarize task nr
    taskInfo = read.csv2(gsub("detailScoring.csv", "taskInfo.csv", files[1]), 
        sep=",", stringsAsFactors = FALSE) %>%
        setNames(c("Depth3.ID", "Depth3", "Depth1", "Depth2"))
    # Here I would like move some core glycan metabolic pathways from carbohydrate metabolism to glycan metabolism
    tasks2mv = c("Synthesis of lactose", "UDP-glucose synthesis", 
        "UDP-galactose synthesis", "UDP-glucuronate synthesis", "GDP-L-fucose synthesis", "GDP-mannose synthesis",
        "UDP-N-acetyl D-galactosamine synthesis", "CMP-N-acetylneuraminate synthesis", "N-Acetylglucosamine synthesis",
        "Glucuronate synthesis (via inositol)", "Glucuronate synthesis (via udp-glucose)")
    taskInfo$Depth1[taskInfo$Depth3 %in% tasks2mv] = "GLYCAN METABOLISM"
    TaskNr1 = taskInfo %>% group_by(Depth1) %>%
        summarise(Depth2 = n_distinct(Depth2))
    TaskNr2 = taskInfo %>% group_by(Depth1) %>%
        summarise(Depth3 = n_distinct(Depth3))
    TaskNr = dplyr::left_join(TaskNr1, TaskNr2)
    TableNames = list(TaskNr) %>% setNames("TaskHierarchy_summary")
    paste0(outdir,"/../Summary/") %>% dir.create(., recursive=T, showWarnings=F)
    openxlsx::write.xlsx(TableNames, 
        file=paste0(outdir,"/../Summary/TaskHierarchy_summary.xlsx"))
    rm(TaskNr1, TaskNr2, TaskNr, TableNames)
    gc()

    # Reshuffle the table to have tasks in the rows and samples in the columns
    cfs = lapply(cfs, function(cf) {
        Samples = gsub(".*_", "", colnames(cf)) %>% unique()
        cf2 =lapply(Samples, function(S) {
            temp = cf[, grep(paste0(S, "$"), colnames(cf))]
            colnames(temp) = gsub("_.*", "", colnames(temp))
            temp
        }) %>% Reduce(rbind, .) %>% 
            dplyr::mutate(Depth3.ID = TaskID) %>%
            dplyr::select(Depth3.ID, everything()) %>% 
            filter(TaskScore!=-1) # remove unscorable rows (-1)
        cf2
    })

    ## Add in task depth information
    cfs = lapply(cfs, function(cf) {
        dplyr::right_join(taskInfo, cf)
        })
    rm(taskInfo)
    gc()

    # Add in Sample information
    cfs = lapply(cancers, function(c){
        if (combined) {
            samples = SampleNames
        } else {
            samples = filter(meta, Cancer_type==c) %>% pull(Sample) %>% intersect(SampleNames,.)
        }
        SampleInfo = data.frame(Sample = samples,
            SampleID = 1:length(samples)
        )
        # SampleInfo = data.frame(Sample = grep(c, SampleNames, value=T),
        #             SampleID = 1:length(grep(c, SampleNames, value=T)))
        cf = cfs[[c]]
        cf4 = dplyr::right_join(SampleInfo, cf) %>%
            dplyr::right_join(meta[,c("Sample", "Cancer_Condition", "Condition", "Cancer_type")], .) %>%
            # dplyr::mutate(Group=gsub("_.*_", "_", Sample)) %>% ##########
            dplyr::select(Condition, everything()) 
        write.table(cf4, paste0(outdir, "/",c,"/detailScoring_new.csv"), sep="\t", 
            quote=F, row.names=F)
        cf5 = dplyr::select(cf4, -EssentialRxnsTask, -ExpressionScoreEssentialRxnsTask,
            -GeneAssociatedToEssentialRxnsTask, -GeneExpressionValue) %>%
            unique()
        cf5
    }) %>% setNames(cancers) %>% Reduce(rbind, .)
    cfs = filter(cfs, Sample%in%samples2keep)
    cfs$ActiveTaskScore = cfs$TaskScore*cfs$BinaryTaskScore
    saveRDS(cfs, paste0(outdir, "/../Summary/cfs.RDS"))

    # Sum up task scores to higher level categories
    ## Task depth1
    Depth = "Depth1"
    w = (cfs$Sample %>% unique() %>% length())/6.5 + 5
    h = ((cfs[[Depth]] %>% unique() %>% length())/15) * 3.5 + 3
    Depth1 = sumUpTaskScores(cfs, Depth=Depth, h=h, w=w,
            OUTDIR=paste0(outdir, "/../Summary/")) %>% 
            dplyr::mutate(TaskScore=MeanTaskScore)
    saveRDS(Depth1, paste0(outdir, "/../Summary/TaskScores_", Depth,".RDS"))
    ## Task depth2
    Depth = "Depth2"
    w = (cfs$Sample %>% unique() %>% length())/6.5 + 5
    h = ((cfs[[Depth]] %>% unique() %>% length())/15) * 3.5 + 2
    Depth2 = sumUpTaskScores(cfs, Depth=Depth, h=h, w=w,
            OUTDIR=paste0(outdir, "/../Summary/"))%>%
            dplyr::mutate(TaskScore=MeanTaskScore)
    saveRDS(Depth2, paste0(outdir, "/../Summary/TaskScores_", Depth,".RDS"))
    ## Task depth3
    Depth = "Depth3"
    L = cfs$Cancer_Condition %>% unique() %>% length() ###############
    # cfs$Condition = gsub(".*_", "", cfs$Group) %>% gsub("NAT","Normal",.) #################
    cfs$Condition = gsub("NAT","Normal", cfs$Condition)
    w = (cfs$Sample %>% unique() %>% length())/5.5 + 7
    h = ((cfs[[Depth]] %>% unique() %>% length())/15) * 3
    pdf(paste0(outdir, "/../Summary/BinaryTaskScores_", Depth,".pdf"),h=h,w=w)
        print(ggplot2::ggplot(cfs, ggplot2::aes_string(y=Depth, x="Sample")) +
        ggplot2::geom_point(shape=18, aes(size=BinaryTaskScore, color=Condition)) +
        ggplot2::labs(size = "BinaryTaskScore") + 
        ggplot2::scale_color_manual(values=c("Tumor"="darkorange", "Normal"="darkturquoise")) +
        ggplot2::xlab("") + 
        ggplot2::ggtitle(paste0("Binary Task Score_", Depth)) + 
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme( # legend.position = "none",
        strip.text.x = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()
    openxlsx::write.xlsx(
        list(
            Depth3 = cfs,
            Depth2 = Depth2,
            Depth1 = Depth1), 
        file=paste0(outdir,"/../Summary/TaskScores_differentDepths.xlsx")
    )

    # Boxplot of task score means at depth2
    w = (cfs$Cancer_Condition %>% unique() %>% length())/3.1*5
    pdf(paste0(outdir,"/../Summary/Boxplots_TaskScoreMeans_depth2.pdf"), w=w, h=24)
        print(ggplot2::ggplot(Depth2, ggplot2::aes(x=Cancer_Condition, y=MeanTaskScore, fill=Cancer_type)) +
        ggplot2::geom_boxplot(position=position_dodge()) +
        ggplot2::facet_wrap(~Depth2, ncol=10, scales="free") +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 14),
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()

    # Boxplot of task score means at depth1
    w = ((cfs$Cancer_Condition %>% unique() %>% length())/3.1*5) * 0.4
    pdf(paste0(outdir,"/../Summary/Boxplots_TaskScoreMeans_depth1.pdf"), w=w, h=8)
        print(ggplot2::ggplot(Depth1, ggplot2::aes(x=Cancer_Condition, y=MeanTaskScore, fill=Cancer_type)) +
        ggplot2::geom_boxplot(position=position_dodge()) +
        ggplot2::facet_wrap(~Depth1, ncol=4, scales="free") +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 14),
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()
    rm(Depth2, Depth1)
    gc()

    # Plot task hierarchy
    getCirPackingPlot_CellFie(cfs, out_dir=paste0(outdir,"/../Summary/"))

    return(cfs)
}




#' Perform limma differential analysis on CellFie-inferred metabolic activities
#'
#' Fit linear models and perform empirical Bayes moderation using limma on
#' CellFie task scores. The function supports unpaired and paired designs and
#' can include an optional covariate alongside the primary `Condition` factor.
#'
#' @param dat numeric matrix or data.frame of CellFie scores (log2 transformed) with features as
#'   rownames and samples as column names. Rows are tested features (e.g. tasks),
#'   columns are samples.
#' @param paired logical scalar. If TRUE a paired analysis is fitted using the
#'   blocking column specified by `paired_variable` (requires at least two
#'   paired blocks). Default: TRUE.
#' @param currentCovariate character scalar or NULL. Optional covariate name
#'   (a column in `meta`) to include in the model in addition to
#'   `Condition` (for example "Batch" or "Age"). If NULL only `Condition`
#'   is used. Default: NULL.
#' @param checkCovariate logical scalar. If TRUE perform basic checks on
#'   `currentCovariate` (existence in `meta` and not constant). Default: FALSE.
#'   if checkCovariate=TRUE, currentCovariate must be one a single covariate, 
#'   this is for testing each covariate individually to decide whether to keep them in the model or not.
#' @param meta data.frame with at least the columns `Sample` and `Condition`.
#'   `Sample` must match the column names of `dat`. Currently `Condition` is
#'   expected to have two categories.
#' @param paired_variable character scalar giving the column name in `meta`
#'   that links paired samples (for example a patient identifier). Default:
#'   "Patient".
#' 
#' @importFrom dplyr %>%
#' 
#' @return A list with components:
#'   - `rslt_of_interest`: Differential result table for Group.
#'   - `rslt_co`: Differential result table for covariant.
#'   - `data`: numeric matrix or data.frame of CellFie scores that were used in limma test.
#'   - `fitted.model`: fitted limma model.
#'
#' @details The function first filters out features lacking sufficient
#'   non-missing observations across groups. If `paired = TRUE` it uses
#'   `limma::duplicateCorrelation` and fits a block model with a consensus
#'   correlation estimate. Column names in the returned design matrix are
#'   cleaned to remove a `Condition` prefix for readability.
#'
#' @examples
#' # dat: matrix of CellFie scores; meta: data.frame with Sample and Condition
#' # res <- limmaTest_CellFie(dat, paired = TRUE, meta = meta,
#' #                         paired_variable = "Patient")
#' @export
#' 
limmaTest_CellFie = function(dat, paired=TRUE, currentCovariate=NULL, 
    checkCovariate=FALSE, meta=meta, paired_variable="Patient") {
    meta = filter(meta, Sample%in%colnames(dat))
    # Make sure the sample are in the same order in both meta and dat:
    dat = dat[, meta$Sample]
    meta = meta %>% dplyr::mutate(across(where(is.factor), droplevels))
    # print("When running differential analysis, we only keep features that have at least 2 values in each group...")
    featureskeep = apply(!is.na(dat), 1, function(row) {
        A=row
        B=meta$Condition
        sum(aggregate(A~B, data=data.frame(A=A,B=B), sum)[,"A"]>1)>1
    })
    dat = dat[featureskeep, ]
    #print(dim(dat))
    #if (!is.null(currentCovariate)) {currentCovariate=NULL}
    baseFormula = ~ Condition
    if (is.null(currentCovariate)) {
      currentFormula = baseFormula
    } else {
      currentFormula = as.formula(paste0("~ ", currentCovariate, " + Condition"))}
    
    design = model.matrix(currentFormula, data=meta) # each column is for coefficient
    group.categories = meta$Condition %>% unique() %>% length()
    colnames(design) = gsub("Condition", "", colnames(design))
    
    if (paired) {
      #Fit with correlated arrays
      dupcor = limma::duplicateCorrelation(dat, design, 
        block=meta[,paired_variable,drop=T])
      fit = limma::lmFit(dat, design, block=meta[,paired_variable,drop=T],
        correlation=dupcor$consensus)
    } else {
      fit = limma::lmFit(dat, design)
    }
    #fit3 = eBayes(fit, robust=TRUE)
    #print("Check whether we should include a factor in the model by checking for how many metabolites this factor is significant in the model: ")
    if (checkCovariate) {
      fit.con_co = limma::eBayes(fit)
      rlst_co = lapply(2:(ncol(design)-group.categories+1), function(i) {
          limma::topTable(fit.con_co, n=Inf, coef=colnames(design)[i])}) %>%
          setNames(colnames(design)[2:(ncol(design)-group.categories+1)])
      #n refers to the number of top-ranked genes to be returned.
    } else {rlst_co = NA}
    
    fit.con = limma::eBayes(fit, robust = TRUE, trend=T)
    rlst_interest = limma::topTable(fit.con, n=Inf, coef=levels(meta$Condition)[2]) %>%
        arrange(-t) 
    
    return(list(rslt_of_interest=rlst_interest, 
        rslt_co=rlst_co, data=dat, fitted.model=fit))
}


#' Run limma differential testing on CellFie task scores and save summaries
#'
#' Convenience wrapper that prepares CellFie task-score matrices at a
#' specified hierarchical depth, runs per-cancer-group limma testing
#' (via [limmaTest_CellFie]), and writes heatmaps and an Excel summary of
#' results to disk.
#'
#' @param depth character scalar. The hierarchical depth column name in
#'   `cfs` to run limma test on (for example "Depth3", "Depth2", or
#'   "Depth1").
#' @param cfs data.frame. CellFie output containing at minimum the columns
#'   `Sample`, `TaskScore`, and the Depth columns (e.g. `Depth1`, `Depth2`,
#'   `Depth3`) referenced by `depth`. If depth = Depth1/2, TaskScore need to be corresponding mean score values.
#' @param meta data.frame. Sample metadata with at least the columns
#'   `Sample`, `Cancer_type`, `Condition`, `Batch`, `Paired` (values
#'   'Paired'/'No'), and `Patient`. `Sample` must match column names in
#'   the reshaped `cfs` matrix.
#' @param OUTDIRV character scalar. Directory in which outputs (heatmaps and
#'   an Excel workbook summarizing limma results) will be written. The
#'   directory is created if it does not exist.
#' @param w numeric. Base width used when computing heatmap figure widths
#'   (default 9.5). The function adjusts width based on the number of
#'   samples.
#' @param height numeric. Figure height used for heatmap output (default
#'   23).
#' 
#' @importFrom dplyr %>%
#'
#' @return A list invisibly returned with elements:
#'   \item{cf5_Depth}{The reshaped data.frame used for analysis (rows = tasks at the chosen depth, columns = samples).}
#'   \item{rlst_list}{A named list of limma results returned by \code{limmaTest_CellFie} for each cancer type (only for types with sufficient data).}
#'
#' @details The function performs the following steps: (1) reshapes `cfs`
#' into a matrix of TaskScore by Sample for the requested `depth`; (2)
#' subsets metadata and decides whether to run a paired analysis based on
#' the `Paired` column and available pairs; (3) runs `limmaTest_CellFie`
#' per cancer type and produces per-cancer heatmaps (with and without
#' sample names); (4) writes a consolidated Excel workbook of the
#' top-table results for each cancer type.
#'
#' @examples
#' \dontrun{
#' # res <- runLimmaCellFie("Depth3", cfs, meta, OUTDIRV = "./results/limma", w = 9.5, height = 23)
#' # str(res)
#' }
#' @export
#' 
runLimmaCellFie = function(depth, cfs, meta, OUTDIRV, w=9.5, height=23) {
    dir.create(OUTDIRV, recursive = TRUE, showWarnings = FALSE)
    rownames(cfs) = NULL
    cfs = cfs %>% dplyr::mutate(TaskScore=log2(TaskScore+0.001))
    cf5_Depth = cfs[, c(depth, "Sample", "TaskScore")] %>% 
        tidyr::spread(key="Sample", value="TaskScore") %>% 
        as.data.frame() %>%
        tibble::column_to_rownames(var=depth) %>% 
        .[, intersect(colnames(.), meta$Sample)]

    meta$Condition = factor(meta$Condition, 
        levels=c("Normal", "Tumor"))

    Type = unique(meta$Cancer_type)
    nonPaired = c("")
    rlst_list_paired = lapply(setdiff(Type, nonPaired), function(ct) {
        print(ct)
        checkCovariate = FALSE
        meta_sub = meta[meta$Cancer_type==ct,]
        sub = cf5_Depth[, intersect(meta_sub$Sample, colnames(cf5_Depth))]
        meta_sub = filter(meta_sub, Sample %in% colnames(sub))
        if (length(unique(meta_sub$Batch))>1) {
          currentCovariate = "Batch"
        } else {currentCovariate = NULL}
        # We need to have at least two pairs in order to run paired analysis
        if ((meta_sub%>%filter(Paired=="Paired") %>% nrow())>2) {
          paired = TRUE
        } else {paired = FALSE}
        ## Remove features with 0 variance
        # invisible(capture.output({library(matrixStats)}))
        keep = matrixStats::rowSds(sub%>%as.matrix(), na.rm = T) > 0
        sub = sub[keep, ]
        if (length(unique(meta_sub$Condition))>1) {
            ## Paired analysis with only group in the model
            rlst = limmaTest_CellFie(dat=sub, paired=paired, 
                currentCovariate=currentCovariate, 
                checkCovariate=checkCovariate,
                meta=meta_sub, 
                paired_variable="Patient")
            sub = sub[rlst$rslt_of_interest %>%
                    tibble::rownames_to_column(var="Task")%>%
                    arrange(logFC) %>% pull(Task) %>% rev(), ] 
            sub = sub[, c(filter(meta_sub,Condition=="Normal")%>%pull(Sample),
                            filter(meta_sub,Condition=="Tumor")%>%pull(Sample))]
            colAnn = data.frame(Group = c(rep("Normal", nrow(filter(meta_sub,Condition=="Normal"))),
                            rep("Tumor", nrow(filter(meta_sub,Condition=="Tumor")))))
            rownames(colAnn) = colnames(sub)
            L = ncol(sub)
            pheatmap::pheatmap(sub,
                    show_colnames=TRUE,
                    color=viridisLite::magma(50),
                    breaks=seq(-2,2,len=50),
                    border_color="black",
                    scale="row",
                    annotation_col=colAnn,
                    main="Task Scores",
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    filename=paste0(OUTDIRV,"/TaskScore_heatmap_",ct,"_wSampleNames.pdf"), 
                    width=w+0.14*L, height=height)
            pheatmap::pheatmap(sub,
                    labels_col=rep("", ncol(sub)),
                    color=viridisLite::magma(50),
                    breaks=seq(-2,2,len=50),
                    border_color="black",
                    scale="row",
                    annotation_col=colAnn,
                    main="Task Scores",
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    filename=paste0(OUTDIRV,"/TaskScore_heatmap_",ct,".pdf"), 
                    width=w+0.14*L, height=height)
        } else {rlst = NA}
        rlst
    }) %>% setNames(setdiff(Type,nonPaired))
    rlst_list_paired = rlst_list_paired[!is.na(rlst_list_paired)]

    ### Save the limma results
    openxlsx::write.xlsx(lapply(names(rlst_list_paired), function(ct) {
            rlst_list_paired[[ct]]$rslt_of_interest %>% 
            tibble::rownames_to_column(var="Task")
        }) %>% setNames(names(rlst_list_paired)), 
        file=paste0(OUTDIRV, "/Limma_differentialAnalysis_result_", depth, ".xlsx"))

    return(list(cf5_Depth=cf5_Depth, rlst_list= rlst_list_paired))
}





#' Generate a task score heatmap across all samples and groups
#'
#' Create a combined heatmap of tasks that are frequently significant
#' (up- or down-regulated) across cancer types using per-group limma
#' results. The function extracts significant tasks from `rlst_list`,
#' orders them by frequency, and plots them across all samples present in
#' `cf5_Depth` with column annotations from `meta`.
#'
#' @param rlst_list named list. Output from `runLimmaCellFie` or
#'   `limmaTest_CellFie`; each element must contain a `rslt_of_interest`
#'   data.frame (topTable).
#' @param cf5_Depth matrix or data.frame. Task-by-sample matrix (rows =
#'   tasks at the chosen depth, columns = samples), typically the
#'   `cf5_Depth` returned by `runLimmaCellFie`.
#' @param OUTDIRV character. Directory where the heatmap files will be
#'   written. The directory is expected to exist or be creatable by the
#'   caller.
#' @param w numeric. Base width used to compute heatmap figure width.
#'   Defaults to 10.
#' @param h numeric. Height used for the heatmap output (defaults to
#'   23).
#' @param adj.P.Val.cutoff numeric. Adjusted p-value threshold used to
#'   select significant tasks (default 0.05).
#' @param meta data.frame. Sample metadata containing at least
#'   `Sample`, `Cancer_type`, and `Condition`, used for column
#'   annotations.
#' 
#' @importFrom dplyr %>%
#'
#' @return Invisibly returns the plotting object (if available) or NULL.
#' The main effect of the function is side-effect: PDF heatmap files are
#' written into `OUTDIRV`, and a consolidated list of significant tasks is
#' used for plotting.
#'
#' @examples
#' \dontrun{
#' makeTaskScoreHeatmap_CellFie(rlst_list, cf5_Depth, OUTDIRV = "./results", w = 10, h = 20, adj.P.Val.cutoff = 0.05, meta)
#' }
#' @export
#' 
makeTaskScoreHeatmap_CellFie = function(rlst_list, cf5_Depth, 
    OUTDIRV, w=10,h=23, adj.P.Val.cutoff = 0.05, meta) {
    sig.feature.list = c(lapply(names(rlst_list), function(ct) {
          rlst_list[[ct]]$rslt_of_interest %>% 
          tibble::rownames_to_column(var="Task") %>%
          filter(adj.P.Val<adj.P.Val.cutoff&logFC>0) %>% pull(Task)
      }) %>% setNames(paste0(names(rlst_list), "_up")),
      lapply(names(rlst_list), function(ct) {
          rlst_list[[ct]]$rslt_of_interest %>% 
          tibble::rownames_to_column(var="Task") %>%
          filter(adj.P.Val<adj.P.Val.cutoff&logFC<0) %>% pull(Task)
      }) %>% setNames(paste0(names(rlst_list), "_down")))

    ups = sig.feature.list[grep("_up", names(sig.feature.list),value=T)] %>% 
        unlist() %>% table() %>%
        .[order(.,decreasing = T)] %>% names()
    downs = sig.feature.list[grep("_down", names(sig.feature.list),value=T)] %>% 
        unlist() %>% table() %>%
        .[order(.,decreasing = T)] %>% names()

    Cancer_type2Sample = setNames(meta$Cancer_type, meta$Sample)
    Condition2Sample = setNames(meta$Condition, meta$Sample)
    colAnn = data.frame(
        Cancer = Cancer_type2Sample[colnames(cf5_Depth)],
        Group = Condition2Sample[colnames(cf5_Depth)] %>% 
            gsub("NAT", "Normal", .) %>% factor(., levels=c("Normal", "Tumor"))) 
    colAnn$Cancer = factor(colAnn$Cancer, levels = unique(colAnn$Cancer))
    rownames(colAnn) = colnames(cf5_Depth)

    L = ncol(cf5_Depth)
    cf5_Depth_2 = cf5_Depth[c(ups,downs)%>%unique(),]
    pheatmap::pheatmap(cf5_Depth_2,
        show_colnames=TRUE,
        color=viridisLite::magma(50),
        breaks=seq(-2,2,len=50),
        border_color="black",
        scale="row",
        annotation_col=colAnn,
        # annotation_colors=ann_colors,
        main="Task Scores",
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        filename=paste0(OUTDIRV,"/TaskScore_heatmap_all_wSampleNames.pdf"), 
        width=w+0.14*L, height=h)
}




#' Make limma dotplot to visualize the significant metabolic tasks/subsystems/systems across cancers
#'
#' Create a dot-plot summarizing differential behavior of tasks (or
#' subsystems/systems) across cancer types using per-group limma results.
#' Each point shows the effect size (t-statistic) and adjusted p-value for
#' a task in a given cancer; point shape encodes direction (Up/Down) and
#' color maps the test statistic.
#'
#' @param rlst_list named list. Output from `runLimmaCellFie` / `limmaTest_CellFie` containing `rslt_of_interest` data.frames per cancer.
#' @param adj.P.Val.cutoff numeric. Adjusted p-value cutoff for selecting significant tasks (e.g. 0.05).
#' @param w numeric or NULL. Figure width in inches; when NULL a sensible width is computed from the number/length of tasks.
#' @param h numeric or NULL. Figure height in inches; when NULL a sensible height is computed from the number of tasks.
#' @param OUTDIRV character. Directory where the dotplot PDF will be written.
#' @param Depth character. Which hierarchical column in `taskHierarchy` to use as the `Task` label (e.g. "Depth3").
#' @param lowfigure logical. If TRUE use compact theming suitable for small figures (default FALSE).
#' @param title character. Plot title (default: "Pan-cancer tasks").
#' @param suffix character. Filename suffix used when saving the plot (default: "Hu2025").
#' @param taskHierarchy data.frame. A mapping table that includes `Depth1`, `Depth2`, `Depth3`, and a `Task` column used to attach hierarchical labels.
#' @param nested logical. If TRUE facets are nested by `Depth1` (requires `ggh4x`).
#' @param filter logical. If TRUE filter tasks to the top frequent ones based on `probs_cutoff`.
#' @param probs_cutoff numeric. Quantile threshold used when `filter = TRUE` to select top tasks (default 0.35).
#' 
#' @importFrom dplyr %>%
#'
#' @return A ggplot object (invisibly returned) and a PDF file written to `OUTDIRV`.
#'
#' @details The function collates per-cancer `topTable` results into a long
#' table, annotates tasks with hierarchical labels from `taskHierarchy`,
#' and filters to tasks that are significant by `adj.P.Val.cutoff`. Points
#' are sized by adjusted p-value (smaller means more significant) and
#' colored by the test statistic `t`. When `filter = TRUE`, tasks are
#' reduced to those occurring above the `probs_cutoff` frequency quantile.
#'
#' @examples
#' \dontrun{
#' # p <- makeLimmaDotplot_CellFie(rlst_list, adj.P.Val.cutoff = 0.05, OUTDIRV = "./results", Depth = "Depth3", taskHierarchy = cfs)
#' # print(p)
#' }
#' @export
#' 
makeLimmaDotplot_CellFie = function(rlst_list,
    adj.P.Val.cutoff, w=NULL, h=NULL, OUTDIRV, 
    Depth="Depth3", lowfigure=F, title = "Pan-cancer tasks", 
    suffix="Hu2025", taskHierarchy=cfs, nested=T, filter=FALSE,
    probs_cutoff=0.35){
    sig.feature.list = c(lapply(names(rlst_list), function(ct) {
        rlst_list[[ct]]$rslt_of_interest %>% 
        tibble::rownames_to_column(var="Task") %>%
        filter(adj.P.Val<adj.P.Val.cutoff&logFC>0) %>% pull(Task)
      }) %>% setNames(paste0(names(rlst_list), "_up")),
      lapply(names(rlst_list), function(ct) {
        rlst_list[[ct]]$rslt_of_interest %>% 
        tibble::rownames_to_column(var="Task") %>%
        filter(adj.P.Val<adj.P.Val.cutoff&logFC<0) %>% pull(Task)
      }) %>% setNames(paste0(names(rlst_list), "_down")))

    ups = sig.feature.list[grep("_up", names(sig.feature.list),value=T)] %>% 
        unlist() %>% table() %>%
        .[order(.,decreasing = T)] %>% names()
    downs = sig.feature.list[grep("_down", names(sig.feature.list),value=T)] %>% 
        unlist() %>% table() %>%
        .[order(.,decreasing = T)] %>% names()

    taskHierarchy$Task = taskHierarchy[[Depth]]

    rlst_long = lapply(seq_along(rlst_list), function(i) {
        rlst = rlst_list[[i]]$rslt_of_interest
        rlst$Cancer = names(rlst_list)[i]
        rlst$Regulation = ifelse(rlst$logFC>0&rlst$adj.P.Val<adj.P.Val.cutoff, "Up", 
            ifelse (rlst$logFC<0&rlst$adj.P.Val<adj.P.Val.cutoff, "Down", "NonSig"))
        rlst %>% tibble::rownames_to_column(var="Task") 
    }) %>% Reduce(rbind, .) %>%
        dplyr::left_join(unique(taskHierarchy[, c("Depth1", "Depth2", "Depth3", "Task")])) %>%
            dplyr::mutate(logFC=round(logFC, 2)) %>%
            dplyr::mutate(Depth2 = stringr::str_sub(Depth2, 1, 3)) %>%
            dplyr::mutate(Depth1 = stringr::str_sub(Depth1, 1, 3)) %>% 
            filter(Task %in% c(ups, downs)) %>%
            filter(adj.P.Val < adj.P.Val.cutoff)

    rlst_long$Task = factor(rlst_long$Task,
            levels=rev(c(ups, downs) %>% unique()))
    rlst_long$t[rlst_long$t>5] = 5
    rlst_long$t[rlst_long$t< -5] = -5

    if (filter) {
        temp = rlst_long %>% dplyr::select(Task, Cancer, Regulation) %>% 
            unique() %>% filter(Regulation %in% c("Up", "Down")) %>%
            pull(Task) %>% table()
        tops = names(temp)[temp > quantile(temp, probs=probs_cutoff)]
        rlst_long = rlst_long %>% filter(Task %in% tops)
    }

    p = ggplot2::ggplot(data=rlst_long, 
            mapping=ggplot2::aes_string(x="Cancer", y="Task")) +
        ggplot2::geom_point(ggplot2::aes(size=adj.P.Val, color=t, shape=Regulation))
    if (nested) p = p+ggh4x::facet_nested(Depth1~., scales = "free_y", space = "free")
    p = p + ggplot2::scale_shape_manual(values = c("Up" = 18, "Down" = 20, "NonSig" = 17)) +
        ggplot2::ggtitle(title) +
        ggplot2::ylab("") + 
        ggplot2::scale_color_gradient2(low="skyblue", mid="snow2", high="purple", midpoint=0) +
        ggplot2::theme_bw(base_size = 15) + 
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 1000)) +
        ggplot2::scale_size(trans="reverse") +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            strip.text.y = ggplot2::element_text(size = 10, angle = 90, hjust=0.5),
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)) +
        ggplot2::labs(caption=paste0("adj.P.ValCutoff = ",adj.P.Val.cutoff,", pAdjustMethod = BH"))

    if (lowfigure) {
      p = p + ggplot2::theme(
              axis.title = ggplot2::element_blank(),
              legend.position = "bottom",
              legend.box = "vertical",         # Stack legends vertically
              legend.box.just = "left",        # Center the whole legend box
              legend.text = ggplot2::element_text(size = 7),
              legend.title = ggplot2::element_text(size = 11),
              legend.spacing.y = ggplot2::unit(0, "cm") # Reduce spacing between legend rows
            ) +
            ggplot2::guides(
              fill  = ggplot2::guide_legend(ncol=1, title.position = "top", byrow = F),
              shape = ggplot2::guide_legend(title.position = "left"),
              size  = ggplot2::guide_legend(title.position = "left"))
    }

    if (is.null(h)) {
        h = ((unique(rlst_long$Task) %>% length())/6) + 4.5
    }
    if (is.null(w)) {
        w = (nchar(unique(as.character(rlst_long$Task))) %>% max())/11 + 4.3
    }
    pdf(paste0(OUTDIRV,"/Dotplot_all_adj.P.ValCutoff", 
                    adj.P.Val.cutoff,"_",suffix,".pdf"), w=w, h=h)
        show(p)
    dev.off()

    return(p)
}


#' Create a volcano plot for limma differential results
#'
#' Produce a publication-ready volcano plot from a limma top-table-like
#' result. Points are colored by significance direction and optionally
#' labeled for the most extreme fold-changes. Point size encodes a
#' coarse depth annotation to visually separate features from different
#' hierarchical levels.
#'
#' @param resi data.frame-like. A limma topTable/result table containing at least the columns
#'   \code{Task}, \code{logFC}, \code{P.Value}, \code{adj.P.Val}, and \code{Depth}.
#' @param NAME character. A short name used for the plot title (and helpful when saving files).
#' @param adj.P.Val.cutoff numeric. Adjusted p-value threshold to define significance (default 0.05).
#' @param logFC.cutoff numeric. Log2 fold-change threshold used to highlight up- and down-regulated features (default 0).
#' @param lfcShrink logical. If \code{TRUE}, annotate the plot title to indicate log-fold-change shrinkage was applied. Default: \code{FALSE}.
#' 
#' @importFrom dplyr %>%
#'
#' @return A \code{ggplot} object (invisibly returned) showing log2 fold-change vs -log10(adjusted p-value).
#'
#' @details The function filters out rows with missing \code{adj.P.Val}, classifies features as "up", "down" or "nonSig" according to \code{adj.P.Val.cutoff} and \code{logFC.cutoff}, and scales point sizes by the \code{Depth} annotation (Depth3 = smallest point). The most extreme fold-changes (top 2% by absolute \code{logFC}) are labeled using \pkg{ggrepel}. Colors are set so up = red, down = blue, nonSig = black.
#'
#' @examples
#' \dontrun{
#' # p <- doVolcano_CellFie(resi, NAME = "MyContrast")
#' # print(p)
#' }
#' @export
#' 
doVolcano_CellFie = function(resi, NAME, 
    adj.P.Val.cutoff=0.05, logFC.cutoff=0, lfcShrink=F) {
    selectedCols = c("Task", "logFC", "P.Value", "adj.P.Val", "Depth")
    toExport = resi[, selectedCols]
    toExport = toExport[order(toExport$P.Value),]

    tmp = as.data.frame(toExport) %>%
        filter(!is.na(adj.P.Val)) %>%
        dplyr::mutate(sig = ifelse(adj.P.Val<adj.P.Val.cutoff&logFC>logFC.cutoff, "up", 
        ifelse(adj.P.Val<adj.P.Val.cutoff&(logFC< -logFC.cutoff), "down", "nonSig")))
    tmp$PointSize = ifelse(tmp$Depth == "Depth3", 1,
                        ifelse(tmp$Depth == "Depth2", 2, 3))
    tmpSignif = tmp %>%
        filter(tmp$logFC>quantile(abs(tmp$logFC),0.98)|(tmp$logFC< -quantile(abs(tmp$logFC),0.98))) %>%
        filter(sig!="nonSig")
    
    p = ggplot2::ggplot(data = tmp, ggplot2::aes(x=logFC, y=-log10(adj.P.Val), col=sig)) + 
        ggplot2::geom_point(alpha = 0.5, ggplot2::aes(size=PointSize)) + 
        ggplot2::xlab("Log2 Fold Change") +
        ggrepel::geom_text_repel(data = tmpSignif, 
            ggplot2::aes(x=logFC, y=-log10(adj.P.Val), label = Task), 
        color="cornsilk4") +
        ggplot2::geom_hline(yintercept = -log10(0.05), lty = 2) +
        ggplot2::geom_vline(xintercept = c(-0, 0), lty = 2) + 
        ggplot2::scale_colour_manual(values = c("down"="blue", "nonSig"="black", "up"="red")) + 
        ggplot2::theme(legend.position="none") +
        ggplot2::ggtitle(ifelse(lfcShrink, paste0(NAME, "_lfcShrink"), paste0(NAME, "")))
    return(p)
}

#' Generate a circular barplot to summarize the sample size per cancer.
#' @param meta meta data table, it must at least have columns, 'Condition' and 'Cancer_type'.
#' @param title character. Plot title (default: "Pan-cancer tasks").
#' @param OUT_DIR character. Directory where the dotplot PDF will be written.
#' @param h numeric or NULL. Figure height in inches
#' @param w numeric or NULL. Figure width in inches
#' 
#' @importFrom dplyr %>%
#' 
#' @export
#' 
SampleSizeCircularBarplot = function(meta, title="Pan-cancer tasks", OUT_DIR, h=4, w=4.5) {
    temp = meta %>% dplyr::count(Cancer_type, Condition) %>%
        tidyr::spread(Condition, n) %>%
        dplyr::mutate(Sum = ifelse(is.na(Normal), 0, Normal) + ifelse(is.na(Tumor), 0, Tumor)) %>% 
        dplyr::mutate(Cancer_type = paste0(Cancer_type, "\n(",Sum,")")) %>%
        dplyr::select(-Sum) %>%
        tidyr::gather(key="Condition", value="Count", -Cancer_type) 
    temp$Condition = factor(temp$Condition, levels=c("Tumor","Normal"))
    plot = temp %>%
        ggplot2::ggplot(., ggplot2::aes(x=Cancer_type, y=Count, fill=Condition)) +
        ggplot2::geom_bar(stat="identity", show.legend=FALSE) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                plot.margin = ggplot2::unit(c(0.5,0.5,0.5,0.5), "cm")) + 
        ggplot2::coord_polar() +
        ggplot2::ggtitle(title) 
    pdf(paste0(OUT_DIR, "SampleSizeCircularBarplot.pdf"), h=h, w=w)
        print(plot)
    dev.off()
}


#' Generate a boxplot to visualize the mean taskscores by cancer types
#'
#' Produce a boxplot showing mean CellFie task scores aggregated at the
#' task level (e.g. `Depth3`) across cancer types. The function annotates
#' cancer labels with the total sample count and saves the plot as a PDF.
#'
#' @param meta data.frame. Sample metadata containing at least the columns
#'   \code{Sample}, \code{Condition} and \code{Cancer_type}. Samples in
#'   \code{meta$Sample} must match column names in the reshaped CellFie
#'   matrix created from \code{cf}.
#' @param cf data.frame. CellFie output containing at minimum the columns
#'   \code{Sample}, \code{TaskScore}, and hierarchical depth columns
#'   (e.g. \code{Depth1}, \code{Depth2}, \code{Depth3}).
#' @param OUT_DIR character. Directory where output files (PDF) will be written.
#' @param w numeric. Figure width in inches (default 2.75).
#' @param h numeric. Figure height in inches (default 3.0).
#' @param title character. Plot title (default: "Pan-cancer tasks").
#' 
#' @importFrom dplyr %>%
#'
#' @return Invisibly returns the ggplot object after writing a PDF to
#'   \code{OUT_DIR}. The primary effect is the written PDF file
#'   \code{meanTaskScores_cancers.pdf}.
#'
#' @details The function left-joins \code{cf} with \code{meta}, computes
#' mean task scores per cancer type at \code{Depth3}, orders cancer types
#' by their mean activity, and renders a boxplot. Cancer labels include
#' the total sample count to make cohort sizes immediately visible.
#'
#' @examples
#' \dontrun{
#' ActivityBoxplot_byCancer(meta, cf, OUT_DIR = "./plots/", w = 3, h = 4, title = "My Study")
#' }
#' @export
#' 
ActivityBoxplot_byCancer = function(
    meta, cf, OUT_DIR, w=2.75, h=3.0, title="Pan-cancer tasks") {
    cf5 = cf %>% dplyr::left_join(meta) %>% 
        dplyr::mutate(Condition=gsub("Normal", "NAT", Condition)) %>%
        dplyr::mutate(Cancer_Condition=paste0(Cancer_type,"_",gsub("Normal","NAT",Condition)))
    
    B_Z_N0 = cf5 %>% group_by(Depth3, Cancer_type) %>%
        dplyr::summarise(meanTaskScore = mean(TaskScore, na.rm=T), .groups = "drop")  

    Cancer_type_orderZ = B_Z_N0 %>% unique() %>% group_by(Cancer_type) %>%
        dplyr::summarise(meanTaskScore = mean(meanTaskScore, na.rm=T), .groups = "drop") %>% 
        arrange(meanTaskScore) %>% pull(Cancer_type) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(Cancer_type=factor(Cancer_type, levels = Cancer_type_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=Cancer_type, y=meanTaskScore)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="Cancer_type", y="Mean Metab. Act.")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/meanTaskScores_cancers.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}


#' Generate a boxplot of mean task scores by system (Depth1)
#'
#' Summarize CellFie mean task scores at the system level (Depth1).
#' The function computes the mean score per system and task (Depth3),
#' orders systems by median activity and renders a boxplot. The plot is
#' saved as a PDF to \code{OUT_DIR} named \code{meanTaskScores_system.pdf}.
#'
#' @param cf data.frame. CellFie output containing at minimum the columns
#'   \code{Sample}, \code{TaskScore}, and hierarchical depth columns
#'   (e.g. \code{Depth1}, \code{Depth2}, \code{Depth3}).
#' @param OUT_DIR character. Directory where output files (PDF) will be written.
#' @param w numeric. Figure width in inches (default 2).
#' @param h numeric. Figure height in inches (default 3.5).
#' @param title character. Plot title (default: "Pan-cancer tasks").
#' 
#' @importFrom dplyr %>%
#'
#' @return Invisibly returns the ggplot object after writing a PDF to
#'   \code{OUT_DIR}. The primary effect is the written PDF file
#'   \code{meanTaskScores_system.pdf}.
#'
#' @examples
#' \dontrun{
#' ActivityBoxplot_bySystem(cf, OUT_DIR = "./plots/", w = 2, h = 3.5, title = "My Study")
#' }
#' @export
#' 
ActivityBoxplot_bySystem = function(
    cf, OUT_DIR, w=2, h=3.5, title="Pan-cancer tasks") {
    B_Z_N0 = cf %>% group_by(Depth1, Depth3) %>%
        dplyr::summarise(meanTaskScore = mean(TaskScore, na.rm=T), .groups = "drop")  

    System_orderZ = B_Z_N0 %>% unique() %>% group_by(Depth1) %>%
        dplyr::summarise(meanTaskScore = mean(meanTaskScore, na.rm=T), .groups = "drop") %>% 
        arrange(meanTaskScore) %>% pull(Depth1) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(System=factor(Depth1, levels = System_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=System, y=meanTaskScore)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="System", y="Mean Metab. Act.")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/meanTaskScores_system.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}


#' Barplot of the number of significantly altered tasks per cancer
#'
#' Summarize and visualize the number of significantly changed tasks
#' detected by limma per cancer cohort. For each element of
#' \code{rlst_list} (a named list produced by \code{runLimmaCellFie} or
#' \code{limmaTest_CellFie}), the function counts tasks with adjusted
#' p-value below \code{sig.cutoff} separately for up- and down-regulated
#' tasks (based on \code{logFC}). Results are presented as a vertical
#' bar chart where up-regulated counts are positive and down-regulated
#' counts are negative. The plot is written to a PDF file in
#' \code{OUT_DIR}.
#'
#' @param rlst_list named list. Each element should be a list or object
#'   containing a data.frame named \code{rslt_of_interest} with at least
#'   the columns \code{adj.P.Val} (adjusted p-value) and \code{logFC}.
#' @param sig.cutoff numeric. Adjusted p-value threshold used to define
#'   significance. Default is \code{0.05}.
#' @param w numeric. Figure width in inches (default 2.75).
#' @param h numeric. Figure height in inches (default 3.0).
#' @param OUT_DIR character. Directory where the PDF output
#'   (\code{SigDiffNr.pdf}) will be written.
#' 
#' @importFrom dplyr %>%
#'
#' @return Invisibly returns the ggplot object used to draw the figure.
#'   The primary side effect is the creation of the PDF file
#'   \code{SigDiffNr.pdf} inside \code{OUT_DIR}.
#'
#' @examples
#' \dontrun{
#' # rlst_list <- list(OV = list(rslt_of_interest = ov_df), BRCA = list(rslt_of_interest = brca_df))
#' SigNrBarplot(rlst_list, sig.cutoff = 0.05, w = 3, h = 4, OUT_DIR = "./plots/")
#' }
#' @export
#' 
SigNrBarplot = function(
    rlst_list, sig.cutoff=0.05, w=2.75, h=3.0, OUT_DIR){
    temp = c(sapply(rlst_list, function(cancer) {
        cancer$rslt_of_interest %>% 
                filter(adj.P.Val<sig.cutoff&logFC>0) %>% nrow()
        }) %>% setNames(paste0(names(rlst_list), "_up")),

       sapply(rlst_list, function(cancer) {
        cancer$rslt_of_interest %>% 
                filter(adj.P.Val<sig.cutoff&logFC<0) %>% nrow()
        }) %>% setNames(paste0(names(rlst_list), "_down")) 
    )
    B3_H = data.frame(Cancer_type = names(temp), 
            Regulation = names(temp),
            Count=as.numeric(temp)) %>%
        dplyr::mutate(Cancer_type=gsub("_.*","",Cancer_type)) %>%
        dplyr::mutate(Regulation=gsub(".*_","",Regulation)) %>%
        dplyr::mutate(Count = ifelse(Regulation=="up",Count,-Count)) %>% arrange(.,Count)
    B3_H$Cancer_type = factor(B3_H$Cancer_type, levels=unique(rev(B3_H$Cancer_type)))

    B3_H = ggplot2::ggplot(B3_H, ggplot2::aes(x=Cancer_type, y=Count, fill=Regulation)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values=c("up"="firebrick","down"="steelblue")) +
        ggplot2::geom_hline(yintercept=0, color="black") +
        ggplot2::labs(title = "Sig. diff. tasks",
            x="Cancer_type", y="# Significant Tasks")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))
    
    pdf(paste0(OUT_DIR, "/SigDiffNr.pdf"), w=w, h=h)
        print(B3_H)
    dev.off()
}


#' @title Visualize Differential Task Effects by System (Hierarchical CellFie Data)
#'
#' @description
#' Generates violin or boxplots showing the distribution of differential task effects
#' (t-statistics) aggregated by system level (Depth1) from hierarchical CellFie 
#' task annotations. Supports both single-level (cohort-based) and nested 
#' (study-cohort) data structures. Effects can optionally be weighted by task 
#' activity level and winsorized to handle outliers. Useful for comparing 
#' metabolic task dysregulation across different biological systems.
#'
#' @param rlst_list A named list of results from \code{limmaTest_CellFie()} or 
#'   \code{runLimmaCellFie()}. Supports two structures:
#'   \itemize{
#'     \item **Single-level**: Each element represents a cohort/cancer type and must 
#'       contain a \code{$rslt_of_interest} data frame with limma results
#'     \item **Hierarchical (depth=5)**: Nested structure with studies at top level,
#'       cancer types nested below, then \code{$rslt_of_interest} data frames
#'   }
#'   Each result data frame must contain at minimum: `t` (t-statistic), 
#'   `logFC`, `AveExpr`, and a significance column (default `P.Value` or `adj.P.Val`).
#'
#' @param cf A data.frame representing the CellFie task hierarchy, containing at 
#'   minimum the columns:
#'   \itemize{
#'     \item `Depth1`: System-level category (highest hierarchical level)
#'     \item `Depth3`: Individual task/pathway identifier (lowest level)
#'   }
#'   Maps task results to their parent system categories for aggregation.
#'   For hierarchical input data, can be a named list of hierarchy tables (one per study).
#'
#' @param OUT_DIR Character string specifying the output directory where the PDF
#'   will be saved. Must be a valid, writable directory path.
#'
#' @param w Numeric. Width of output PDF in inches. Default: 2.
#'
#' @param h Numeric. Height of output PDF in inches. Default: 3.5.
#'
#' @param title Character string for the plot title (typically dataset or analysis name).
#'   Default: "Pan-cancer tasks".
#'
#' @param sig.statistic Character string naming the significance column in result 
#'   data frames (e.g., "P.Value", "adj.P.Val"). Default: "P.Value".
#'
#' @param sig.cutoff Numeric. Significance threshold for filtering tasks.
#'   Only tasks with `sig.statistic < sig.cutoff` are included.
#'   Default: 0.05.
#'
#' @param weight.effect.by.gene Logical. When TRUE, weights the differential effect 
#'   (t-statistic) by the inverse rank of task activity level (AveExpr), emphasizing 
#'   effects in highly active tasks. Default: FALSE.
#'
#' @param winsorize.effect Logical. When TRUE, clips extreme weighted effect values 
#'   to the specified probability quantiles. Default: TRUE.
#'
#' @param prob Numeric. Probability level for symmetrical winsorization at both tails 
#'   (e.g., 0.05 clips at 5th and 95th percentiles). Ignored if 
#'   `winsorize.effect=FALSE`. Default: 0.05.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Detects data structure depth (single-level vs. hierarchical)
#'   \item Filters tasks for significance using the specified threshold
#'   \item Maps tasks to systems using the hierarchy table (`Depth3` to `Depth1`)
#'   \item Optionally weights effects by inverse activity rank
#'   \item Optionally winsorizes extreme weighted effects
#'   \item Orders systems by median effect across all tasks
#'   \item Creates visualization with:
#'     \itemize{
#'       \item X-axis: system categories (ordered by median effect)
#'       \item Y-axis: weighted differential effects
#'       \item Single-level data: boxplots
#'       \item Hierarchical data: violin plots with study-based color fill
#'     }
#' }
#' Violin plots effectively display the full distribution of task effects 
#' within each system, revealing both typical and extreme dysregulation patterns.
#'
#' @importFrom dplyr %>% left_join select group_by summarise arrange pull mutate
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot position_dodge 
#'   scale_fill_manual labs theme_minimal theme element_text
#'
#' @return Invisibly returns the ggplot2 object used to generate the visualization.
#'   Primary side effect is creation of "diffEffectBoxplot_system.pdf" in `OUT_DIR`.
#'
#' @examples
#' \dontrun{
#'   rlst_list <- list(
#'     OV = list(rslt_of_interest = ov_limma_results),
#'     BRCA = list(rslt_of_interest = brca_limma_results)
#'   )
#'   diffEffectBoxplot_bySystem(rlst_list, cf, OUT_DIR = "./plots/", 
#'                              w = 3, h = 4, title = "Metabolic Systems")
#' }
#' @export
#' 
diffEffectBoxplot_bySystem = function(
    rlst_list, cf, 
    OUT_DIR, w=2, h=3.5, title="Pan-cancer tasks", sig.statistic="P.Value", sig.cutoff=0.05,
    weight.effect.by.gene = FALSE, winsorize.effect = T, prob=0.05) {
    depth = function(x) {
        if (!is.list(x)) return(0)
        1 + max(sapply(x, depth))
        }
    if (depth(rlst_list)==5) {
        B_Z_N0 = lapply(names(rlst_list), function(study) {
            rlst = rlst_list[[study]]
            temp2 = lapply(names(rlst), function(c){
                temp = rlst[[c]]$rslt_of_interest
                # temp = temp[temp[[sig.statistic]] < sig.cutoff, ]
                temp$sig = ifelse(temp[[sig.statistic]]<sig.cutoff, TRUE, FALSE)
                temp$Cancer = c
                if (!"Depth3" %in% colnames(temp)) {
                    temp = as.data.frame(temp) %>% 
                        tibble::rownames_to_column(var="Depth3") 
                }
                temp$Source = study
                temp
            }) %>% Reduce(rbind, .) %>% as.data.frame()
            if (!"Depth1" %in% colnames(temp2)) {
                temp2 = temp2 %>% dplyr::left_join(cf[[study]] %>% 
                    dplyr::select(Depth1, Depth3) %>% unique())
            }
            temp2
            }) %>% Reduce(rbind, .) %>% as.data.frame()
    } else {
        B_Z_N0 = lapply(names(rlst_list), function(c){
            temp = rlst_list[[c]]$rslt_of_interest
            # temp = temp[temp[[sig.statistic]] < sig.cutoff, ]
            temp$sig = ifelse(temp[[sig.statistic]]<sig.cutoff, TRUE, FALSE)
            temp$Cancer = c
            as.data.frame(temp) %>% tibble::rownames_to_column(var="Depth3") 
        }) %>% Reduce(rbind, .) %>% as.data.frame() %>% 
        dplyr::left_join(cf %>% dplyr::select(Depth1, Depth3) %>% unique())
    }
    B_Z_N0$effect = B_Z_N0$t
    if (weight.effect.by.gene) {
        B_Z_N0$AveExpr = rev(rank(B_Z_N0$AveExpr, ties.method = "average")) #2^(B_Z_N0$AveExpr)
        # B_Z_N0$AveExpr_scaled = 
        #     (B_Z_N0$AveExpr-min(B_Z_N0$AveExpr))/(max(B_Z_N0$AveExpr)-min(B_Z_N0$AveExpr))
        B_Z_N0$effect_weighted = B_Z_N0$effect*B_Z_N0$AveExpr
    } else {
        B_Z_N0$effect_weighted = B_Z_N0$effect
    }

    if (winsorize.effect) {
        B_Z_N0$effect_weighted[B_Z_N0$effect_weighted>quantile(B_Z_N0$effect_weighted, probs=(1-prob))] = 
            quantile(B_Z_N0$effect_weighted, probs=(1-prob))
        B_Z_N0$effect_weighted[B_Z_N0$effect_weighted<quantile(B_Z_N0$effect_weighted, probs=prob)] = 
            quantile(B_Z_N0$effect_weighted, probs=prob) 
    }

    System_orderZ = B_Z_N0 %>% unique() %>% group_by(Depth1) %>%
        dplyr::summarise(meanEffect = median(effect_weighted, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(Depth1) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(System=factor(Depth1, levels = System_orderZ)) 
    if (depth(rlst_list)==5) {
        B_Z_N = B_Z_N %>% 
            ggplot2::ggplot(., ggplot2::aes(x=System, y=effect_weighted, fill=Source)) +
            ggplot2::geom_violin(position = position_dodge(width = 0.85),  # reduce spacing
                width = 2, trim=FALSE, alpha=0.75, col=NA)+
            # geom_jitter(size = 0.1, color = "grey", alpha = 0.5, 
            #             position = position_dodge(width = 0.85)) +
            # geom_point(data = subset(B_Z_N, sig == TRUE),
            #             aes(x = System, y = effect_weighted),
            #             color = "red", size = 0.1,
            #             position = position_dodge(width = 0.85)) +
            ggplot2::scale_fill_manual(values = 
                c("#E69F00", "#56B4E9","#999999","#F0E442","#009E73","#0072B2","#D55E00","#CC79A7"))
    } else {
        B_Z_N = B_Z_N %>%  
            ggplot2::ggplot(., ggplot2::aes(x=System, y=effect_weighted)) +
            ggplot2::geom_violin(trim=FALSE, fill="#E69F00", color=NA, alpha=0.6, width=1.1) #+
            # geom_jitter(size = 0.1, color = "grey", alpha = 0.5, 
            #             position = position_dodge(width = 0.85)) +
            # geom_point(data = subset(B_Z_N, sig == TRUE),
            #             aes(x = System, y = effect_weighted),
            #             color = "red", size = 0.1,
            #             position = position_dodge(width = 0.85))
    }

    B_Z_N = B_Z_N + ggplot2::labs(title = title,
            x="System", y="t")+
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
        ggplot2::theme_minimal() 
    if (depth(rlst_list)==5) {
        B_Z_N = B_Z_N + 
            ggplot2::theme(legend.position="right",
                axis.text.x = ggplot2::element_text(angle = 135, hjust=1, vjust=1))
    } else {
        B_Z_N = B_Z_N +
        ggplot2::theme(legend.position="right",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=1))
    }

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_system.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}


#' Boxplot of limma differential effect (logFC) by cancer type
#'
#' Create a boxplot summarizing the differential effect logFC (weighted by activity level) for each
#' cancer cohort. The function expects limma results in \code{rlst_list}
#' (output from \code{runLimmaCellFie} or \code{limmaTest_CellFie}) where
#' each element contains a data.frame \code{rslt_of_interest} with at
#' least the columns \code{t}, \code{adj.P.Val} and \code{logFC}. It
#' orders cancer types by median effect and writes a PDF plot to
#' \code{OUT_DIR}.
#'
#' @param rlst_list named list. Typical input is the output from
#'   \code{runLimmaCellFie} or \code{limmaTest_CellFie}; each element must
#'   contain a data.frame named \code{rslt_of_interest} with at least
#'   the columns \code{t}, \code{adj.P.Val}, and \code{logFC}.
#' @param OUT_DIR character. Directory where the PDF output file will be written.
#' @param w numeric. Figure width in inches (default 2.75).
#' @param h numeric. Figure height in inches (default 3).
#' @param title character. Plot title (default: "Pan-cancer tasks").
#' 
#' @importFrom dplyr %>%
#'
#' @return Invisibly returns the ggplot object used to draw the figure.
#'   The main side-effect is the written PDF file
#'   \code{diffEffectBoxplot_cancer.pdf} inside \code{OUT_DIR}.
#'
#' @examples
#' \dontrun{
#' diffEffectBoxplot_byCancer(rlst_list, OUT_DIR = "./plots/", w = 3, h = 4,
#'                           title = "Study Cancer Types")
#' }
#' @export
#' 
diffEffectBoxplot_byCancer = function(
    rlst_list, OUT_DIR, w=2.75, h=3.0, title="Pan-cancer tasks") {
    B_Z_N0 = lapply(names(rlst_list), function(cancer){
        temp = rlst_list[[cancer]]$rslt_of_interest
        temp$Cancer_type = cancer
        temp %>% tibble::rownames_to_column(var="Depth3") 
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>% unique()

    B_Z_N0$AveExpr = rev(rank(B_Z_N0$AveExpr, ties.method = "average")) #2^(B_Z_N0$AveExpr)
    # B_Z_N0$AveExpr_scaled = 
    #     (B_Z_N0$AveExpr-min(B_Z_N0$AveExpr))/(max(B_Z_N0$AveExpr)-min(B_Z_N0$AveExpr))
    B_Z_N0$logFC_weighted = B_Z_N0$logFC*B_Z_N0$AveExpr
    B_Z_N0$logFC_weighted[B_Z_N0$logFC_weighted>quantile(B_Z_N0$logFC_weighted, probs=0.95)] = 
        quantile(B_Z_N0$logFC_weighted, probs=0.95)
    B_Z_N0$logFC_weighted[B_Z_N0$logFC_weighted<quantile(B_Z_N0$logFC_weighted, probs=0.05)] = 
        quantile(B_Z_N0$logFC_weighted, probs=0.05) 

    Cancer_orderZ = B_Z_N0 %>% unique() %>% group_by(Cancer_type) %>%
        dplyr::summarise(meanEffect = median(logFC_weighted, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(Cancer_type) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(Cancer=factor(Cancer_type, levels = Cancer_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=Cancer, y=logFC_weighted)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="Cancer", y="Wei. Diff. effect")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_cancer.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}

# library
# library(tidyverse)
# library(viridis)

#' Create circular stacked barplot (radial barplot)
#'
#' Generates a circular (polar) stacked barplot grouped by `group`, where each
#' bar corresponds to an `individual` and segments are stacked by `observation`.
#' The plot is saved as a PDF file at `file.path(outdir, paste0(name, ".pdf"))`.
#'
#' @param data A data.frame containing at minimum the columns: `group`,
#'   `individual`, `observation`, and `value` (numeric). Additional columns are ignored.
#' @param outdir Character. Output directory where the PDF will be saved. The
#'   directory will be created if it does not exist.
#' @param name Character. Base name for the output PDF file (default: "SampleSizes").
#' @param h Numeric. Height (in inches) of the saved PDF (default: 10).
#' @param w Numeric. Width (in inches) of the saved PDF (default: 10).
#'
#' @return Invisibly returns the created `ggplot` object (invisible) after
#'   saving the PDF to disk.
#'
#' @details The function adds a small number of empty bars between groups for
#'   spacing, computes label angles to keep text readable, draws grid lines and
#'   group titles, and uses a viridis color scale for the stacked segments.
#'
#' @examples
#' # getCircularBarplot(df, outdir = "plots", name = "samples")
#'
#' @importFrom dplyr %>%
#' 
#' @export
#' 
getCircularBarplot = function(data, outdir, name="SampleSizes", h=10, w=10) {
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar = 2
    nObsType = nlevels(as.factor(data$observation))
    data$group = as.factor(data$group)
    to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
    colnames(to_add) = colnames(data)
    # to_add[(nrow(to_add) + 1):(nrow(to_add) + n), ] = NA
    to_add$group = rep(levels(data$group), each=empty_bar*nObsType )
    data = rbind(data, to_add)
    data = data %>% arrange(group, individual)
    # data$id = rep( seq(1, nrow(data)/nObsType) , each=nObsType)
    mapping = setNames(1:length(unique(paste0(data$individual, data$group))), 
        unique(paste0(data$individual, data$group)))
    data$id = mapping[paste0(data$individual, data$group)]

    # Get the name and the y position of each label
    data$value = as.numeric(data$value)
    label_data = data %>% group_by(id, individual) %>% summarize(tot=sum(value, na.rm=T))
    number_of_bar = nrow(label_data)
    angle = 90 - 360 * (label_data$id-0.5) /number_of_bar     
    # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust = ifelse( angle < -90, 1, 0)
    label_data$angle = ifelse(angle < -90, angle+180, angle)
    
    # prepare a data frame for base lines
    base_data = data %>% 
        group_by(group) %>% 
        summarize(start=min(id), end=max(id) - empty_bar) %>% 
        rowwise() %>% 
        dplyr::mutate(title=mean(c(start, end)))
    
    # prepare a data frame for grid (scales)
    grid_data = base_data
    grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 2
    # grid_data$start = grid_data$start - 1
    grid_data = grid_data[-1,]
    
    # Make the plot
    data$Condition = data$observation
    step = (floor(max(na.omit(label_data$tot))/200))*50
    p = ggplot2::ggplot(data) +      
        # Add the stacked bar
        ggplot2::geom_bar(ggplot2::aes(x=as.factor(id), y=value, fill=Condition), stat="identity", alpha=0.5) +
        viridis::scale_fill_viridis(discrete = TRUE, na.translate = FALSE) +
        # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
        ggplot2::geom_segment(data=grid_data, 
            aes(x = end, y = 0, xend = start, yend = 0), 
            colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
        ggplot2::geom_segment(data=grid_data, 
            aes(x = end, y = step, xend = start, yend = step), 
            colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
        ggplot2::geom_segment(data=grid_data, 
            aes(x = end, y = step*2, xend = start, yend = step*2), 
            colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
        ggplot2::geom_segment(data=grid_data, 
            aes(x = end, y = step*3, xend = start, yend = step*3), 
            colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
        ggplot2::geom_segment(data=grid_data, 
            aes(x = end, y = step*4, xend = start, yend = step*4), 
            colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
        # Add text showing the value of each 100/75/50/25 lines
        ggplot2::annotate("text", x = rep(max(data$id),5), 
            y = c(0, step, step*2, step*3, step*4), 
            label = as.character(c(0, step, step*2, step*3, step*4)),
            color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
        ggplot2::ylim(-150,max(label_data$tot, na.rm=T)) +
        # coord_cartesian(ylim = c(-150, max(label_data$tot, na.rm = TRUE))) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "right",
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(rep(0,4), "cm") 
        ) +
        ggplot2::coord_polar() +
        # Add labels on top of each bar
        ggplot2::geom_text(data=label_data,
            ggplot2::aes(x=id, y=tot, label=individual, hjust=hjust), 
            color="black", fontface="bold", alpha=0.6, size=5, 
            angle= label_data$angle, inherit.aes = FALSE ) +
        # Add base line information
        ggplot2::geom_segment(data=base_data, 
            ggplot2::aes(x = start, y = -5, xend = end, yend = -5), 
            colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
        ggplot2::geom_text(data=base_data, 
            ggplot2::aes(x = title, y = -18, label=group), 
            hjust=c(1,1,0,0), colour = "black", 
            alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

    pdf(paste0(outdir,"/", name ,".pdf"), h=h, w=w)
    print(p)
    dev.off()
}







#' @title Extract and Annotate CellFie Tasks by System
#'
#' @description
#' Converts gene entrez IDs to gene symbols using the mygene database query service,
#' then filters CellFie output tasks to extract only those belonging to specified
#' metabolic systems (e.g., Lipid Metabolism, Energy Metabolism). This is useful
#' for focusing downstream analyses on specific metabolic categories.
#'
#' @param taskInfo2 A data frame containing CellFie output, typically generated by
#'   \code{processCellFieOutput()}. Must contain at minimum:
#'   \itemize{
#'     \item `GeneAssociatedToEssentialRxnsTask`: Column with gene entrez IDs
#'       (numeric or character) to be converted to gene symbols
#'     \item `Depth1`: Column containing system/pathway classification
#'       (higher-level metabolic category)
#'   }
#'   Additional columns from CellFie output are preserved in the result.
#'
#' @param System Character vector specifying one or more system names to filter by.
#'   These should match values in the `Depth1` column (e.g., "Carbohydrate metabolism",
#'   "Lipid metabolism", "Energy metabolism").
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts unique entrez gene IDs from the input table
#'   \item Queries the mygene service to retrieve gene symbols for each ID
#'   \item Adds a new column `GeneAssociatedToEssentialRxnsTask_symbol` with
#'     the converted gene symbols mapped back to the original entries
#'   \item Filters rows to keep only those where `Depth1` matches one of the
#'     specified systems
#'   \item Returns the filtered and annotated table
#' }
#' This function requires internet connectivity for mygene queries.
#'
#' @importFrom mygene queryMany
#' @importFrom dplyr filter %>%
#'
#' @return A filtered data frame containing only tasks from the specified systems,
#'   with an added column `GeneAssociatedToEssentialRxnsTask_symbol` containing
#'   the gene symbols corresponding to the entrez IDs.
#'
#' @examples
#' \dontrun{
#'   # Extract lipid metabolism tasks
#'   lipid_tasks <- getSystem(taskInfo2, System = "Lipid metabolism")
#'   
#'   # Extract multiple systems
#'   energy_lipid <- getSystem(taskInfo2, System = c("Energy metabolism", "Lipid metabolism"))
#' }
#'
#' @export
#' 
getSystem = function(taskInfo2, System) {
    entrezID = unique(taskInfo2$GeneAssociatedToEssentialRxnsTask)
    res = mygene::queryMany(
        entrezID,
        scopes = "entrezgene",
        fields = "symbol",
        species = "human"
    )
    mapping = setNames(res$symbol, res$query)
    taskInfo2$GeneAssociatedToEssentialRxnsTask_symbol = 
        mapping[as.character(taskInfo2$GeneAssociatedToEssentialRxnsTask)]
    head(taskInfo2)
    system = taskInfo2 %>% filter(Depth1%in%System)
    return(system)
}


#' @title Create Heatmap of Differential Feature Expression Across Cohorts
#'
#' @description
#' Generates a comprehensive heatmap visualizing log2 fold-changes (logFC) of features
#' (genes, metabolites, tasks) across multiple cohorts or conditions. The heatmap includes
#' significance annotations (asterisks for FDR < 0.05) and task/pathway annotations in
#' the row margin. Features are sorted by total logFC across cohorts, and logFC values
#' are capped at specified thresholds for improved color gradient visualization.
#'
#' @param limma_rslt A named list of limma differential expression results, with one
#'   element per cohort or condition. Each element should be a data frame containing:
#'   \itemize{
#'     \item `Feature`: Feature identifiers (gene names, metabolite IDs, task names, etc.)
#'     \item `logFC`: log2 fold-change values
#'     \item `adj.P.Val`: Adjusted p-values (FDR) for significance testing
#'   }
#'   List names become cohort labels displayed on the heatmap x-axis.
#'
#' @param taskInfo A data frame containing feature annotation information, including:
#'   \itemize{
#'     \item `GeneAssociatedToEssentialRxnsTask_symbol`: Gene or task identifier
#'       (must match values in the feature row names)
#'     \item Column specified by `depth` parameter: Feature grouping/classification
#'       (e.g., pathway, system, or metabolic category)
#'   }
#'   Used to annotate rows in the heatmap with feature classification.
#'
#' @param outdir Character string specifying the output directory where the PDF
#'   heatmap will be saved. Default: "." (current directory).
#'
#' @param h Numeric or NULL. Height of output PDF in inches. When NULL, automatically
#'   computed based on number of features. Default: NULL.
#'
#' @param w Numeric or NULL. Width of output PDF in inches. When NULL, automatically
#'   computed based on number of cohorts and annotation categories. Default: NULL.
#'
#' @param max_logFC Numeric. Upper and lower bounds for logFC values displayed in
#'   the heatmap. LogFC values exceeding these bounds are clipped (winsorized) to
#'   improve color gradient visualization and focus on moderate effect sizes.
#'   Default: 3.
#'
#' @param depth Character string naming the column in `taskInfo` to use for row
#'   annotations (e.g., "Depth3" for CellFie task hierarchy). Default: "Depth3".
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Collects logFC values across all cohorts into a wide-format matrix
#'   \item Collects adjusted p-values in a parallel matrix for significance testing
#'   \item Sorts features by total logFC across cohorts (descending)
#'   \item Creates significance indicator matrix (asterisks for adj.P.Val < 0.05)
#'   \item Builds task/pathway annotation matrix for rows using specified `depth` column
#'   \item Clips logFC values to ±max_logFC range (winsorization)
#'   \item Generates heatmap with:
#'     \itemize{
#'       \item Color gradient: blue (negative logFC) → white (zero) → red (positive logFC)
#'       \item Asterisks overlaid for significant features
#'       \item Row annotations showing feature classification/pathway
#'       \item Features sorted by total effect magnitude
#'     }
#'   \item Automatically scales PDF dimensions based on data size
#' }
#' The heatmap effectively summarizes differential expression patterns across
#' multiple cohorts and highlights which features are consistently dysregulated.
#'
#' @importFrom dplyr select full_join filter mutate %>%
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @importFrom pheatmap pheatmap
#'
#' @return Invisibly returns the heatmap object. Primary side effect is creation
#'   of a PDF file in `outdir` containing the feature expression heatmap.
#'
#' @examples
#' \dontrun{
#'   # Create a heatmap of differential features across three cohorts
#'   visualizeDiffFeatures(
#'     limma_rslt = list(
#'       Cohort1 = cohort1_limma_results,
#'       Cohort2 = cohort2_limma_results,
#'       Cohort3 = cohort3_limma_results
#'     ),
#'     taskInfo = feature_annotations,
#'     outdir = "./plots/",
#'     max_logFC = 2.5,
#'     depth = "Depth3"
#'   )
#' }
#'
#' @export
#' 
visualizeDiffFeatures = function(
    limma_rslt=limma_rslt_list, 
    outdir=outDirZhou, 
    taskInfo=systemZhou,
    h=NULL, w=NULL,
    max_logFC=3,
    depth="Depth3") {

    logFCdat = lapply(names(limma_rslt), function(l) {
        lr = limma_rslt[[l]] %>% dplyr::select(Feature, logFC) %>%
            setNames(c("Feature", l))
    }) %>% Reduce(dplyr::full_join, .) %>% as.data.frame() %>%
        tibble::column_to_rownames(var="Feature")
    logFCdat = logFCdat[rev(order(rowSums(logFCdat))), ]

    padj_dat = lapply(names(limma_rslt), function(l) {
        lr = limma_rslt[[l]] %>% dplyr::select(Feature, adj.P.Val) %>%
            setNames(c("Feature", l))
    }) %>% Reduce(dplyr::full_join, .) %>% as.data.frame() %>%
        tibble::column_to_rownames(var="Feature")
    padj_dat = padj_dat[rownames(logFCdat), ]

    sig_mat = matrix("", nrow=nrow(logFCdat), ncol=ncol(logFCdat))
    sig_mat[padj_dat<0.05] = "*"

    taskInfo$Depth = taskInfo[[depth]]
    ann_row = taskInfo %>% dplyr::select(Depth, GeneAssociatedToEssentialRxnsTask_symbol) %>%
        na.omit() %>% unique() %>% 
        filter(GeneAssociatedToEssentialRxnsTask_symbol %in% as.character(rownames(logFCdat))) %>%
        dplyr::mutate(value = 1) %>%  # mark presence with 1
        tidyr::pivot_wider(names_from = Depth, 
            values_from = value,
            values_fill = list(value = 0)) %>%
        tibble::column_to_rownames(var="GeneAssociatedToEssentialRxnsTask_symbol")
    ann_row = ann_row[rownames(logFCdat), ]
    ann_colors = 
        replicate(ncol(ann_row), c("0" = "snow2", "1" = "#d95f02"), 
        simplify = FALSE) %>% setNames(colnames(ann_row))
    
    logFCdat[abs(logFCdat)>max_logFC] = max_logFC
    max_abs = max(abs(logFCdat), na.rm = TRUE)
    # Define breaks centered at 0
    breaks = seq(-max_abs, max_abs, length.out = 101)
    p1 = pheatmap::pheatmap(logFCdat,
            color = colorRampPalette(c("blue", "white", "red"))(100),
            breaks = breaks,
            # annotation_col = ann_row,   # top bar
            display_numbers = sig_mat,
            annotation_row = ann_row,
            annotation_colors = ann_colors,
            show_rownames=TRUE, 
            show_colnames=TRUE,
            fontsize = 10,
            border_color = "white",
            annotation_legend = FALSE,
            main=paste0("logFC (*: padj<0.05)"))
    if (is.null(w)) {
        w = 0.165*(ncol(logFCdat)+length(ann_colors))+2
    }
    if (is.null(h)) {
        h = 0.12*nrow(logFCdat)+0.05*max(nchar(colnames(ann_row))) +2.5
    }
    pdf(paste0(outdir, "/logFC_tumorVSnormal", depth, ".pdf"),
        h=h, w=w)
    print(p1)
    dev.off()
}


#' @title Perform System-Level Differential Analysis Across Multiple CellFie Cohorts
#'
#' @description
#' Comprehensive pipeline for extracting specific metabolic system genes/tasks from 
#' CellFie output and performing paired differential expression analysis across multiple 
#' cohorts (e.g., Hu2025, Zhou2020). The function filters genes associated with a 
#' specified metabolic system, conducts limma testing on paired samples (tumor vs. normal),
#' and returns consolidated differential analysis results with feature annotation. 
#' Supports flexible data source paths for reproducibility and modularity.
#'
#' @param System Character string specifying the metabolic system to analyze
#'   (e.g., "LIPIDS METABOLISM", "CARBOHYDRATES METABOLISM", "ENERGY METABOLISM"). 
#'   Must match a value in the `Depth1` column of CellFie output. Default: "LIPIDS METABOLISM".
#'
#' @param outDirHu Character string specifying the output directory for Hu2025
#'   cohort results. System tasks and differential analysis results will be saved here.
#'   Directory is created if it does not exist.
#'
#' @param outDirZhou Character string specifying the output directory for Zhou2020
#'   cohort results. System tasks and differential analysis results will be saved here.
#'   Directory is created if it does not exist.
#'
#' @param taskInfo2_path Character string specifying the path to CellFie task information
#'   RDS file for Hu2025 cohort (contains Depth1, Depth3 hierarchy and gene associations).
#'   Default: "../results/CellFie_Hu2025/CellFieOut/taskInfo.RDS".
#'
#' @param detailScoring_new_path2 Character string specifying the path to CellFie output
#'   CSV file for Zhou2020 cohort with task annotations and gene associations.
#'   Default: "../results/CellFie_Zhou2020/CellFieOut/all/detailScoring_new.csv".
#'
#' @param NdatHu_path Character string specifying the path to normalized expression matrix
#'   RDS file for Hu2025 cohort (samples as columns, genes as rows).
#'   Default: "../results/CellFie_Hu2025/Pdat2.RDS".
#'
#' @param metaHu_path Character string specifying the path to sample metadata RDS file
#'   for Hu2025 cohort (must contain Patient, Condition, Cancer_type columns).
#'   Default: "../results/CellFie_Hu2025/data/Hu2025/Meta_hu2025.RDS".
#'
#' @param NdatZhou_path Character string specifying the path to normalized expression matrix
#'   RDS file for Zhou2020 cohort.
#'   Default: "../results/CellFie_Zhou2020/Pdat2.RDS".
#'
#' @param metaZhou_path Character string specifying the path to sample metadata RDS file
#'   for Zhou2020 cohort.
#'   Default: "../results/CellFie_Zhou2020/data/Zhou2020/meta.RDS".
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item **System Extraction**: 
#'     \itemize{
#'       \item Loads CellFie task information for each cohort from specified paths
#'       \item Filters tasks to extract only those belonging to the specified system
#'       \item Removes duplicate entries and saves system task information as RDS and XLSX files
#'       \item Maps entrez gene IDs to gene symbols using mygene service (via getSystem)
#'     }
#'   \item **Data Preparation**:
#'     \itemize{
#'       \item Loads normalized expression data and sample metadata from specified paths
#'       \item Deduplicates genes (keeps first occurrence)
#'       \item Filters genes to include only those associated with system tasks
#'       \item Identifies paired samples (same patient with multiple conditions)
#'     }
#'   \item **Paired Differential Analysis**:
#'     \itemize{
#'       \item Runs paired limma t-tests for each cancer type within each cohort
#'       \item Only cancer types with >4 paired samples and 2+ conditions are analyzed
#'       \item Tests: limma with paired design (Patient as blocking factor)
#'       \item Maps entrez gene IDs to gene symbols for interpretability
#'     }
#'   \item **Result Consolidation**:
#'     \itemize{
#'       \item Saves system task information (RDS and XLSX) for each cohort
#'       \item Saves limma results (RDS and XLSX) for each cohort
#'       \item Exports metadata (RDS) for reproducibility
#'       \item Returns consolidated list with all results
#'     }
#' }
#' Flexible path parameters allow easy adaptation to different project structures and 
#' alternative CellFie outputs. Missing data is handled gracefully with NA results for 
#' cancer types not meeting analysis criteria.
#'
#' @importFrom dplyr filter mutate select %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom openxlsx write.xlsx
#'
#' @return A named list containing:
#'   \itemize{
#'     \item `limma_rslt_listHu`: Named list of limma results for each cancer type 
#'       in Hu2025 cohort (names = cancer types)
#'     \item `limma_rslt_listZhou`: Named list of limma results for each cancer type 
#'       in Zhou2020 cohort
#'     \item `systemHu`: Extracted system task information for Hu2025 (data frame with 
#'       gene associations and symbols)
#'     \item `systemZhou`: Extracted system task information for Zhou2020
#'   }
#'   Each element in the limma results lists contains a data frame with:
#'   \itemize{
#'     \item `Feature`: Gene symbol (human-readable identifier)
#'     \item `EntrezID`: Original entrez gene ID (row name from expression matrix)
#'     \item `logFC`: log2 fold-change (tumor vs. normal)
#'     \item `AveExpr`: Average expression level across samples
#'     \item `t`: t-statistic from paired limma test
#'     \item `P.Value`: Unadjusted p-value
#'     \item `adj.P.Val`: Adjusted p-value (FDR-corrected)
#'   }
#'
#' @section Output Files:
#' Creates output files in `outDirHu` and `outDirZhou`:
#'   \itemize{
#'     \item `systemTaskInfoHu.RDS/XLSX`: System task information for Hu2025
#'     \item `systemTaskInfoZhou.RDS/XLSX`: System task information for Zhou2020
#'     \item `limma_rslt_listHu.RDS/XLSX`: Differential analysis results for Hu2025
#'     \item `limma_rslt_listZhou.RDS/XLSX`: Differential analysis results for Zhou2020
#'     \item `metaHu.RDS`: Sample metadata for Hu2025 (copy from source)
#'     \item `metaZhou.RDS`: Sample metadata for Zhou2020 (copy from source)
#'   }
#'
#' @section Data Requirements:
#' Metadata tables must contain:
#'   \itemize{
#'     \item `sample.submitter_id` or `sample_id`: Sample/column identifier
#'     \item `Condition` or `sample_type`: Sample type (typically "Solid Tissue Normal" vs. "Primary Tumor")
#'     \item `Patient`: Patient identifier for pairing
#'     \item `Cancer_type`: Cancer type classification
#'   }
#' Expression matrices should have:
#'   \itemize{
#'     \item "genes" column with gene identifiers
#'     \item Sample identifiers matching metadata
#'   }
#'
#' @examples
#' \dontrun{
#'   # Analyze lipid metabolism system with default paths
#'   results <- systemDiffAnalysis(
#'     System = "LIPIDS METABOLISM",
#'     outDirHu = "./results/SystemAnalysis_Hu/",
#'     outDirZhou = "./results/SystemAnalysis_Zhou/"
#'   )
#'   
#'   # Analyze with custom data paths
#'   results <- systemDiffAnalysis(
#'     System = "ENERGY METABOLISM",
#'     outDirHu = "./analysis/energy_hu/",
#'     outDirZhou = "./analysis/energy_zhou/",
#'     taskInfo2_path = "/custom/path/taskInfo.RDS",
#'     NdatHu_path = "/custom/path/expression_data.RDS"
#'   )
#'   
#'   # Extract and use results
#'   hu_results <- results$limma_rslt_listHu
#'   zhou_results <- results$limma_rslt_listZhou
#'   
#'   # Visualize or perform downstream analysis
#'   visualizeDiffFeatures(limma_rslt = hu_results, taskInfo = results$systemHu)
#' }
#'
#' @export
#' 
systemDiffAnalysis = function(
    System="LIPIDS METABOLISM",
    outDirHu, outDirZhou,
    taskInfo2_path = "../results/CellFie_Hu2025/CellFieOut/taskInfo.RDS",
    detailScoring_new_path2 = "../results/CellFie_Zhou2020/CellFieOut/all/detailScoring_new.csv",
    NdatHu_path = "../results/CellFie_Hu2025/Pdat2.RDS",
    metaHu_path = "../results/CellFie_Hu2025/data/Hu2025/Meta_hu2025.RDS",
    NdatZhou_path = "../results/CellFie_Zhou2020/Pdat2.RDS",
    metaZhou_path = "../results/CellFie_Zhou2020/data/Zhou2020/meta.RDS"){
    print("Extract the specific METABOLISM genes ...")
    # Hu2025
    taskInfo2 = readRDS(taskInfo2_path)
    systemHu = getSystem(taskInfo2, System=System)
    systemHu = systemHu %>% unique()
    saveRDS(systemHu, paste0(outDirHu, "/systemTaskInfoHu.RDS"))
    openxlsx::write.xlsx(list(Hu2025=systemHu),
        paste0(outDirHu, "systemTaskInfoHu.xlsx"))
    rm(taskInfo2, OUTDIR)
    gc()

    # Zhou2020
    taskInfo2 = 
        read.csv(detailScoring_new_path2, sep="\t", header=TRUE) 
    head(taskInfo2)
    systemZhou = getSystem(taskInfo2, System=System)
    systemZhou = systemZhou %>% unique()
    saveRDS(systemZhou, paste0(outDirZhou, "/systemTaskInfoZhou.RDS"))
    openxlsx::write.xlsx(list(Zhou2020=systemZhou),
        paste0(outDirZhou, "systemTaskInfoZhou.xlsx"))
    rm(taskInfo2, OUTDIR)
    gc()

    print("Paired differential analysis ... ")
    NdatHu = readRDS(NdatHu_path)  
    systemHu = readRDS(paste0(outDirHu, "/systemTaskInfoHu.RDS")) 
    metaHu = readRDS(metaHu_path) 
    saveRDS(metaHu, paste0(outDirHu, "/metaHu.RDS"))

    NdatHu = NdatHu[!duplicated(NdatHu$genes),]
    rownames(NdatHu) = NULL
    NdatHu = NdatHu %>% tibble::column_to_rownames(var="genes")
    NdatHu[] = lapply(NdatHu, as.numeric)
    sharedGenes = 
        intersect(as.character(unique(systemHu$GeneAssociatedToEssentialRxnsTask)), rownames(NdatHu))
    patientPair = which(table(metaHu$Patient)==2) %>% names()
    limma_rslt_listHu = 
        lapply(unique(metaHu$Cancer_type), function(c) {
            print(c)
            meta_sub = metaHu %>% filter(Patient%in%patientPair) %>% filter(Cancer_type==c)
            if ((nrow(meta_sub)>4) & length(unique(meta_sub$Condition))==2) {
                limmaTest_CellFie(
                    dat=NdatHu[sharedGenes,], 
                    paired=TRUE, 
                    currentCovariate=NULL, 
                    checkCovariate=FALSE, 
                    meta=meta_sub, 
                    paired_variable="Patient")  
            } else {NA}
        }) %>%setNames(unique(metaHu$Cancer_type))
    limma_rslt_listHu = limma_rslt_listHu[!is.na(limma_rslt_listHu)]
    rm(NdatHu, metaHu, sharedGenes, OUTDIR, patientPair)
    gc()

    mappingHu = setNames(systemHu$GeneAssociatedToEssentialRxnsTask_symbol,
        systemHu$GeneAssociatedToEssentialRxnsTask)
    limma_rslt_listHu = lapply(limma_rslt_listHu, function(x) {
            x$rslt_of_interest %>% 
            tibble::rownames_to_column(var="EntrezID")%>%
            dplyr::mutate(Feature=as.character(mappingHu[EntrezID])) %>%
            dplyr::select(Feature, everything())})
    saveRDS(limma_rslt_listHu, paste0(outDirHu, "limma_rslt_listHu.RDS"))
    openxlsx::write.xlsx(
        limma_rslt_listHu, 
        file = paste0(outDirHu, "limma_reslut.xlsx"))

    NdatZhou = readRDS(NdatZhou_path)  
    systemZhou = readRDS(paste0(outDirZhou, "/systemTaskInfoZhou.RDS")) 
    metaZhou = readRDS(metaZhou_path)
    saveRDS(metaZhou, paste0(outDirZhou, "/metaZhou.RDS"))

    NdatZhou = NdatZhou[!duplicated(NdatZhou$genes),]
    rownames(NdatZhou) = NULL
    NdatZhou = NdatZhou %>% tibble::column_to_rownames(var="genes")
    NdatZhou[] = lapply(NdatZhou, as.numeric)
    sharedGenes = 
        intersect(as.character(unique(systemZhou$GeneAssociatedToEssentialRxnsTask)), rownames(NdatZhou))
    patientPair = which(table(metaZhou$Patient)==2) %>% names()
    limma_rslt_listZhou = 
        lapply(unique(metaZhou$Cancer_type), function(c) {
            print(c)
            meta_sub = metaZhou %>% filter(Patient%in%patientPair) %>% filter(Cancer_type==c)
            if ((nrow(meta_sub)>4) & length(unique(meta_sub$Condition))==2) {
                limmaTest_CellFie(
                    dat=NdatZhou[sharedGenes,], 
                    paired=TRUE, 
                    currentCovariate=NULL, 
                    checkCovariate=FALSE, 
                    meta=meta_sub, 
                    paired_variable="Patient")  
            } else {NA}
        }) %>%setNames(unique(metaZhou$Cancer_type))
    limma_rslt_listZhou = limma_rslt_listZhou[!is.na(limma_rslt_listZhou)]
    rm(NdatZhou, metaZhou, sharedGenes, OUTDIR, patientPair)
    gc()

    print("Save differential analysis results ...")
    mappingZhou = 
        setNames(systemZhou$GeneAssociatedToEssentialRxnsTask_symbol,
        systemZhou$GeneAssociatedToEssentialRxnsTask)
    limma_rslt_listZhou = lapply(limma_rslt_listZhou,function(x) {
            x$rslt_of_interest %>% 
            tibble::rownames_to_column(var="EntrezID")%>%
            dplyr::mutate(Feature=as.character(mappingZhou[EntrezID])) %>%
            dplyr::select(Feature, everything())})
    saveRDS(limma_rslt_listZhou, paste0(outDirZhou, "limma_rslt_listZhou.RDS"))
    openxlsx::write.xlsx(
        limma_rslt_listZhou, 
        file = paste0(outDirZhou, "limma_reslut.xlsx"))

    return(list(limma_rslt_listZhou = limma_rslt_listZhou,
            limma_rslt_listHu = limma_rslt_listHu,
            systemHu = systemHu,
            systemZhou = systemZhou))
}





#' @title Create Effect Size Boxplot for Metabolic Systems
#'
#' @description
#' Generates a combined violin and boxplot showing the distribution of differential
#' expression effect sizes (logFC, t-statistics, or other metrics) aggregated by
#' metabolic system across all features within a dataset. Systems are ordered by
#' median effect size, and extreme values are winsorized to improve visualization.
#' Useful for comparing metabolic dysregulation magnitude across different pathways.
#'
#' @param rslt_lists A named list of results from \code{systemDiffAnalysis()} or
#'   similar pipeline. Each element represents a metabolic system and contains nested
#'   lists with elements named "limma_rslt_listHu" or "limma_rslt_listZhou" 
#'   (depending on dataset). Each inner list contains differential expression results
#'   for one or more cancer types, with each element containing a data frame with
#'   columns `Feature` and the effect statistic (e.g., `logFC`, `t`).
#'
#' @param dataset Character string specifying which dataset to use: "Hu" (selects
#'   "limma_rslt_listHu") or "Zhou" (selects "limma_rslt_listZhou").
#'   Default: "Hu".
#'
#' @param effect.statistic Character string naming the column in differential
#'   expression results to use as the effect size metric (e.g., "logFC", "t",
#'   "AveExpr"). Default: "logFC".
#'
#' @param winsorize.quantil Numeric. Probability level for symmetrical winsorization
#'   of extreme effect values at both tails (e.g., 0.05 clips at 5th and 95th
#'   percentiles). Reduces extreme value influence on visualization. Default: 0.05.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Selects the appropriate dataset's limma results based on the `dataset` parameter
#'   \item Extracts effect sizes from each system and cancer type
#'   \item Aggregates effects across all features and cancer types within each system
#'   \item Orders systems by median effect size (ascending)
#'   \item Winsorizes extreme values to specified quantiles
#'   \item Creates a combined visualization with:
#'     \itemize{
#'       \item X-axis: metabolic systems (ordered by median effect)
#'       \item Y-axis: effect sizes (logFC or specified statistic)
#'       \item Violin plots (orange) showing full distribution
#'       \item Overlaid boxplots (grey) showing quartiles
#'       \item Red dashed line at zero for reference
#'     }
#'   \item Saves PDF output (hard-coded to `outDirSummary` global variable)
#' }
#' System names are simplified by replacing "METABOLISM" with "M." for compact
#' axis labels.
#'
#' @importFrom dplyr select mutate group_by summarize arrange pull %>% 
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs geom_hline 
#'   theme_minimal theme element_text
#'
#' @return Invisibly returns NULL. Primary side effect is creation of a PDF file
#'   in the global `outDirSummary` directory named:
#'   "effectSize_<effect.statistic>_system_<dataset>.pdf"
#'
#' @section Note:
#' This function depends on the global variable `outDirSummary` for output path.
#' Ensure this variable is defined in the calling environment.
#'
#' @examples
#' \dontrun{
#'   # Create effect size boxplot for Hu dataset
#'   outDirSummary <- "./results/Summary/"
#'   results <- systemDiffAnalysis(System = "LIPIDS METABOLISM", 
#'                                 outDirHu = "./hu/", outDirZhou = "./zhou/")
#'   
#'   makeEffectSizeBoxplot_system(results, dataset = "Hu", effect.statistic = "logFC")
#'   
#'   # Also works with t-statistics
#'   makeEffectSizeBoxplot_system(results, dataset = "Zhou", effect.statistic = "t")
#' }
#'
#' @export
#' 
makeEffectSizeBoxplot_system = function(
    rslt_lists, dataset="Hu", 
    effect.statistic="logFC",
    winsorize.quantil=0.05) {
    if (dataset=="Hu") {
        limma_rslt = "limma_rslt_listHu"
    } else if (dataset=="Zhou") {
        limma_rslt = "limma_rslt_listZhou"
    } else {stop("dataset must be 'Hu' or 'Zhou'")}

    dat4plot = lapply(names(rslt_lists), function(sname) {
        s = rslt_lists[[sname]]
        lapply(s[[limma_rslt]], function(c) {
            c$logFC = c[[effect.statistic]]
            c %>% dplyr::select(Feature, logFC)
        }) %>% Reduce(rbind, .) %>%
        dplyr::mutate(System = sname)
    }) %>% Reduce(rbind, .)

    dat4plot$System = gsub("METABOLISM", "M.", dat4plot$System)
    order = dat4plot %>% group_by(System) %>% 
        summarize(median=median(logFC)) %>% 
        arrange(median) %>% pull(System) %>% unique()
    dat4plot$System = factor(dat4plot$System, levels=order)
    quantileLower = quantile(dat4plot$logFC, probs=winsorize.quantil)
    quantileUpper = quantile(dat4plot$logFC, probs=1-winsorize.quantil)
    dat4plot$logFC[dat4plot$logFC < quantileLower] = quantileLower
    dat4plot$logFC[dat4plot$logFC > quantileUpper] = quantileUpper

    pdf(paste0(outDirSummary, "/effectSize_",effect.statistic,"_system_", dataset, ".pdf"), h=4, w= 3)
    print(dat4plot %>% ggplot2::ggplot(., ggplot2::aes(x=System, y=logFC)) +
        ggplot2::geom_violin(trim=FALSE, fill="#E69F00", color="#E69F00")+
        ggplot2::geom_boxplot(width=0.4, fill="grey", alpha = 0.4, color="black", linewidth = 0.1, outlier.size = 0.5) + 
        ggplot2::labs(title = "",
            x="System", y=effect.statistic)+
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()
}





#' @title Create Effect Size Boxplot for Cancer Types
#'
#' @description
#' Generates a boxplot showing the distribution of differential expression effect sizes 
#' (logFC, t-statistics, or other metrics) across different cancer types. Effects are 
#' aggregated across all metabolic systems and features within each cancer type. Cancer 
#' types are ordered by median effect size, and extreme values are winsorized to improve 
#' visualization. Useful for comparing overall metabolic dysregulation magnitude between 
#' different cancer types.
#'
#' @param rslt_lists A named list of results from \code{systemDiffAnalysis()} or similar 
#'   pipeline. Each element represents a metabolic system and contains nested lists with 
#'   elements named "limma_rslt_listHu" or "limma_rslt_listZhou" (depending on dataset). 
#'   Each inner list contains differential expression results for multiple cancer types.
#'   The function extracts cancer type names from the first system's results.
#'
#' @param dataset Character string specifying which dataset to use: "Hu" (selects 
#'   "limma_rslt_listHu") or "Zhou" (selects "limma_rslt_listZhou"). 
#'   Default: "Hu".
#'
#' @param effect.statistic Character string naming the column in differential 
#'   expression results to use as the effect size metric (e.g., "logFC", "t", 
#'   "AveExpr"). Default: "logFC".
#'
#' @param winsorize.quantil Numeric. Probability level for symmetrical winsorization 
#'   of extreme effect values at both tails (e.g., 0.05 clips at 5th and 95th 
#'   percentiles). Reduces extreme value influence on visualization. Default: 0.05.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Selects the appropriate dataset's limma results based on the `dataset` parameter
#'   \item Extracts cancer type names from the first system's results
#'   \item Aggregates effect sizes across all metabolic systems and features for each cancer type
#'   \item Orders cancer types by median effect size (ascending)
#'   \item Winsorizes extreme values to specified quantiles
#'   \item Dynamically scales PDF width based on number of cancer types
#'   \item Creates boxplot with:
#'     \itemize{
#'       \item X-axis: cancer types (ordered by median effect)
#'       \item Y-axis: effect sizes (logFC or specified statistic)
#'       \item Orange boxplots showing quartile distributions
#'       \item Red dashed line at zero for reference
#'     }
#'   \item Saves PDF output (hard-coded to `outDirSummary` global variable)
#' }
#' Unlike \code{makeEffectSizeBoxplot_system()}, this function aggregates across
#' all systems to show cancer-specific metabolic dysregulation patterns.
#'
#' @importFrom dplyr select mutate group_by summarize arrange pull %>%
#' @importFrom ggplot2 ggplot aes geom_boxplot labs geom_hline theme_minimal theme element_text
#'
#' @return Invisibly returns NULL. Primary side effect is creation of a PDF file 
#'   in the global `outDirSummary` directory named: 
#'   "effectSize_<effect.statistic>_Cancer_<dataset>.pdf"
#'
#' @section Note:
#' This function depends on the global variable `outDirSummary` for output path. 
#' Ensure this variable is defined in the calling environment. PDF width is 
#' automatically scaled as 0.15 inches per cancer type plus 1 inch base width.
#'
#' @examples
#' \dontrun{
#'   # Create effect size boxplot for cancer types in Hu dataset
#'   outDirSummary <- "./results/Summary/"
#'   results <- systemDiffAnalysis(System = "LIPIDS METABOLISM", 
#'                                 outDirHu = "./hu/", outDirZhou = "./zhou/")
#'   
#'   makeEffectSizeBoxplot_cancer(results, dataset = "Hu", effect.statistic = "logFC")
#'   
#'   # Also works with t-statistics
#'   makeEffectSizeBoxplot_cancer(results, dataset = "Zhou", effect.statistic = "t")
#' }
#'
#' @export
#' 
makeEffectSizeBoxplot_cancer = function(
    rslt_lists, dataset="Hu", 
    effect.statistic="logFC",
    winsorize.quantil=0.05) {
    if (dataset=="Hu") {
        limma_rslt = "limma_rslt_listHu"
    } else if (dataset=="Zhou") {
        limma_rslt = "limma_rslt_listZhou"
    } else {stop("dataset must be 'Hu' or 'Zhou'")}

    cancers = names(rslt_lists[[1]][[limma_rslt]])
    dat4plot = lapply(cancers, function(c) {
        lapply(names(rslt_lists), function(sname) {
            cdat = rslt_lists[[sname]][[limma_rslt]][[c]]
            cdat$logFC = cdat[[effect.statistic]]
            cdat %>% dplyr::select(Feature, logFC)
        }) %>% Reduce(rbind, .) %>%
        dplyr::mutate(Cancer = c)
    }) %>% Reduce(rbind, .)

    order = dat4plot %>% group_by(Cancer) %>% 
        summarize(median=median(logFC)) %>% 
        arrange(median) %>% pull(Cancer) %>% unique()
    dat4plot$Cancer = factor(dat4plot$Cancer, levels=order)
    quantileLower = quantile(dat4plot$logFC, probs=winsorize.quantil)
    quantileUpper = quantile(dat4plot$logFC, probs=1-winsorize.quantil)
    dat4plot$logFC[dat4plot$logFC < quantileLower] = quantileLower
    dat4plot$logFC[dat4plot$logFC > quantileUpper] = quantileUpper

    w = length(order)*0.15 + 1
    pdf(paste0(outDirSummary, "/effectSize_",effect.statistic,"_Cancer_", dataset, ".pdf"), h=4, w= w)
    print(dat4plot %>% ggplot2::ggplot(., ggplot2::aes(x=Cancer, y=logFC)) +
        # ggplot2::geom_violin(trim=FALSE, fill="#E69F00", color="white")+
        ggplot2::geom_boxplot(width=0.75, fill="#E69F00", outlier.size = 0.5, line.width=0.3) + 
        ggplot2::labs(title = "",
            x="Cancer", y=effect.statistic)+
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()
}




#' @title Create Expression Level Boxplot for Metabolic Systems
#'
#' @description
#' Generates a combined violin and boxplot showing the distribution of average feature 
#' expression levels (genes, metabolites, or tasks) aggregated by metabolic system. 
#' Features are extracted from differential expression results and their mean expression 
#' values are computed from normalized expression data. Systems are ordered by median 
#' expression level, and extreme values are winsorized. Useful for assessing baseline 
#' abundance/activity of different metabolic systems in the dataset.
#'
#' @param NormExpr A data frame or matrix of normalized expression values with gene/feature 
#'   identifiers in a column named "genes" (which will be converted to row names) and 
#'   samples as columns. Typically output from normalization pipelines. Can include 
#'   duplicate gene entries, which are automatically deduplicated (first occurrence kept).
#'
#' @param rslt_lists A named list of results from \code{systemDiffAnalysis()} or similar 
#'   pipeline. Each element represents a metabolic system and contains nested lists with 
#'   elements named "limma_rslt_listHu" or "limma_rslt_listZhou" (depending on dataset). 
#'   Each inner element should contain a data frame with columns `Feature`, `EntrezID`, 
#'   and other limma outputs.
#'
#' @param dataset Character string specifying which dataset to use: "Hu" (selects 
#'   "limma_rslt_listHu") or "Zhou" (selects "limma_rslt_listZhou"). 
#'   Default: "Hu".
#'
#' @param winsorize.quantil Numeric. Probability level for symmetrical winsorization 
#'   of extreme mean expression values at both tails (e.g., 0.05 clips at 5th and 95th 
#'   percentiles). Reduces extreme value influence on visualization. Default: 0.05.
#'
#' @param outDirSummary Character string specifying the output directory where the PDF 
#'   will be saved. Must be a valid, writable directory path.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Deduplicates genes in the expression matrix (keeps first occurrence)
#'   \item Converts gene column to row names
#'   \item Selects the appropriate dataset's differential expression results
#'   \item For each system and cancer type combination:
#'     \itemize{
#'       \item Extracts feature IDs from limma results
#'       \item Subsets normalized expression data to these features
#'       \item Computes mean expression across all samples for each feature
#'       \item Maps entrez IDs to feature symbols
#'     }
#'   \item Aggregates mean expression values across all features and cancer types per system
#'   \item Orders systems by median expression level (ascending)
#'   \item Winsorizes extreme values to specified quantiles
#'   \item Creates combined visualization with:
#'     \itemize{
#'       \item X-axis: metabolic systems (ordered by median expression)
#'       \item Y-axis: mean feature expression levels
#'       \item Orange violin plots showing full distribution
#'       \item Overlaid boxplots for quartile visualization
#'     }
#'   \item Saves PDF output to `outDirSummary`
#' }
#' System names are simplified by replacing "METABOLISM" with "M." for compact axis labels.
#'
#' @importFrom dplyr select mutate group_by summarize arrange pull %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_minimal theme element_text
#'
#' @return Invisibly returns NULL. Primary side effect is creation of a PDF file 
#'   in `outDirSummary` named "FeatureAbundance_system_<dataset>.pdf"
#'
#' @examples
#' \dontrun{
#'   # Create expression level boxplot for systems in Hu dataset
#'   outDirSummary <- "./results/Summary/"
#'   results <- systemDiffAnalysis(System = "LIPIDS METABOLISM", 
#'                                 outDirHu = "./hu/", outDirZhou = "./zhou/")
#'   
#'   # Assume NdatHu is normalized expression data
#'   makeExprBoxplot_system(NormExpr = NdatHu, rslt_lists = results, 
#'                          dataset = "Hu", outDirSummary = outDirSummary)
#'   
#'   makeExprBoxplot_system(NormExpr = NdatZhou, rslt_lists = results, 
#'                          dataset = "Zhou", outDirSummary = outDirSummary)
#' }
#'
#' @export
#' 
makeExprBoxplot_system = function(
    NormExpr=NdatHu, 
    rslt_lists,
    dataset="Hu", 
    winsorize.quantil=0.05,
    outDirSummary) {

    if (dataset=="Hu") {
        limma_rslt = "limma_rslt_listHu"
    } else if (dataset=="Zhou") {
        limma_rslt = "limma_rslt_listZhou"
    } else {stop("dataset must be 'Hu' or 'Zhou'")}
    NormExpr = NormExpr[!duplicated(NormExpr$genes), ] 
    rownames(NormExpr) = NULL
    NormExpr = NormExpr %>% 
        tibble::column_to_rownames(var="genes")
    dat4plot = lapply(names(rslt_lists), function(sname) {
        cancers = names(rslt_lists[[sname]][[limma_rslt]])
        lapply(cancers, function(c) {
            features = rslt_lists[[sname]][[limma_rslt]][[c]][["EntrezID"]]
            dat = NormExpr[features, ]
            rownames(dat) = NULL
            rownames(dat) = rslt_lists[[sname]][[limma_rslt]][[c]][["Feature"]]
            data.frame(Feature=rownames(dat), meanExpr = rowMeans(dat), System=sname)
        }) %>% Reduce(rbind, .)
    }) %>% Reduce(rbind, .)
    dat4plot$System = gsub("METABOLISM", "M.", dat4plot$System)
    order = dat4plot %>% group_by(System) %>% 
        summarize(median=median(meanExpr)) %>% 
        arrange(median) %>% pull(System) %>% unique()
    dat4plot$System = factor(dat4plot$System, levels=order)
    quantileLower = quantile(dat4plot$meanExpr, probs=winsorize.quantil)
    quantileUpper = quantile(dat4plot$meanExpr, probs=1-winsorize.quantil)
    dat4plot$meanExpr[dat4plot$meanExpr < quantileLower] = quantileLower
    dat4plot$meanExpr[dat4plot$meanExpr > quantileUpper] = quantileUpper

    pdf(paste0(outDirSummary, "/FeatureAbundance_system_", dataset, ".pdf"), h=4, w= 2)
    print(dat4plot %>% ggplot2::ggplot(., ggplot2::aes(x=System, y=meanExpr)) +
        ggplot2::geom_violin(trim=FALSE, fill="#E69F00", color="white", width=1)+
        ggplot2::geom_boxplot(width=0.1, line.width=0.1) + 
        ggplot2::labs(title = "",
            x="System", y="meanExpr")+
        # ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0)))
    dev.off()
}

 
#' @title Build and Visualize Co-expression Network for Metabolic System
#'
#' @description
#' Constructs a co-expression network for genes/proteins associated with a specific 
#' metabolic system or pathway. The function filters high-correlation gene pairs 
#' (Spearman correlation > 0.6), identifies hub genes based on network connectivity 
#' (hubness), and generates a publication-ready network visualization with node sizes 
#' scaled by hubness and hub nodes highlighted distinctly.
#'
#' @param dat A data frame or matrix of expression data with genes/proteins as rows 
#' and samples as columns. Row names should contain gene identifiers.
#' Default: \code{NdatHu}
#'
#' @param taskInfo A data frame containing task/metabolic system annotations. 
#' Default: \code{taskInfo2}
#'
#' @param System A character string specifying the metabolic system name 
#' (e.g., "Glycolysis", "TCA Cycle"). Required.
#'
#' @param outdir A character string specifying the output directory path where 
#' the network visualization PDF will be saved. The function will create 
#' a \code{Network/} subdirectory if it doesn't exist.
#'
#' @param suffix A character string to append to the output file name for 
#' distinguishing multiple runs. Default: \code{""}
#'
#' @param framework A character string specifying the analysis framework. 
#' Options: \code{"CellFie"} (default) or KEGG-based analysis.
#' Default: \code{"CellFie"}
#'
#' @param keggTable An optional data frame containing KEGG pathway mappings 
#' with columns \code{class} (pathway/system name) and \code{gene_symbol}. 
#' If \code{NULL} and framework is not "CellFie", the function retrieves 
#' KEGG data automatically. Default: \code{NULL}
#'
#' @param size.index A numeric multiplier controlling overall network plot dimensions 
#' (height and width). Useful for large networks. Default: \code{1}
#'
#' @param vertex.size.index A numeric multiplier for scaling vertex (node) sizes 
#' based on hubness/connectivity. Larger values produce more size variation. 
#' Default: \code{10}
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Retrieves genes/proteins associated with the specified metabolic system
#'   \item Calculates Spearman correlation matrix across samples
#'   \item Identifies hub genes (top 3% or minimum 3 genes) with highest connectivity
#'   \item Filters edges using correlation threshold (default: |r| > 0.6)
#'   \item Builds igraph network object with edge weights and node hubness scores
#'   \item Visualizes using Fruchterman-Reingold layout
#'   \item Node colors: red for hubs, orange for non-hubs
#'   \item Node sizes scaled by within-module connectivity (hubness)
#'   \item Edge widths and colors reflect correlation strength and sign
#' }
#'
#' @return
#' A numeric vector named \code{kWithin} representing the hubness 
#' (within-module degree Z-score approximation) for each gene in the network. 
#' Higher values indicate genes with more connections to other genes in the module.
#'
#' @examples
#' \dontrun{
#'   # Visualize co-expression network for glycolysis pathway
#'   hubness <- getSystemNetwork(
#'     dat = expression_data,
#'     taskInfo = pathway_annotations,
#'     System = "Glycolysis",
#'     outdir = "./results/",
#'     framework = "CellFie",
#'     size.index = 1.5
#'   )
#'   
#'   # View top hub genes
#'   head(sort(hubness, decreasing = TRUE), n = 10)
#' }
#'
#' @import igraph
#' @import dplyr
#' @import scales
#' @importFrom tibble column_to_rownames
#'
#' @seealso
#' \code{\link{getSystem}} for retrieving system-specific genes,
#' \code{\link{getdb_metabolism}} for metabolic pathway databases
#'
#' @author Yuan Li
#'
#' @export
#'
getSystemNetwork = function(
    dat=NdatHu, taskInfo=taskInfo2, System, outdir, suffix="",
    framework="CellFie", keggTable=NULL, size.index=1,
    vertex.size.index=10, meta=NULL, 
    r2_threshold=0.8) {
    
    estimate_scale_free_fit <- function(cor_mat, powers = 1:20, n_bins = 10) {
        cor_mat[is.na(cor_mat)] <- 0
        results <- data.frame(Power = powers, R2 = NA, Slope = NA)
        
        for (i in seq_along(powers)) {
            beta <- powers[i]
            
            # 1. Calculate Weighted Adjacency Matrix
            adj_mat <- abs(cor_mat)^beta
            
            # 2. Calculate connectivity (kWithin)
            kWithin <- rowSums(adj_mat) - 1
            
            # Avoid zero connectivity errors by adding an infinitesimal value
            kWithin[kWithin == 0] <- 1e-5 
            
            # 3. Discretize kWithin into logarithmic frequency bins
            # This mirrors the internal math of the standard WGCNA package
            binned_k <- cut(log10(kWithin), breaks = n_bins)
            bin_counts <- table(binned_k)
            
            # Calculate P(k) - the probability distribution of connectivity
            p_k <- as.numeric(bin_counts) / sum(bin_counts)
            
            # Calculate the midpoints of the connectivity bins
            bin_breaks <- seq(min(log10(kWithin)), max(log10(kWithin)), length.out = n_bins + 1)
            bin_midpoints <- 10^(bin_breaks[-1] - diff(bin_breaks)/2)
            
            # 4. Filter out empty bins to prevent log(0) errors
            valid <- p_k > 0
            if (sum(valid) > 2) {
            log_midpoints <- log10(bin_midpoints[valid])
            log_p_k <- log10(p_k[valid])
            
            # 5. Fit the linear model: log(P(k)) ~ log(k)
            fit <- lm(log_p_k ~ log_midpoints)
            
            results$R2[i] <- summary(fit)$r.squared
            results$Slope[i] <- coef(fit)[2]
            }
        }
        return(results)
    }

    if (framework == "CellFie") {
        systemHu = getSystem(taskInfo, System=System)
        sub = dat[rownames(dat) %in% systemHu$GeneAssociatedToEssentialRxnsTask,]
        mapping = setNames(
            systemHu$GeneAssociatedToEssentialRxnsTask_symbol, 
            as.character(systemHu$GeneAssociatedToEssentialRxnsTask))
        sub$symbol = as.character(mapping[rownames(sub)])
        sub = sub[!duplicated(sub$symbol), ]
        rownames(sub) = NULL
        sub = sub %>% tibble::column_to_rownames(var="symbol")
    } else if (!is.null(keggTable)){
            keggTable_sub = filter(keggTable, class==System)
            sub = dat[dat$symbol %in% keggTable_sub$gene_symbol,]
            sub = sub[!duplicated(sub$symbol), ]
            rownames(sub) = NULL
            sub = sub %>% tibble::column_to_rownames(var="symbol")
    } else {
            kegg_metab_dbs = 
                getdb_metabolism(databaseDIR=paste0(outdir, "/KEGG_metabolism/"))
            keggTable_sub = filter(kegg_metab_dbs$kegg_metab_db, class==System)
            sub = dat[dat$symbol %in% keggTable_sub$gene_symbol,]
            sub = sub[!duplicated(sub$symbol), ]
            rownames(sub) = NULL
            sub = sub %>% tibble::column_to_rownames(var="symbol")
    }

    # Calcaulte co-expression matrix: genes x samples and find hub features
    cor_mat = (lapply(unique(meta$Cancer_type), function(c) {
        samples = meta %>% 
            filter(Cancer_type==c) %>% 
            pull(Sample) %>% intersect(colnames(sub))
        cor(t(sub[,samples]), method = "spearman")
    }) %>% Reduce("+", .)) / length(unique(meta$Cancer_type))
    # cor_mat = cor(t(sub), method = "spearman")

    # IMPLEMENT BETA POWER: Calculate Weighted Adjacency Matrix and select optimal power (the lowest power among those R2>0.8)
    # 1. Run your topology fit function (from the previous step)
    topology_fits <- estimate_scale_free_fit(cor_mat, powers = 1:20)
    # 2. Define your target scale-free R2 threshold (0.80 is standard)
    # r2_threshold <- 0.80
    # 3. Filter for powers that meet the threshold and have a negative slope
    valid_powers <- subset(topology_fits, R2 >= r2_threshold & Slope < 0)
    if (nrow(valid_powers) > 0) {
    # Option A: Choose the lowest power that achieved an R2 >= 0.80
    optimal_beta <- min(valid_powers$Power)
    message(paste("Success: Automated selection chose beta =", optimal_beta, "based on R2 >=", r2_threshold))
    } else {
    # Option B: Fallback if no power hits 0.80 (picks the absolute highest R2)
    optimal_beta <- topology_fits$Power[which.max(topology_fits$R2)]
    warning(paste("Warning: No power reached R2 =", r2_threshold, ". Selecting highest available R2 at beta =", optimal_beta))
    }
    # 4. Pass the automated power directly into your network calculation
    kWithin = rowSums(abs(cor_mat)^optimal_beta) - 1

    hub.nr = max(3, 0.01*length(kWithin))
    hubs = names(sort(kWithin, decreasing = TRUE))[1:hub.nr]

    # library(igraph)
    # Filter the cor matrix to only keep those cor >0.6
    # correlation threshold
    thr = 0.6
    cor_mat[is.na(cor_mat)] = 0
    edges = which(
        abs(cor_mat) > thr & lower.tri(cor_mat),
        arr.ind = TRUE
    )
    # You’re extracting the pairs of genes whose absolute correlation is above your threshold, 
    # but only from the lower triangle to avoid duplicates, and returning their positions as row/column indices. 
    # These pairs represent the edges in your co-expression network.
    edge_df = data.frame(
    from   = rownames(cor_mat)[edges[,1]],
    to     = colnames(cor_mat)[edges[,2]],
    weight = cor_mat[edges]
    )
    # Build a filtered graph
    g = igraph::graph_from_data_frame(edge_df, directed = FALSE)

    # add node attributes
    V(g)$kWithin = kWithin[V(g)$name] #hubness
    V(g)$isHub   = V(g)$name %in% hubs
    # # You’re making a smaller graph containing only your hub genes plus their direct neighbors
    # g_viz = induced_subgraph(
    #   g,
    #   vids = unique(unlist(ego(g, nodes = hubs, order = 1)))
    # )

    set.seed(42)
    paste0(outdir, "/Network/") %>% dir.create(recursive=T, showWarnings=F)
    h = w = size.index*length(V(g))/10 + 5
    pdf(paste0(outdir, "/Network/Network_hub_", gsub(" ", "", System),suffix,"_power",optimal_beta, ".pdf"), h=h,w=w)
    # library(scales)
    edge.color = ifelse(E(g)$weight > 0, "snow2", "skyblue")
    E(g)$weight = abs(E(g)$weight)
    plot(
        g,
        layout = layout_with_fr,
        
        vertex.size = vertex.size.index * (V(g)$kWithin - min(V(g)$kWithin)) / (max(V(g)$kWithin) - min(V(g)$kWithin)),  
        # normalize size nicely
        vertex.color = ifelse(V(g)$isHub, scales::alpha("#D73027",0.6), scales::alpha("orange", 0.6)),  # nicer red and blue shades
        vertex.frame.color = "gold",   # subtle border around nodes
        
        vertex.label = V(g)$name,
        vertex.label.cex = 1,
        vertex.label.color = "black",
        vertex.label.family = "Helvetica",
        
        edge.width = 1 + 4 * abs(E(g)$weight),  # base width plus scaled by weight
        edge.color = edge.color,  # clean blue/red colors
        edge.curved = 0.1,  # slight curve to edges for better visibility
        
        margin = c(0, 0, 0, 0),  # reduce plot margin
        main = "Co-expression Network: Hubs and Neighbors"
    )
    dev.off()

    # # Only keep postively correlated edges
    # edge_df <- subset(edge_df, weight > thr)
    # V(g)$betweenness <- betweenness(g)
    # V(g)$eigen <- eigen_centrality(g)$vector
    # library(ggraph)
    # set.seed(1)
    # ggraph(g_viz, layout = "fr") +
    #   geom_edge_link(aes(width = abs(weight), color = weight > 0), alpha = 0.7) +
    #   geom_node_point(aes(size = kWithin, color = isHub)) +
    #   geom_node_text(aes(label = name), repel = TRUE) +
    #   theme_void()
    return(kWithin)
}





#' @title Survival Analysis via Cox Proportional Hazards Regression
#'
#' @description
#' Performs univariate and multivariable Cox proportional hazards regression to assess 
#' the association between metabolic pathway activity (GSVA scores) and patient survival. 
#' For each pathway, the function fits Cox models with progressively simpler covariate 
#' structures (from full model including stage to univariate models) and tests the 
#' proportional hazards assumption. Results include hazard ratios, confidence intervals, 
#' p-values, and adjusted p-values.
#'
#' @param meta A data frame containing patient metadata with required columns:
#'   \itemize{
#'     \item `sample.submitter_id`: Sample identifiers (must match column names in `GSVAscores`)
#'     \item `vital_status`: Patient status ("Dead" or alive; used to code event indicator)
#'     \item `days_to_death`: Days to death (numeric, may contain NA)
#'     \item `days_to_last_follow_up`: Days to last follow-up (numeric, may contain NA)
#'     \item `age_at_diagnosis`: Patient age at diagnosis (optional, used for covariate adjustment)
#'     \item `gender`: Patient gender (optional, used for covariate adjustment)
#'     \item `ajcc_pathologic_stage`: AJCC cancer stage (optional, used in full model)
#'   }
#'   Rows with missing survival time data are automatically excluded.
#'
#' @param GSVAscores A data frame or matrix containing GSVA pathway activity scores with:
#'   \itemize{
#'     \item First column named "Features": Pathway/feature identifiers
#'     \item Remaining columns: Sample scores with column names matching `meta$sample.submitter_id`
#'       (note: hyphens may be converted to periods in column names)
#'   }
#'   Each row represents one pathway; columns represent samples/patients.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Filters metadata to include only patients with survival time data
#'   \item Matches samples between metadata and GSVA scores
#'   \item Creates survival object (Surv) from vital status and follow-up times
#'   \item For each pathway:
#'     \itemize{
#'       \item Attempts to fit a full Cox model with pathway score + age + gender + stage
#'       \item Falls back to progressively simpler models if fitting fails:
#'         (pathway score + age + gender) → (pathway score + age) → (pathway score only)
#'       \item Tests proportional hazards assumption using Schoenfeld residuals
#'       \item Extracts hazard ratio (HR), 95% CI, and p-values
#'     }
#'   \item Applies false discovery rate (FDR) correction to p-values across all pathways
#'   \item Returns consolidated results table sorted by pathway
#' }
#' Survival time is defined as days to death if available, otherwise days to last follow-up.
#' Event indicator is 1 for deceased patients, 0 for censored (alive) patients.
#'
#' @importFrom survival Surv coxph cox.zph
#' @importFrom dplyr mutate
#'
#' @return A data frame with one row per pathway containing:
#'   \itemize{
#'     \item `Feature`: Pathway identifier
#'     \item `HR`: Hazard ratio (risk of event per unit increase in pathway score)
#'     \item `CI_lower`: Lower bound of 95% confidence interval for HR
#'     \item `CI_upper`: Upper bound of 95% confidence interval for HR
#'     \item `p_value`: Unadjusted p-value from Cox model Wald test
#'     \item `SchoenfieldResidualTest_pvalue`: p-value from proportional hazards assumption test
#'     \item `padj`: False discovery rate (FDR)-adjusted p-value
#'   }
#'   HR > 1 indicates increased risk with higher pathway activity; HR < 1 indicates protective effect.
#'   Rows with failed model fitting are included with NULL values for model parameters.
#'
#' @examples
#' \dontrun{
#'   # Load TCGA metadata and GSVA results
#'   meta <- readRDS("TCGA_metadata.RDS")
#'   gsva_scores <- readRDS("TCGA_GSVA_scores.RDS")
#'   
#'   # Perform Cox regression analysis
#'   survival_results <- metabolicSurvival_coxreg(meta, gsva_scores)
#'   
#'   # Filter for significant pathways (FDR < 0.05)
#'   sig_pathways <- survival_results[survival_results$padj < 0.05, ]
#'   
#'   # Visualize results
#'   head(survival_results)
#' }
#'
#' @export
#' 
metabolicSurvival_coxreg = function(meta, GSVAscores) {
    meta = meta[!(is.na(meta$days_to_death)&is.na(meta$days_to_last_follow_up)),]
    meta = meta[!is.na(meta$sample.submitter_id), ]
    shared = intersect(meta$sample.submitter_id, colnames(GSVAscores))
    meta = meta[meta$sample.submitter_id %in% shared, ]
    Gscores = GSVAscores[, c("Features",  gsub("-", ".", meta$sample.submitter_id))]
    # Create survival time: death time if dead, otherwise last follow up time
    survival_time = 
        ifelse(!is.na(meta$days_to_death), meta$days_to_death, meta$days_to_last_follow_up)
    # Event: 1 if dead, 0 if alive
    status = ifelse(meta$vital_status == "Dead", 1, 0)
    # Assuming:
    # - surv_time: numeric vector of survival times (time to event or censor)
    # - status: event indicator (1=event, 0=censored)
    # - score: continuous predictor vector
    surv_obj = survival::Surv(time = survival_time, event = status)
    return(lapply(1:nrow(Gscores), function(i) {
        # print(i)
        feature = Gscores$Features[i]
        # Function to safely try fitting a Cox model
        try_cox = function(formula) {
            tryCatch({
                cox_model = survival::coxph(formula)
                ph_test = survival::cox.zph(cox_model) # assumption test
                list(cox_model = cox_model, ph_test=ph_test)
            },
                error = function(e) NULL
            )
        }

        # Try full model first
        fit = try_cox(surv_obj ~ 
            as.numeric(Gscores[i,2:ncol(Gscores)])+ meta$age_at_diagnosis+meta$gender+meta$ajcc_pathologic_stage)

        # If error (fit is NULL), try simpler model
        if (is.null(fit)) {
        fit = try_cox(surv_obj ~ 
            as.numeric(Gscores[i,2:ncol(Gscores)])+ meta$age_at_diagnosis+meta$gender)
        }

        if (is.null(fit)) {
        fit = try_cox(surv_obj ~ 
            as.numeric(Gscores[i,2:ncol(Gscores)])+ meta$age_at_diagnosis)
        }

        if (is.null(fit)) {
        fit = try_cox(surv_obj ~ 
            as.numeric(Gscores[i,2:ncol(Gscores)]))
        }

        # Check result
        if (is.null(fit)) {
            # stop("All models failed to fit.")
            summary_table = as.data.frame(list(
                Feature = NULL,
                survival_table = NULL,
                test_statistic = NULL,
                degrees_freedom = NULL,
                p_value = NULL
            ))
        } else {
            sum_model = summary(fit$cox_model)
            # ph_test = survival::cox.zph(cox_model), # assumption test
            # sum_model = summary(cox_model)
            data.frame(
                Feature = feature,
                HR = exp(sum_model$coefficients[1]),
                CI_lower = sum_model$conf.int[,"lower .95"][1],
                CI_upper = sum_model$conf.int[,"upper .95"][1],
                p_value = sum_model$coefficients[,"Pr(>|z|)"][1],
                SchoenfeldResidualTest_pvalue = fit$ph_test$table[1, "p"]
            )
        }
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>%
    dplyr::mutate(padj = p.adjust(p_value, method="fdr")))
}




#' @title Kaplan-Meier Survival Analysis with Log-Rank Test
#'
#' @description
#' Performs Kaplan-Meier survival analysis stratifying patients by metabolic pathway 
#' activity levels. Pathways are dichotomized into "High" and "Low" activity groups 
#' using the 35th and 65th percentiles of GSVA scores (middle 30% excluded). The function 
#' generates survival curves for each pathway, performs log-rank tests to assess 
#' differences between groups, and saves PDF plots. Results include test statistics, 
#' observed vs. expected events, and FDR-adjusted p-values.
#'
#' @param meta A data frame containing patient metadata with required columns:
#'   \itemize{
#'     \item `sample.submitter_id`: Sample identifiers (must match column names in `GSVAscores`)
#'     \item `vital_status`: Patient status ("Dead" or alive; used to code event indicator)
#'     \item `days_to_death`: Days to death (numeric, may contain NA)
#'     \item `days_to_last_follow_up`: Days to last follow-up (numeric, may contain NA)
#'   }
#'   Rows with missing survival time data and missing sample IDs are automatically excluded.
#'
#' @param GSVAscores A data frame or matrix containing GSVA pathway activity scores with:
#'   \itemize{
#'     \item First column named "Features": Pathway/feature identifiers
#'     \item Remaining columns: Sample scores with column names matching `meta$sample.submitter_id`
#'       (note: hyphens may be converted to periods in column names)
#'   }
#'   Each row represents one pathway; columns represent samples/patients.
#'
#' @param outdir Character string specifying the output directory where Kaplan-Meier 
#'   curve PDFs will be saved. Must be a valid, writable directory path.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Filters metadata to include only patients with survival time data and matching sample IDs
#'   \item Creates survival object (Surv) from vital status and follow-up times
#'   \item For each pathway:
#'     \itemize{
#'       \item Dichotomizes pathway scores into "High" (≥65th percentile) and "Low" 
#'         (≤35th percentile), excluding middle 30%
#'       \item Fits Kaplan-Meier curves separately for each group
#'       \item Performs log-rank test (survdiff) to compare survival between groups
#'       \item Generates and saves survival curve plot with group-specific colors
#'         (red for High, blue for Low activity)
#'       \item Extracts test statistics: chi-square, degrees of freedom, p-value
#'     }
#'   \item Applies false discovery rate (FDR) correction to p-values across all pathways
#'   \item Returns summary table with test results for all pathways
#' }
#' Survival time is defined as days to death if available, otherwise days to last follow-up.
#' Event indicator is 1 for deceased patients, 0 for censored (alive) patients.
#' 
#' Two loops perform identical operations: the first generates plots only without storing 
#' results, the second collects results table output.
#'
#' @importFrom survival Surv survfit survdiff
#' @importFrom dplyr mutate
#'
#' @return A data frame with one row per pathway containing:
#'   \itemize{
#'     \item `Feature`: Pathway identifier
#'     \item `survival_table`: Nested data frame with columns:
#'       \itemize{
#'         \item `group`: Group name ("High" or "Low")
#'         \item `N`: Number of subjects per group
#'         \item `Observed`: Observed number of events
#'         \item `Expected`: Expected number of events under null hypothesis
#'       }
#'     \item `test_statistic`: Chi-square test statistic from log-rank test
#'     \item `degrees_freedom`: Degrees of freedom (typically 1 for 2-group comparison)
#'     \item `p_value`: Unadjusted p-value from log-rank test
#'     \item `padj`: False discovery rate (FDR)-adjusted p-value
#'   }
#'   p-value < 0.05 indicates significant difference in survival between High and Low 
#'   activity groups. Significant pathways may be prognostic biomarkers.
#'
#' @section Output Files:
#' PDF files are saved in `outdir` with naming pattern: 
#' "KaplanMeierCurves_<pathway_name>.pdf" where spaces and slashes are removed from 
#' pathway names.
#'
#' @examples
#' \dontrun{
#'   # Load TCGA metadata and GSVA results
#'   meta <- readRDS("TCGA_metadata.RDS")
#'   gsva_scores <- readRDS("TCGA_GSVA_scores.RDS")
#'   
#'   # Perform Kaplan-Meier survival analysis
#'   km_results <- metabolicSurvival_Kaplan_Meier(meta, gsva_scores, 
#'                                                 outdir = "./plots/KaplanMeier/")
#'   
#'   # Filter for significant pathways (FDR < 0.05)
#'   sig_pathways <- km_results[km_results$padj < 0.05, ]
#'   
#'   head(km_results)
#' }
#'
#' @export
#' 
metabolicSurvival_Kaplan_Meier = function(meta, GSVAscores, outdir) {
    meta = meta[!(is.na(meta$days_to_death)&is.na(meta$days_to_last_follow_up)),]
    meta = meta[!is.na(meta$sample.submitter_id), ]
    shared = intersect(meta$sample.submitter_id, colnames(GSVAscores))
    meta = meta[meta$sample.submitter_id %in% shared, ]
    Gscores = GSVAscores[, c("Features",  gsub("-", ".", meta$sample.submitter_id))]
    # Create survival time: death time if dead, otherwise last follow up time
    survival_time = 
        ifelse(!is.na(meta$days_to_death), meta$days_to_death, meta$days_to_last_follow_up)
    # Event: 1 if dead, 0 if alive
    status = ifelse(meta$vital_status == "Dead", 1, 0)
    # Assuming:
    # - surv_time: numeric vector of survival times (time to event or censor)
    # - status: event indicator (1=event, 0=censored)
    # - score: continuous predictor vector
    surv_obj = survival::Surv(time = survival_time, event = status)
    lapply(1:nrow(Gscores), function(i) {
            # print(i)
            feature = Gscores$Features[i]
            # Function to safely try fitting a Cox model
            try_cox = function(formula) {
                tryCatch({
                    fitted = survival::survfit(formula)
                    diff_res = survival::survdiff(formula)
                    list(fitted=fitted, diff_res=diff_res)
                },
                    error = function(e) NULL
                )
            }
            quantiles = 
                quantile(Gscores[i,2:ncol(Gscores)], probs = c(0.25, 0.75), na.rm = TRUE)
            group0 = 
                ifelse(Gscores[i,2:ncol(Gscores)] >= quantiles[2], "High",
                ifelse(Gscores[i,2:ncol(Gscores)] <= quantiles[1], "Low", NA))
            group = group0[!is.na(group0)]
            surv_obj_sub = surv_obj[!is.na(group0)]
            meta_sub = meta[!is.na(group0),]
            meta_sub$group=group
            # Try full model first
            # fit = try_cox(formula = surv_obj_sub ~ 
            #     gorup+ meta_sub$age_at_diagnosis+meta_sub$gender+meta_sub$ajcc_pathologic_stage)

            # # If error (fit is NULL), try simpler model
            # if (is.null(fit)) {
            # fit = try_cox(formula = surv_obj_sub ~ 
            #     group+ meta_sub$age_at_diagnosis+meta_sub$gender)
            # }

            # if (is.null(fit)) {
            # fit = try_cox(formula = surv_obj_sub ~ 
            #     group+ meta_sub$age_at_diagnosis)
            # }

            # if (is.null(fit)) {
            fit = try_cox(formula = surv_obj_sub ~ group)
            # }

            # Check result
            if (is.null(fit)) {
                # stop("All models failed to fit.")
                print(paste0("All models failed to fit for ", feature))
            } else {
                fitted = fit$fitted
                diff_res = fit$diff_res
                pdf(paste0(outdir, "/KaplanMeierCurves_",gsub(" |/","",feature),".pdf"))
                print(plot(fitted, col = c("red", "blue"), xlab = "Time", ylab = "Survival Probability"))
                print(legend("topright", legend = levels(factor(meta_sub$group)), col = c("red", "blue"), lty = 1))
                dev.off()
            }
        })

    rslt = lapply(1:nrow(Gscores), function(i) {
        # print(i)
        feature = Gscores$Features[i]
        # Function to safely try fitting a Cox model
        try_cox = function(formula) {
            tryCatch({
                fitted = survival::survfit(formula)
                diff_res = survival::survdiff(formula)
                list(fitted=fitted, diff_res=diff_res)
            },
                error = function(e) NULL
            )
        }
        quantiles = 
            quantile(Gscores[i,2:ncol(Gscores)], probs = c(0.35, 0.65), na.rm = TRUE)
        group0 = 
            ifelse(Gscores[i,2:ncol(Gscores)] >= quantiles[2], "High",
            ifelse(Gscores[i,2:ncol(Gscores)] <= quantiles[1], "Low", NA))
        group = group0[!is.na(group0)]
        surv_obj_sub = surv_obj[!is.na(group0)]
        meta_sub = meta[!is.na(group0),]
        meta_sub$group=group
        # Try full model first
        # fit = try_cox(formula = surv_obj_sub ~ 
        #     gorup+ meta_sub$age_at_diagnosis+meta_sub$gender+meta_sub$ajcc_pathologic_stage)

        # # If error (fit is NULL), try simpler model
        # if (is.null(fit)) {
        # fit = try_cox(formula = surv_obj_sub ~ 
        #     group+ meta_sub$age_at_diagnosis+meta_sub$gender)
        # }

        # if (is.null(fit)) {
        # fit = try_cox(formula = surv_obj_sub ~ 
        #     group+ meta_sub$age_at_diagnosis)
        # }

        # if (is.null(fit)) {
        fit = try_cox(formula = surv_obj_sub ~ group)
        # }

        # Check result
        if (is.null(fit)) {
            # stop("All models failed to fit.")
            summary_table = as.data.frame(list(
                Feature = NULL,
                survival_table = NULL,
                test_statistic = NULL,
                degrees_freedom = NULL,
                p_value = NULL
            ))
        } else {
            fitted = fit$fitted
            diff_res = fit$diff_res
            # Extract table of groups
            surv_table = data.frame(
                group = names(diff_res$n),
                N = diff_res$n,
                Observed = diff_res$obs,
                Expected = diff_res$exp
            )
            # Extract test statistics
            chisq = diff_res$chisq
            df = length(diff_res$n) - 1
            p_value = 1 - pchisq(chisq, df)

            # Combine everything for reporting
            summary_table = as.data.frame(list(
                Feature = feature,
                survival_table = surv_table,
                test_statistic = chisq,
                degrees_freedom = df,
                p_value = p_value
            ))
        }
        
        # print(plot(fitted, col = c("red", "blue"), xlab = "Time", ylab = "Survival Probability"))
        # print(legend("topright", legend = levels(factor(meta_sub$group)), col = c("red", "blue"), lty = 1))

        summary_table
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>%
        dplyr::mutate(padj = p.adjust(p_value, method="fdr"))
    
    return(rslt)
}



#' @title Create Heatmap of Hazard Ratios from Survival Analysis
#'
#' @description
#' Generates a comprehensive heatmap visualizing hazard ratios (HR) and their statistical 
#' significance from Cox proportional hazards or Kaplan-Meier survival analysis across 
#' multiple cancer types. Pathways are annotated by metabolic class and ordered by total 
#' HR across cohorts. Significance is marked with asterisks (default: FDR < 0.1). 
#' HR values are capped at ±3 for improved visualization. Useful for identifying 
#' metabolic pathways with consistent or divergent prognostic associations.
#'
#' @param metabolicSur_rslts A named list of survival analysis results, with one element 
#'   per cancer type. Each element should be a data frame containing:
#'   \itemize{
#'     \item `Feature`: Pathway identifier (gene set name)
#'     \item `HR`: Hazard ratio from Cox model or Kaplan-Meier analysis
#'     \item `p_value`: Unadjusted p-value
#'     \item `padj`: Adjusted p-value (typically FDR-corrected)
#'   }
#'   List names become cohort/cancer type labels displayed on the heatmap x-axis.
#'   Only pathways with p_value < 0.05 are included in visualization.
#'
#' @param kegg_metab_db_table A data frame containing pathway metadata with at least 
#'   two columns:
#'   \itemize{
#'     \item `gs_name`: Pathway identifier (must match row names in survival results)
#'     \item `class`: Metabolic class/category for row annotations
#'   }
#'   Used to annotate rows in the heatmap with metabolic classification and order pathways.
#'
#' @param outDir Character string specifying the output directory where the PDF heatmap 
#'   will be saved. Must be a valid, writable directory path.
#'
#' @param w Numeric or NULL. Width of output PDF in inches. When NULL, automatically 
#'   computed based on number of cancer types and pathway names. Default: NULL.
#'
#' @param h Numeric or NULL. Height of output PDF in inches. When NULL, automatically 
#'   computed based on number of pathways and class label length. Default: NULL.
#'
#' @param pvalue_cutoff Numeric. Threshold for marking significant pathways with asterisks 
#'   (uses padj column). Default: 0.1.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Filters survival results to include only pathways with p_value < 0.05
#'   \item Collects HR values across all cancer types into a wide-format matrix
#'   \item Collects adjusted p-values in a parallel matrix for significance testing
#'   \item Sorts pathways by total HR across cancer types (descending)
#'   \item Retrieves pathway metadata and orders by metabolic class
#'   \item Creates significance indicator matrix (asterisks for padj < pvalue_cutoff)
#'   \item Clips HR values to ±3 range (winsorization) for visualization
#'   \item Generates heatmap with:
#'     \itemize{
#'       \item Color gradient: blue (HR < 1, protective) → white (HR = 1, neutral) → red (HR > 1, risk)
#'       \item Asterisks overlaid for significant pathways
#'       \item Row annotations showing metabolic class/pathway category
#'       \item No row clustering (maintained metabolic class ordering)
#'       \item Dynamic dimension scaling based on data size
#'     }
#'   \item Saves PDF output to `outDir`
#' }
#' HR > 1 indicates increased risk of death with higher pathway activity (hazard).
#' HR < 1 indicates decreased risk (protective effect).
#' Pathways consistently above/below HR = 1 across multiple cohorts suggest robust 
#' prognostic biomarkers.
#'
#' @importFrom dplyr select full_join filter arrange
#' @importFrom tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#'
#' @return Invisibly returns NULL. Primary side effect is creation of a PDF file 
#'   in `outDir` named "HR_survival.pdf" containing the hazard ratio heatmap.
#'
#' @examples
#' \dontrun{
#'   # Load survival analysis results
#'   cox_results <- readRDS("Cox_regression_results.RDS")  # list of data frames by cancer type
#'   kegg_db <- readRDS("KEGG_metabolic_database.RDS")
#'   
#'   # Create hazard ratio heatmap
#'   visualizeDiffFeatures_sur(metabolicSur_rslts = cox_results,
#'                             kegg_metab_db_table = kegg_db,
#'                             outDir = "./plots/",
#'                             pvalue_cutoff = 0.05)
#'   
#'   # With custom dimensions
#'   visualizeDiffFeatures_sur(metabolicSur_rslts = cox_results,
#'                             kegg_metab_db_table = kegg_db,
#'                             outDir = "./plots/",
#'                             w = 8, h = 12,
#'                             pvalue_cutoff = 0.1)
#' }
#' 
#' @export
#' 
visualizeDiffFeatures_sur = function(
    metabolicSur_rslts=metabolicSur_rslts, 
    kegg_metab_db_table = kegg_metab_dbs$kegg_metab_db_table,
    outDir, w=NULL, h=NULL,
    pvalue_cutoff=0.1) {
    metabolicSur_rslts_dat = lapply(names(metabolicSur_rslts), function(c) {
        cancer = metabolicSur_rslts[[c]]
        cancer = cancer[cancer$p_value <0.05, ]
        # cancer = cancer[cancer$padj <0.05, ]
        if (nrow(cancer)>0) {
            cancer$Cancer_type = c
        }
        cancer
    }) %>% Reduce(rbind, .) %>% arrange(-HR)

    logFCdat = lapply(names(metabolicSur_rslts), function(l) {
        lr = metabolicSur_rslts[[l]] %>% dplyr::select(Feature, HR) %>%
            setNames(c("Feature", l))
    }) %>% Reduce(dplyr::full_join, .) %>% as.data.frame() %>%
        tibble::column_to_rownames(var="Feature")
    logFCdat = logFCdat[rev(order(rowSums(logFCdat))), ]

    padj_dat = lapply(names(metabolicSur_rslts), function(l) {
        # lr = metabolicSur_rslts[[l]] %>% dplyr::select(Feature, p_value) %>%
        lr = metabolicSur_rslts[[l]] %>% dplyr::select(Feature, padj) %>%
            setNames(c("Feature", l))
    }) %>% Reduce(dplyr::full_join, .) %>% as.data.frame() %>%
        tibble::column_to_rownames(var="Feature")
    padj_dat = padj_dat[rownames(logFCdat), ]

    sig_mat = matrix("", nrow=nrow(logFCdat), ncol=ncol(logFCdat))
    sig_mat[padj_dat<pvalue_cutoff] = "*"

    ann_row = unique(kegg_metab_db_table[, c("gs_name", "class")]) %>%
        filter(gs_name %in% rownames(padj_dat)) %>%
        arrange(class) %>%
        tibble::column_to_rownames(var="gs_name")
    logFCdat = logFCdat[rownames(ann_row), ,drop=FALSE]
    padj_dat = padj_dat[rownames(ann_row), ,drop=FALSE]

    # ann_colors = 
    #     replicate(ncol(ann_row), c("0" = "snow2", "1" = "#d95f02"), 
    #     simplify = FALSE) %>% setNames(colnames(ann_row))
    
    logFCdat[abs(logFCdat)>3] = 3
    # max_abs = max(abs(logFCdat), na.rm = TRUE)
    # # Define breaks centered at 0
    breaks = seq(-1, 3, length.out = 101)
    p1 = pheatmap::pheatmap(logFCdat,
            color = colorRampPalette(c("blue", "white", "red"))(100),
            breaks = breaks,
            # annotation_col = ann_row,   # top bar
            display_numbers = sig_mat,
            number_fontsize = 50,
            annotation_row = ann_row,
            cluster_rows = FALSE,
            # annotation_colors = ann_colors,
            show_rownames=TRUE, 
            show_colnames=TRUE,
            fontsize = 10,
            border_color = "white",
            annotation_legend = TRUE,
            main=paste0("HR (*: padj<",pvalue_cutoff,")"))
    if (is.null(w)) {
        w = 0.165*ncol(logFCdat)+0.05*max(nchar(rownames(logFCdat)))+5
    }
    if (is.null(h)) {
        h = 0.12*nrow(logFCdat)+0.05*max(nchar(colnames(ann_row))) +2.5
    }
    pdf(paste0(outDir, "/HR_survival.pdf"), h=h, w=w)
    print(p1)
    dev.off()
}



# Functions for anlayzing SMTdb
#' Import SMTdb files
#'
#' Read SMTdb text files organized under cancer-type/sample directories and
#' return a nested list: top-level keys are cancer types, second-level keys
#' are sample IDs, values are data.frames (or NA when a file cannot be read).
#'
#' @title Import SMTdb files
#' @param INDIR Character. Directory containing SMTdb files organized as
#'   <Cancer_type>/<Sample_ID>/*.txt.
#' @param pattern Character. File name (or pattern) to import, e.g. "st_neb.txt".
#' @param meta Data.frame. Metadata containing at least columns `cancer` and
#'   `datasetID`. Optionally `sample_type` when selecting tumor-only samples.
#' @return A named list (by cancer) of named lists (by sample). Each element
#'   is a data.frame read from the corresponding file or `NA` if missing or
#'   unreadable.
#' @examples
#' # importFiles_SMTdb("/path/to/SMTdb", pattern = "st_gene_exp_count.txt", meta = meta_df)
#' @export
#' @importFrom dplyr filter pull %>%
#' 
importFiles_SMTdb = function(INDIR, pattern="st_gene_exp_count.txt", meta = NULL) {
    if (missing(INDIR) || !dir.exists(INDIR)) stop("INDIR must be an existing directory")
    if (is.null(meta) || !all(c("cancer","datasetID") %in% colnames(meta))) {
        stop("meta must be provided and contain columns 'cancer' and 'datasetID'")
    }
    cancers = meta$cancer %>% unique()
    # cancers = list.dirs(INDIR, full.names=FALSE, recursive=FALSE)
    lapply(cancers, function(c) {
        cancer_dir = paste0(INDIR, "/", c, "/")
        # file_paths = 
        #     list.files(cancer_dir, 
        #     pattern=pattern, full.names=TRUE, recursive=TRUE)
        # file_samples = list.dirs(cancer_dir, full.names=F, recursive=FALSE)
        if (pattern=="st_neb.txt") {
            file_samples = 
                meta %>% 
                filter(cancer==c & slice_type=="cancer") %>%
                pull(datasetID) %>% unique()
        } else {
            file_samples = 
                meta %>% filter(cancer==c) %>% pull(datasetID) %>% unique()
        }
        file_paths = paste0(INDIR, "/", c, "/", file_samples, "/", pattern)

        lapply(seq_along(file_paths), function(i) {
            file = file_paths[i]
            if (!file.exists(file)) {
                warning("Missing file: ", file)
                return(NA)
            }
            tryCatch(
                read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL) %>% as.data.frame(),
                error = function(e) {
                    warning("Failed to read file: ", file, " -> ", conditionMessage(e))
                    NA
                }
            )
        }) %>% setNames(file_samples)
    }) %>% setNames(cancers)
}


#' Make pseudobulk data from spatial transcriptomics data.
#' 
#' Aggregate the gene count across all Malignant/Non-malignant spots for each cancer sample.
#' 
#' @param SpatialNeighborhood The spatial neighborhood data got with importFiles_SMTdb().
#' @param STexpr The spatial transcriptomcis expression data got with importFiles_SMTdb().
#' 
#' @return A list of pseudobulk RNA counts, with each component represent a cancer type.
#' 
#' @examples 
#' # getPseudoBulk(SpatialNeighborhood, STexpr)
#' 
#' @export
#' 
getPseudoBulk_malignant = function(SpatialNeighborhood, STexpr) {
    lapply(names(SpatialNeighborhood), function(c) {
        print(c)
        sp_nei = SpatialNeighborhood[[c]]
        st_expr = STexpr[[c]]
        lapply(names(sp_nei), function(sample) {
            sp_nei_sample = sp_nei[[sample]] %>% na.omit() 
            st_expr_sample = st_expr[[sample]] %>% na.omit() %>%
                tibble::column_to_rownames(var="row.names")
            st_expr_sample = st_expr_sample[rowSums(st_expr_sample, na.rm=T)>1,]
            st_expr_sample = st_expr_sample[, colSums(st_expr_sample, na.rm=T)>0]
            Malignant_spots = sp_nei_sample %>% 
                filter(Location%in%c("Malignant", "Mal")) %>% pull(row.names) %>%
                gsub("-", ".", .) %>%
                intersect(., colnames(st_expr_sample))
            Non_malignant_spots = sp_nei_sample %>% 
                filter(Location%in%c("Non-malignant", "nMal")) %>% pull(row.names) %>%
                gsub("-", ".", .) %>%
                intersect(., colnames(st_expr_sample))
                # we only keep cancers with at least 3 malignant samples
            if (length(Malignant_spots)>=3) {
                Malignant =  rowSums(st_expr_sample[, Malignant_spots], na.rm=T)} else {
                    Malignant = NA
                }
                # we only keep cancers with at least 3 stromal samples
            if (length(Non_malignant_spots)>=3) {
                Non_malignant =  rowSums(st_expr_sample[, Non_malignant_spots], na.rm=T)} else {
                    Non_malignant = NA
                }
            data.frame(Feature=rownames(st_expr_sample),
                Malignant =  Malignant,
                Non_malignant = Non_malignant
            ) %>% setNames(c("Feature", paste0(c("Malignant", "nonMalignant"), "_", sample)))
        }) %>% Reduce(full_join, .)
    }) %>% setNames(names(SpatialNeighborhood))
}


#' Compute GSVA scores for metabolic pathways across spatial samples
#'
#' Calculates gene set variation analysis (GSVA) scores for KEGG metabolic
#' pathways across pseudo-bulk samples from multiple spatial regions. The
#' function normalizes expression data, performs optional feature filtering,
#' and computes pathway enrichment scores for each sample.
#'
#' @param pseudoBulk_cancer A named list of numeric matrices where each element
#'   represents a spatial region (e.g., tissue slice). Each matrix contains
#'   normalized expression counts with features (genes) as rows and samples
#'   (formatted as "GroupType_SpatialID") as columns.
#' @param Feature.filtering Logical scalar. If \code{TRUE}, features are
#'   filtered to retain only those with counts > \code{min.count} in at least
#'   \code{samplProp2rm} proportion of samples within each group (default: \code{TRUE}).
#' @param min.count Numeric scalar. Minimum count threshold for feature filtering
#'   (default: 5).
#' @param samplProp2rm Numeric scalar. Proportion of samples per group that must
#'   exceed \code{min.count} for feature retention. Given as a fraction between
#'   0 and 1 (default: 0.1).
#' @param kegg_metab_db A list of KEGG metabolic pathways mapping pathway
#'   names/IDs to vectors of member gene symbols.
#' @param OUTDIR Character scalar. Directory where output files will be written.
#'   Results are saved in subdirectories named after each spatial region.
#'
#' @return A named list where names correspond to the input list names
#'   (spatial regions). Each element is a data frame of GSVA scores with
#'   metabolic pathways as rows and samples as columns. Regions with fewer than
#'   3 samples per group are excluded (returns \code{NA} and are filtered out).
#'
#' @details
#' For each spatial region in the input list, the function:
#' \enumerate{
#'   \item Parses sample group and slice information from column names
#'   \item Optionally filters features based on expression thresholds
#'   \item Performs TMM normalization using edgeR
#'   \item Saves normalized counts to CSV
#'   \item Computes GSVA pathway enrichment scores using GSVA::gsva()
#'   \item Saves GSVA scores to CSV
#' }
#'
#' Normalized counts and GSVA scores are written to:
#' - \code{OUTDIR/GSVAmetabolicScores/{region}/NormalizedCount.csv}
#' - \code{OUTDIR/GSVAmetabolicScores/{region}/GSVAmetabolicScores.csv}
#'
#' @importFrom dplyr filter mutate across
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom GSVA gsvaParam gsva
#'
#' @export
#' 
getGSVAscores = function(
    pseudoBulk_cancer, Feature.filtering=T, 
    min.count = 5, samplProp2rm = 0.1, kegg_metab_db, OUTDIR) {
    lapply(names(pseudoBulk_cancer), function(c) {
        dat = pseudoBulk_cancer[[c]] 
        rownames(dat) = NULL
        dat = dat %>% as.data.frame() %>%
            tibble::column_to_rownames(var="Feature") #%>% na.omit()
        meta0 = data.frame(Slice = colnames(dat) %>% 
                    gsub("nonMalignant_", "", .) %>%
                    gsub("Malignant_", "", .),
                Group = gsub("_.*", "", colnames(dat)),
                Sample = colnames(dat))
        # We need to make sure there are more than 2 samples in each condition:
        if (sum(table(meta0$Group)>2)==2) {
            meta0$Group = factor(meta0$Group, levels=c("nonMalignant", "Malignant")) # Reference first!
            paste0(OUTDIR, "/GSVAmetabolicScores/", c, "/") %>%
                dir.create(., recursive=TRUE, showWarnings=FALSE)
                meta0 = filter(meta0, Sample%in%colnames(dat))

            # Make sure the sample are in the same order in both meta and dat:
            dat = dat[, meta0$Sample]
            meta0 = meta0 %>% mutate(across(where(is.factor), droplevels))
            mapping = setNames(meta0$Group, meta0$Sample)
            print("Only keep features that have at least min.nr number of values that are >5 in each group...")
            min.nr = round(max((ncol(dat)/2)*samplProp2rm, 2), 0)
            if (Feature.filtering) {
                B = mapping[colnames(dat)] %>% as.character()
                featureskeep = apply(dat>min.count, 1, function(row) {
                    A = row
                    sum(aggregate(A~B, data=data.frame(A=A,B=B), sum, na.rm=T)[,"A"]>=min.nr)>1
                })
                dat = dat[featureskeep, ]
            }
            temp = replace(dat, is.na(dat), 0)
            dge = edgeR::DGEList(counts=temp)
            dge = edgeR::calcNormFactors(dge, method = "TMMwsp")
            effect.lib.sizes = (dge$samples$lib.size) * (dge$samples$norm.factors)
            normalizedCount = (t(t(dat)/effect.lib.sizes))*median(effect.lib.sizes)
            write.table(
                normalizedCount %>% as.data.frame() %>%
                        tibble::rownames_to_column(var="Features"),
                paste0(OUTDIR, "/GSVAmetabolicScores/", c, "/NormalizedCount.csv"),
                quote = F, row.names = F, col.names = T, sep="\t")

            # Calcualte GSVA per sample
            # BiocManager::install("GSVA")
            normalizedCount = (t(t(temp)/effect.lib.sizes))*median(effect.lib.sizes)
            normalizedCount = as.matrix(normalizedCount)
            normalizedCount = matrix(as.numeric(normalizedCount),
                nrow = nrow(normalizedCount),
                ncol = ncol(normalizedCount), dimnames = dimnames(normalizedCount))

            param = GSVA::gsvaParam(normalizedCount, kegg_metab_db, minSize=5, maxSize=500)
            gsva_scores = 
                tryCatch(GSVA::gsva(param, verbose=TRUE) %>% 
                as.data.frame(), error = function(e) NA)
            write.table(
                gsva_scores %>% as.data.frame() %>%
                        tibble::rownames_to_column(var="Features"),
                paste0(OUTDIR, "/GSVAmetabolicScores/", c, "/GSVAmetabolicScores.csv"),
                quote = F, row.names = F, col.names = T, sep="\t")
            return(as.data.frame(gsva_scores))
        } else {return(NA)}
    }) %>% setNames(names(pseudoBulk_cancer)) %>% .[!is.na(.)]
} 


#' Differential analysis of GSVA metabolic pathway scores using limma
#'
#' Performs differential testing on GSVA-derived metabolic pathway scores
#' between sample groups using limma. Supports optional sample pairing and
#' additional covariates in the linear model.
#'
#' @param GSVAscores A named list of data frames where each element represents
#'   a spatial region. Each data frame contains GSVA scores (metabolic pathways
#'   as rows, samples as columns with group identifier prefix as "Group_ID").
#' @param currentCovariate Optional character scalar. Name of an additional
#'   covariate to include in the model alongside \code{Group}. If \code{NULL},
#'   only \code{Group} is used (default: \code{NULL}).
#' @param paired Logical scalar. If \code{TRUE}, a duplicate correlation model
#'   is fitted to account for repeated measures/paired samples. Set to
#'   \code{FALSE} for unpaired analysis (default: \code{TRUE}).
#' @param OUTDIR Character scalar. Directory where output results files
#'   (limma statistics) will be written in subdirectories.
#'
#' @return A named list where names correspond to the input list names
#'   (spatial regions). Each element is a data frame of limma topTable results
#'   sorted by t-statistic (descending), containing:
#'   \item{logFC}{Log2 fold-change}
#'   \item{AveExpr}{Average log2 expression}
#'   \item{t}{t-statistic}
#'   \item{P.Value}{Raw p-value}
#'   \item{adj.P.Val}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{B}{B-statistic (log-odds of differential expression)}
#'
#' @details
#' For each spatial region in the input list, the function:
#' \enumerate{
#'   \item Extracts group and sample information from column names
#'   \item Constructs a design matrix from \code{Group} and optional \code{currentCovariate}
#'   \item Fits a linear model with limma::lmFit()
#'   \item If \code{paired = TRUE}, applies duplicate correlation correction
#'   \item Applies empirical Bayes moderation with trend (limma::eBayes)
#'   \item Extracts differential statistics for the second factor level (typically the test group)
#'   \item Writes results to \code{OUTDIR/Limma_GSVA/{region}/limmaResult.csv}
#' }
#'
#' Sample column names must follow the format "GroupType_SampleID" (e.g.,
#' "Malignant_001", "nonMalignant_002"). The group variable is extracted by
#' removing everything after the first underscore.
#'
#' @importFrom limma duplicateCorrelation lmFit eBayes topTable
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr arrange
#'
#' @export
#' 
GSVAlimmaTest_inner = function(
    GSVAscores, currentCovariate=NULL, paired=T, OUTDIR) {
    lapply(names(GSVAscores), function(c) {
        g = GSVAscores[[c]]
        meta0 = data.frame(# Sample = colnames(g),
                        Group = gsub("_.*", "", colnames(g)),
                        Sample = gsub("nonMalignant_|Malignant_", "", colnames(g)))
        paired_variable = "Sample"
        baseFormula = ~ Group
        if (is.null(currentCovariate)) {
            currentFormula = baseFormula
        } else {
            currentFormula = as.formula(paste0("~ ", currentCovariate, " + Group"))}

        design = model.matrix(currentFormula, data=meta0) #
        colnames(design) = gsub("Group", "", colnames(design))
        if (paired) { #
            #Fit with correlated arrays
            dupcor = limma::duplicateCorrelation(g, design, #
                block=meta0[, paired_variable, drop=T]) #
            fit = limma::lmFit(g, design, block=meta0[,paired_variable,drop=T], #
                correlation=dupcor$consensus)
        } else { fit = limma::lmFit(g, design) }

        # # Limma trend to refine gene-level variance estimates
        # gssizes = sapply(kegg_metab_db[rownames(fit$coefficients)], function(gs) {
        #     length(gs)
        # })
        # fit.con = limma::eBayes(fit, robust = TRUE, trend=gssizes)
        fit.con = limma::eBayes(fit, robust = TRUE, trend=TRUE)
        rlst_interest =
            limma::topTable(fit.con, n=Inf, coef=levels(meta0$Group)[2]) %>%
            arrange(-t)
        paste0(OUTDIR, "/Limma_GSVA/", c, "/") %>% dir.create(., recursive=T, showWarnings=F)
        write.table(
            rlst_interest %>% as.data.frame() %>%
                    tibble::rownames_to_column(var="Features"),
            paste0(OUTDIR, "/Limma_GSVA/", c, "/limmaResult.csv"), #
            quote = F, row.names = F, col.names = T, sep="\t")
        return(rlst_interest)
    }) %>% setNames(names(GSVAscores))
}

