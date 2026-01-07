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
 #'
 #' @return A list with components:
 #' \item{rslt_of_interest}{A data frame of limma topTable results for the
 #'   comparison of interest (including statistics and adjusted p-values).}
 #' \item{gsva}{The matrix of GSVA scores (pathways x samples).}
 #' \item{fitted.model}{The limma fit object (returned invisibly).}
 #'
 GSVAlimmaTest = function(dat, paired = TRUE, Feature.filtering = TRUE,
     currentCovariate = NULL, meta = meta, paired_variable = "Patient",
     OUTDIR, kegg_metab_db) {
    meta = filter(meta, Sample%in%colnames(dat))
    # Make sure the sample are in the same order in both meta and dat:
    dat = dat[, meta$Sample]
    meta = meta %>% mutate(across(where(is.factor), droplevels))
    mapping = setNames(meta$Group, meta$Sample)
    print("Only keep features that have at least min.nr number of values that are >5 in each group...")
    min.nr = round(max((ncol(dat)/2)*0.15, 2), 0)
    if (Feature.filtering) {
        B = mapping[colnames(dat)] %>% as.character()
        featureskeep = apply(dat>5, 1, function(row) {
            A = row
            sum(aggregate(A~B, data=data.frame(A=A,B=B), sum, na.rm=T)[,"A"]>=min.nr)>1
        })
        dat = dat[featureskeep, ]
    }

    # library(limma)
    # conda install conda-forge::r-locfit
    # conda install bioconda::bioconductor-edger
    # BiocManager::install("edgeR")
    # library(edgeR)
    Group = factor(meta$Group)
    # Create a DEGList from the edgeR package
    dge = edgeR::DGEList(counts=dat)
    # Normalize the data
    dge = edgeR::calcNormFactors(dge)
    # Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.
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
    dat2 = limma::voom(dge, design, plot=T)
    # What is voom doing?
    # Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
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
    # library(GSVA)
    # library(tibble)
    # library(dplyr)
    normalizedCount = as.matrix(normalizedCount)
    normalizedCount = matrix(as.numeric(normalizedCount),
        nrow = nrow(normalizedCount),
        ncol = ncol(normalizedCount), dimnames = dimnames(normalizedCount))
    param = GSVA::gsvaParam(normalizedCount, kegg_metab_db, minSize=5, maxSize=500)
    gsva_scores = GSVA::gsva(param, verbose=TRUE)
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
    gssizes = sapply(kegg_metab_db[rownames(fit$coefficients)], function(gs) {
        length(gs)
    })
    fit.con = limma::eBayes(fit, robust = TRUE, trend=gssizes)
    rlst_interest =
        limma::topTable(fit.con, n=Inf, coef=levels(meta$Group)[2]) %>%
        arrange(-t)
    
    return(list(rslt_of_interest=rlst_interest, gsva=gsva_scores, fitted.model=fit))

}


#' @import dplyr
#' @param databaseDIR a diretory to save the KEGG metabolic pathway database
#'
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
        paste0(outdir, "/KEGG_metabolism/KEGG_metabolism.csv"),
        header = TRUE, stringsAsFactors = FALSE) %>%
        mutate(class = gsub("Metabolism; ", "", class))
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
 #' @return Invisibly returns the generated \code{ggplot2} object. The primary
 #'   side-effect is writing the PDF file to \code{outdir} when provided.
 #'
 getCirPackingPlot_GSVA = function(kegg_metab_db_table, outdir) {
    pathway_hierachy =  kegg_metab_db_table %>% mutate(Depth0="Metabolism") %>%
        mutate(Depth2=gs_name) %>% mutate(Depth1=class) %>%
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
            mutate(shortName=gsub("etabolism",".",paste0(name,"(",size,")"))) %>%
            mutate(shortName=gsub("iosynthesis","iosyn.", shortName)) %>%
            # mutate(shortName=gsub("secondary metabolites", "sec. mets.", shortName)) %>%
            mutate(shortName = gsub(" m.", "", shortName)) %>%
            mutate(shortName = gsub(" biosyn. and m.", "", shortName)) %>%
            mutate(shortName = gsub("M. of ", "", shortName)) %>%
            mutate(shortName = gsub(" and m.", "", shortName)) %>%
            mutate(shortName = gsub(" biosyn. and", "", shortName)) %>%
            mutate(shortName = gsub(" biodegrad. and", "", shortName)) %>%
            mutate(shortName = gsub("Biosyn. of ", "", shortName))
    vertices$shortName[!vertices$name%in%unique(pathway_hierachy$Depth1)] = NA
    mygraph = igraph::graph_from_data_frame(edges, vertices=vertices)
    p = ggraph::ggraph(mygraph, layout = 'circlepack', weight=size) +
        ggraph::geom_node_circle(aes(col=depth)) +
        ggraph::geom_node_label(aes(label=shortName),repel=T,size=4,color="darkgreen") +
        ggplot2::theme_void() +
        viridis::scale_color_viridis() +
        ggplot2::theme(legend.position="FALSE") +
        ggplot2::ggtitle("                                           TCGA")

    pdf(paste0(outdir, "/KEGG_metabolism/TaskSummary.pdf"), h=5, w=4.85)
        print(p)
    dev.off()
}


 #' Visualize GSVA pathway means across groups (and normals)
 #'
 #' Generate side-by-side boxplots that summarize GSVA pathway activity means
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
 #'
 #' @return Invisibly returns a list containing the two ggplot2 objects (all
 #'   samples and normal-only samples) and the path to the saved PDF (if
 #'   `OUTDIR` is provided). A PDF named `GSVA_samples.pdf` is written to
 #'   `OUTDIR` when `OUTDIR` is non-NULL.
 #'
 visualizeGSVSscoresGroup = function(GSVA_limma_rslt_gsva, paired_data, OUTDIR) {
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
    P1 = dat %>% mutate(Cancer_type = factor(Cancer_type, levels=Cancer_type_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Cancer_type, y=means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(title = "TCGA",
                x="Cancer_type", y="Mean Metab. Act.")+
            ggplot2::theme_minimal() +
            ggplto2::theme(legend.position="none",
                axis.text.x = element_text(angle = 90, hjust=1, vjust=0))

    # GSVA scores for healthy control samples only
    dat = lapply(seq_along(GSVA_limma_rslt_gsva), function(j) {
        normal_samples = paired_data[[j]] %>%
            filter(sample_type=="Solid Tissue Normal") %>% pull(sample.submitter_id)
        gsva = GSVA_limma_rslt_gsva[[j]] %>% .[, normal_samples]
        data.frame(
            means = (rowSums(gsva, na.rm=TRUE)/ncol(gsva)) %>% as.numeric(),
            Cancer_type = names(GSVA_limma_rslt_gsva)[j])
        }) %>% Reduce(rbind, .)
    Cancer_type_order = dat %>% group_by(Cancer_type) %>%
            summarise(Means=mean(means)) %>%
            arrange(Means) %>% pull(Cancer_type) %>% unique()
    P2 = dat %>% mutate(Cancer_type = factor(Cancer_type, levels=Cancer_type_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Cancer_type, y=means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(title = "TCGA",
                x="Cancer_type", y="Mean Norm. Metab. Act.")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUTDIR, "/GSVA_samples.pdf"), w=5.5, h=3.0)
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
 #' @return A pdf file, Pathway_size_class.pdf in the OUTDIR.
 #' @examples
 #' # visualizeGSsizeClass(GSVA_limma_rslt_gsva, genes4GSVA, kegg_metab_db, kegg_metab_db_table, OUTDIR)
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
        mutate(Class = gsub("etabolism",".", Class)) %>%
        mutate(Class = gsub("thesis", ".", Class)) %>%
        # mutate(Class = gsub("secondary metabolites", "sec. mets.", Class)) %>%
        mutate(Class = gsub("biodegradation", "biodegrad.", Class)) %>%
        mutate(Class = gsub(" m.", "", Class)) %>%
        mutate(Class = gsub(" biosyn. and m.", "", Class)) %>%
        mutate(Class = gsub("M. of ", "", Class)) %>%
        mutate(Class = gsub(" and m.", "", Class)) %>%
        mutate(Class = gsub(" biosyn. and", "", Class)) %>%
        mutate(Class = gsub(" biodegrad. and", "", Class)) %>%
        mutate(Class = gsub("Biosyn. of ", "", Class)) %>%
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

 #' Visualize GSVA pathway means by metabolic class
 #'
 #' Produce comparative boxplots of mean GSVA pathway activities grouped by
 #' metabolic class. The function creates two panels: overall cohort means and
 #' means computed using matched normal samples only (when available in
 #' `paired_data`). Outputs are saved to `OUTDIR` when provided.
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
 #'
 #' @return Invisibly returns a list with the ggplot2 objects for the two
 #'   panels (overall and normals-only) and, when \code{OUTDIR} is provided,
 #'   the path to the saved PDF file. The primary side-effect is creation of a
 #'   PDF file in \code{OUTDIR} named \file{GSVA_pathway_class.pdf}.
 #'
 visualizeGSVSscoresClass = function(GSVA_limma_rslt, kegg_metab_db_table, paired_data, OUTDIR) {
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
                    mutate(Class = gsub("etabolism",".", Class)) %>%
                    mutate(Class = gsub("thesis", ".", Class)) %>%
                    # mutate(Class = gsub("secondary metabolites", "sec. mets.", Class)) %>%
                    mutate(Class = gsub("biodegradation", "biodegrad.", Class)) %>%
                    mutate(Class = gsub(" m.", "", Class)) %>%
                    mutate(Class = gsub(" biosyn. and m.", "", Class)) %>%
                    mutate(Class = gsub("M. of ", "", Class)) %>%
                    mutate(Class = gsub(" and m.", "", Class)) %>%
                    mutate(Class = gsub(" biosyn. and", "", Class)) %>%
                    mutate(Class = gsub(" biodegrad. and", "", Class)) %>%
                    mutate(Class = gsub("Biosyn. of ", "", Class))
        }) %>% Reduce(rbind, .)
        Class_order = dat %>% group_by(Class) %>%
            summarise(Means=mean(Means)) %>% arrange(Means) %>% pull(Class) %>% unique()
    P1 = dat %>% mutate(Class = factor(Class, levels=Class_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Class, y=Means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(# title = "TCGA",
                x="", y="Mean Metab. Act.")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    dat = lapply(seq_along(GSVA_limma_rslt), function(j) {
        s = names(GSVA_limma_rslt)[j]
        normal_samples = paired_data[[s]] %>%
            filter(sample_type=="Solid Tissue Normal") %>% pull(sample.submitter_id)
        gsva = GSVA_limma_rslt[[s]]$gsva %>% .[, normal_samples]
        pathway_means = rowSums(gsva,na.rm=TRUE)/ncol(gsva)
        kegg_metab_db_table =
            kegg_metab_dbs$kegg_metab_db_table %>% dplyr::select(gs_name, class) %>%
            unique()
        kegg_metab_db_mapping = setNames(kegg_metab_db_table$class, kegg_metab_db_table$gs_name)
        data.frame(Pathway = names(pathway_means),
                    Means = as.numeric(pathway_means),
                    Class = kegg_metab_db_mapping[names(pathway_means)]) %>%
                    mutate(Class = gsub("etabolism",".",Class)) %>%
                    mutate(Class = gsub("thesis", ".", Class)) %>%
                    # mutate(Class = gsub("secondary metabolites", "sec. mets.", Class)) %>%
                    mutate(Class = gsub("biodegradation", "biodegrad.", Class)) %>%
                    mutate(Class = gsub(" m.", "", Class)) %>%
                    mutate(Class = gsub(" biosyn. and m.", "", Class)) %>%
                    mutate(Class = gsub("M. of ", "", Class)) %>%
                    mutate(Class = gsub(" and m.", "", Class)) %>%
                    mutate(Class = gsub(" biosyn. and", "", Class)) %>%
                    mutate(Class = gsub(" biodegrad. and", "", Class)) %>%
                    mutate(Class = gsub("Biosyn. of ", "", Class))
        }) %>% Reduce(rbind, .)
        Class_order = dat %>% group_by(Class) %>%
            summarise(Means=mean(Means)) %>% arrange(Means) %>% pull(Class) %>% unique()
    P2 = dat %>% mutate(Class = factor(Class, levels=Class_order)) %>%
            ggplot2::ggplot(., ggplot2::aes(x=Class, y=Means)) +
            ggplot2::geom_boxplot(fill="green4") +
            ggplot2::labs(# title = "TCGA",
                x="", y="Mean Norm. Metab.")+
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position="none",
                axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUTDIR, "/GSVA_pathway_class.pdf"), w=5.5, h=3.5)
        cowplot::plot_grid(P1, P2)
    dev.off()
}


 #' Visualize number of significant pathways per group
 #'
 #' Create a bar chart summarizing the number of significantly different
 #' pathways per group (up- and down-regulated) based on limma results
 #' produced by \code{GSVAlimmaTest}. The plot is saved as a PDF in
 #' \code{OUTDIR}.
 #'
 #' @param OUTDIR Character scalar. Directory where the output PDF
 #'   (\file{Sig_pathway_nr.pdf}) will be written.
 #' @param GSVA_limma_rslt A named list of results returned by
 #'   \code{GSVAlimmaTest()} (one element per group). Each element must
 #'   contain a component \code{rslt_of_interest} (a data frame produced by
 #'   \code{limma::topTable}) with at least \code{logFC} and \code{adj.P.Val}
 #'   columns.
 #'
 #' @return Invisibly returns the path to the saved PDF file; the primary
 #'   side-effect is the creation of \file{Sig_pathway_nr.pdf} in \code{OUTDIR}.
 #'
 visualizeSigNr = function(OUTDIR, GSVA_limma_rslt) {
    pdf(paste0(OUTDIR, "/Sig_pathway_nr.pdf"), w=3, h=3.5)
        sig_nr = lapply(names(GSVA_limma_rslt), function(cancer) {
            rslt = GSVA_limma_rslt[[cancer]]$rslt_of_interest %>%
                filter(adj.P.Val<0.05)
            data.frame(up = rslt %>% filter(logFC>0) %>% nrow(),
                down = -(rslt %>% filter(logFC<0) %>% nrow())) %>%
                t() %>% as.data.frame() %>% setNames(cancer)
            }) %>% Reduce(cbind, .) %>%
                tibble::rownames_to_column(var="Regulation") %>%
                tidyr::gather(key="Cancer_type", value="Count", -Regulation) %>%
                arrange(-Count) %>%
                mutate(Cancer_type=factor(Cancer_type,levels=unique(Cancer_type)))
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


#' Dot-plot summary of GSVA metabolic activity-based limma results across cancer types
#'
#' Produce a dot-plot that summarizes pathway-level limma results across
#' multiple cancer cohorts (TCGA). For each cohort (element of
#' \code{rlst_list}) the function selects significant pathways by
#' adjusted p-value, annotates directionality (Up/Down) based on
#' \code{logFC}, joins pathway metadata from \code{taskHierarchy}, and
#' draws a dot for each pathway-cohort pair sized by adjusted p-value and
#' colored by the limma t-statistic. The plot is written to a PDF in
#' \code{OUTDIRV}.
#'
#' @param rlst_list list. A list of data.frames or tibbles produced by
#'   limma (one element per cohort). Each element should contain at
#'   least the columns \code{logFC}, \code{adj.P.Val} and \code{t}.
#' @param adj.P.Val.cutoff numeric. Adjusted p-value threshold used to
#'   select significant pathways (default: user-supplied). Only rows with
#'   \code{adj.P.Val < adj.P.Val.cutoff} are plotted.
#' @param w numeric or NULL. Plot width in inches; a sensible default is
#'   computed when \code{NULL}.
#' @param h numeric or NULL. Plot height in inches; a sensible default is
#'   computed when \code{NULL}.
#' @param OUTDIRV character. Directory where the generated PDF dotplot
#'   will be saved.
#' @param Depth character. Column name in \code{taskHierarchy} to use as
#'   the pathway label (default: \code{"gs_name")).
#' @param lowfigure logical. When \code{TRUE} produce a compact layout
#'   with legends placed below the plot (default: \code{FALSE}).
#' @param title character. Plot title (default: \code{"Dysregulated pathways_TCGA"}).
#' @param suffix character. Filename suffix appended to the saved PDF
#'   (default: \code{"TCGA").
#' @param taskHierarchy data.frame. Mapping table containing pathway
#'   metadata; must include the pathway label column named by
#'   \code{gs_name} and a \code{class} column used for optional nesting.
#' @param nested logical. If \code{TRUE} use nested facets to group
#'   pathways by their class (default: \code{TRUE}).
#' @param top logical. If \code{TRUE} keep only frequently significant
#'   pathways (default: \code{FALSE}).
#' @param color_low character. Color used for negative t-statistics
#'   (default \code{"skyblue").
#' @param color_high character. Color used for positive t-statistics
#'   (default \code{"purple").
#'
#' @return Invisibly returns the ggplot object used to render the
#'   figure. The primary side-effect is the PDF file
#'   \code{Dotplot_all_adj.P.ValCutoff<cutoff>_<suffix>.pdf} written into
#'   \code{OUTDIRV}.
#'
#' @examples
#' \dontrun{
#' makeLimmaDotplot_TCGA(rlst_list, adj.P.Val.cutoff = 0.05, OUTDIRV = "./plots/")
#' }
#' @export

makeLimmaDotplot_TCGA = function(rlst_list, adj.P.Val.cutoff, w=NULL, h=NULL, OUTDIRV, 
    Depth="gs_name", lowfigure=F, title = "Dysregulated pathways_TCGA", 
    suffix="TCGA", taskHierarchy=cfs, nested=T, top=FALSE,
    color_low="skyblue", color_high="purple"){
    taskHierarchy$Task = taskHierarchy[[Depth]]
    rlst_long = lapply(seq_along(rlst_list), function(i) {
        rlst = rlst_list[[i]]
        rlst$Cancer = names(rlst_list)[i]
        rlst$Regulation = 
            ifelse(rlst$logFC>0&rlst$adj.P.Val<adj.P.Val.cutoff, "Up", 
            ifelse (rlst$logFC<0&rlst$adj.P.Val<adj.P.Val.cutoff, "Down", "NonSig"))
        rlst %>% tibble::rownames_to_column(var="Task") %>% 
            filter(adj.P.Val<adj.P.Val.cutoff)
    }) %>% Reduce(rbind, .) %>%
        dplyr::left_join(taskHierarchy) %>%
            mutate(logFC=round(logFC, 2)) %>%
            mutate(Class = gsub("Carbohydrate metabolism","CAR",class)) %>%
            mutate(Class = gsub("Lipid metabolism", "LIP", Class)) %>%
            mutate(Class = gsub("Metabolism of cofactors and vitamins", "COF", Class)) %>%
            mutate(Class = gsub("Energy metabolism", "ENE", Class)) %>%
            mutate(Class = gsub("Amino acid metabolism", "AMI", Class)) %>%
            mutate(Class = gsub("Nucleotide metabolism", "NUC", Class)) %>%
            mutate(Class = gsub("Biosynthesis of other secondary metabolites", "OSM", Class)) %>%
            mutate(Class = gsub("Metabolism of other amino acids", "OAA", Class)) %>%
            mutate(Class = gsub("Glycan biosynthesis and metabolism", "GLY", Class)) %>%
            mutate(Class = gsub("Metabolism of terpenoids and polyketides", "TER", Class)) %>%
            mutate(Class = gsub("Xenobiotics biodegradation and metabolism" , "XEN", Class))

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
        ggplot2::geom_point(ggplot2::aes(size=adj.P.Val, color=t, shape=Regulation))
    if (nested) p = p + ggh4x::facet_nested(Class~., scales = "free_y", space = "free")
    p = p + ggplot2::scale_shape_manual(values = c("Up" = 18, "Down" = 20)) +
        ggplot2::ggtitle(title) +
        ggplot2::ylab("") +  ggplot2::labs(size="adj.P.Val") +
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
        ggplot2::labs(caption=paste0("adj.P.ValCutoff = ",adj.P.Val.cutoff,", pAdjustMethod = BH"))

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
    pdf(paste0(OUTDIRV,"/Dotplot_all_adj.P.ValCutoff",
        adj.P.Val.cutoff,"_",suffix,".pdf"),w=w,h=h)
        show(p)
    dev.off()
}



 #' Differential effect boxplot by metabolic class
 #'
 #' Generate a boxplot summarizing differential effects (t-statistics) for
 #' KEGG metabolic pathways aggregated by metabolic class. The function
 #' collects limma test results across input groups, maps pathways to classes
 #' via the provided hierarchy table, and writes a PDF boxplot to
 #' \code{OUT_DIR}.
 #'
 #' @param rlst_list A named list of results returned by \code{GSVAlimmaTest()}
 #'   (one element per group). Each element must contain a component
 #'   \code{rslt_of_interest} (a data frame produced by \code{limma::topTable})
 #'   with at least the columns \code{t}, \code{logFC} and \code{adj.P.Val}.
 #' @param taskHierarchy A data.frame or tibble with at least two columns:
 #'   \code{gs_name} (pathway name or ID) and \code{class} (the higher-level
 #'   metabolic class). This table is used to map pathways to classes for
 #'   grouping and plotting.
 #' @param OUT_DIR Character scalar. Directory where the output PDF
 #'   (\file{diffEffectBoxplot_system.pdf}) will be written.
 #' @param w Numeric. Output figure width (in inches).
 #' @param h Numeric. Output figure height (in inches).
 #' @param title Character. Plot title.
 #'
 #' @return Invisibly returns the path to the written PDF; primary side-effect
 #'   is creation of \file{diffEffectBoxplot_system.pdf} in \code{OUT_DIR}.
 #'
 diffEffectBoxplot_bySystem_GSVA = function(
     rlst_list, taskHierarchy, OUT_DIR, w = 2, h = 3.5, title = "Zhou2020") {
    B_Z_N0 = lapply(names(rlst_list), function(cancer){
        temp = rlst_list[[cancer]]$rslt_of_interest
        temp$Cancer = cancer
        temp %>% tibble::rownames_to_column(var="gs_name") 
        }) %>% Reduce(rbind, .) %>% as.data.frame() %>% 
        dplyr::left_join(taskHierarchy %>% unique()) %>%
        mutate(class = gsub("etabolism",".",class)) %>%
        mutate(class = gsub("thesis", ".", class)) %>%
        mutate(class = gsub("biodegradation", "biodegrad.", class)) %>%
        mutate(class = gsub(" m.", "", class)) %>%
        mutate(class = gsub(" biosyn. and m.", "", class)) %>%
        mutate(class = gsub("M. of ", "", class)) %>%
        mutate(class = gsub(" and m.", "", class)) %>%
        mutate(class = gsub(" biosyn. and", "", class)) %>%
        mutate(class = gsub(" biodegrad. and", "", class)) %>%
        mutate(class = gsub("Biosyn. of ", "", class))

    System_orderZ = B_Z_N0 %>% unique() %>% group_by(class) %>%
        dplyr::summarise(meanEffect = median(t, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(class) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(class=factor(class, levels = System_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=class, y=t)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="Class", y="Diff effect")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_system.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}


 #' Differential effect boxplot by cancer/group
 #'
 #' Combine limma test outputs from multiple groups and plot the distribution
 #' of differential effects (t-statistics) per group (for example, cancer
 #' types). The function extracts the top-table results provided by
 #' \code{GSVAlimmaTest()}, computes a per-group median effect ordering and
 #' writes a boxplot PDF to the specified output directory.
 #'
 #' @param rlst_list Named list. Results returned by \code{GSVAlimmaTest()}
 #'   (one element per group). Each element must contain a component
 #'   \code{rslt_of_interest}, a data frame (as produced by
 #'   \code{limma::topTable}) that includes at minimum the columns
 #'   \code{t}, \code{logFC} and \code{adj.P.Val}.
 #' @param OUT_DIR character(1). Directory where the output PDF will be
 #'   written. The plot file is named \file{diffEffectBoxplot_cancer.pdf}.
 #' @param w numeric(1). Figure width in inches (default: 2.75).
 #' @param h numeric(1). Figure height in inches (default: 3.0).
 #' @param title character(1). Plot title (default: "Zhou2020").
 #'
 #' @return Invisibly returns the path to the generated PDF; the primary
 #'   side-effect is creation of \file{diffEffectBoxplot_cancer.pdf} in
 #'   \code{OUT_DIR}.
 #'
diffEffectBoxplot_byCancer_GSVA = function(
    rlst_list, OUT_DIR, w=2.75, h=3.0, title="Zhou2020") {
    B_Z_N0 = lapply(names(rlst_list), function(cancer){
        temp = rlst_list[[cancer]]$rslt_of_interest
        temp$Cancer_type = cancer
        temp %>% tibble::rownames_to_column(var="gs_name") 
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>% unique()

    Cancer_orderZ = B_Z_N0 %>% unique() %>% group_by(Cancer_type) %>%
        dplyr::summarise(meanEffect = median(t, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(Cancer_type) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(Cancer=factor(Cancer_type, levels = Cancer_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=Cancer, y=t)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="Cancer", y="Diff effect")+
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
 #' @return A data.frame with one row per sample and category containing the
 #'   aggregated statistics: \code{MeanTaskScore}, \code{MeanActiveTaskScore},
 #'   \code{SumTaskScore}, \code{SumActiveTaskScore}, \code{FractionOfActiveTasks}
 #'   and \code{NrOfActiveTasks}. Invisibly returns the same data.frame; the
 #'   primary side-effect is writing the PDF to \code{OUTDIR}.
 #'
 #' @examples
 #' # depths <- sumUpTaskScores(cfs, Depth = "Depth1", h = 12, w = 30, OUTDIR = "./out")
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
 #' @return Invisibly returns the summarized data frame of gene counts per
 #'   task and system. The main side-effect is writing
 #'   \file{Task_size_system.pdf} to \code{OUTDIR}.
 #'
 #' @examples
 #' # plotTaskSizes(taskInfo, OUTDIR = "./out")
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
    P = geneCounts %>% mutate(System = gsub("METABOLISM","M.", Depth1)) %>%
            ggplot(., aes(x=System, y=geneCount)) +
            geom_boxplot(fill="green4") +
            labs(# title = "TCGA",
                x="", y="Task size")+
            theme_minimal() +
            theme(legend.position="none",
                axis.text.x = element_text(angle = 90, hjust=1, vjust=0))
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
 #' @return Invisibly returns the ggplot2 object used to render the figure.
 #'   The primary side-effect is writing \file{TaskSummary.pdf} to
 #'   \code{out_dir}.
 #'
 getCirPackingPlot_CellFie = function(cfs, out_dir) {
    edges = cfs%>%mutate(Depth0="Metabolism")%>%.[,c("Depth0","Depth1")]%>%
        unique()%>%setNames(c("from", "to")) %>%
        rbind(cfs[, c("Depth1","Depth2")]%>%unique()%>%setNames(c("from", "to"))) %>%
        rbind(cfs[, c("Depth2","Depth3")]%>%unique()%>%setNames(c("from","to")))
    vertices = cfs%>%mutate(Depth0="Metabolism")%>%.[,c("Depth0","Depth3")]%>%unique()%>%
        pull(Depth0)%>%table()%>%as.data.frame()%>%setNames(c("name", "size")) %>%
        rbind(cfs[,c("Depth1","Depth3")]%>%unique()%>%pull(Depth1)%>%table()%>%
        as.data.frame()%>%setNames(c("name", "size"))) %>%
        rbind(cfs[,c("Depth2","Depth3")]%>%unique()%>%pull(Depth2)%>%table()%>%
        as.data.frame()%>%setNames(c("name", "size"))) %>%
        rbind(cfs%>%pull(Depth3)%>%unique()%>%table()%>%as.data.frame()%>%
        setNames(c("name","size"))) %>%
        mutate(size = size) %>%
        mutate(shortName=gsub("METABOLISM","M.",paste0(name,"(",size,")")))
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
 #' @return A single data.frame (invisibly) containing the concatenated,
 #'   post-processed CellFie results across groups. Side-effects: multiple
 #'   summary files and PDF figures are written to \code{file.path(outdir, "../Summary/")},
 #'   including \file{TaskHierarchy_summary.xlsx}, \file{cfs.RDS}, task score
 #'   tables for each depth, and several PDF plots.
 #'
 #' @examples
 #' # processCellFieOutput(combined = FALSE, outdir = "./CellFieOut/", SampleNames = samples, meta = meta, samples2keep = samples)
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
            mutate(Depth3.ID = TaskID) %>%
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
            # mutate(Group=gsub("_.*_", "_", Sample)) %>% ##########
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
            mutate(TaskScore=MeanTaskScore)
    saveRDS(Depth1, paste0(outdir, "/../Summary/TaskScores_", Depth,".RDS"))
    ## Task depth2
    Depth = "Depth2"
    w = (cfs$Sample %>% unique() %>% length())/6.5 + 5
    h = ((cfs[[Depth]] %>% unique() %>% length())/15) * 3.5 + 2
    Depth2 = sumUpTaskScores(cfs, Depth=Depth, h=h, w=w,
            OUTDIR=paste0(outdir, "/../Summary/"))%>%
            mutate(TaskScore=MeanTaskScore)
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
        print(ggplot2::ggplot(Depth1, aes(x=Cancer_Condition, y=MeanTaskScore, fill=Cancer_type)) +
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
#' @param dat numeric matrix or data.frame of CellFie scores with features as
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

limmaTest_CellFie = function(dat, paired=TRUE, currentCovariate=NULL, 
    checkCovariate=FALSE, meta=meta, paired_variable="Patient") {
    meta = filter(meta, Sample%in%colnames(dat))
    # Make sure the sample are in the same order in both meta and dat:
    dat = dat[, meta$Sample]
    meta = meta %>% mutate(across(where(is.factor), droplevels))
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
    
    fit.con = limma::eBayes(fit, robust = TRUE)
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
#'   `Depth3`) referenced by `depth`.
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

runLimmaCellFie = function(depth, cfs, meta, OUTDIRV, w=9.5, height=23) {
    dir.create(OUTDIRV, recursive = TRUE, showWarnings = FALSE)
    rownames(cfs) = NULL
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
        invisible(capture.output({library(matrixStats)}))
        keep = rowSds(sub%>%as.matrix(), na.rm = T) > 0
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
            mutate(logFC=round(logFC, 2)) %>%
            mutate(Depth2 = stringr::str_sub(Depth2, 1, 3)) %>%
            mutate(Depth1 = stringr::str_sub(Depth1, 1, 3)) %>% 
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

doVolcano_CellFie = function(resi, NAME, 
    adj.P.Val.cutoff=0.05, logFC.cutoff=0, lfcShrink=F) {
  selectedCols = c("Task", "logFC", "P.Value", "adj.P.Val", "Depth")
  toExport = resi[, selectedCols]
  toExport = toExport[order(toExport$P.Value),]

  tmp = as.data.frame(toExport) %>%
    filter(!is.na(adj.P.Val)) %>%
    mutate(sig = ifelse(adj.P.Val<adj.P.Val.cutoff&logFC>logFC.cutoff, "up", 
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
SampleSizeCircularBarplot = function(meta, title="Pan-cancer tasks", OUT_DIR, h=4, w=4.5) {
    temp = meta %>% dplyr::count(Cancer_type, Condition) %>%
        tidyr::spread(Condition, n) %>%
        dplyr::mutate(Sum = ifelse(is.na(Normal), 0, Normal)+ifelse(is.na(Tumor), 0, Tumor)) %>% 
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

ActivityBoxplot_byCancer = function(meta, cf, OUT_DIR, w=2.75, h=3.0, title="Pan-cancer tasks") {
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
#' @return Invisibly returns the ggplot object after writing a PDF to
#'   \code{OUT_DIR}. The primary effect is the written PDF file
#'   \code{meanTaskScores_system.pdf}.
#'
#' @examples
#' \dontrun{
#' ActivityBoxplot_bySystem(cf, OUT_DIR = "./plots/", w = 2, h = 3.5, title = "My Study")
#' }
#' @export

ActivityBoxplot_bySystem = function(cf, OUT_DIR, w=2, h=3.5, title="Pan-cancer tasks") {
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

SigNrBarplot = function(rlst_list, sig.cutoff=0.05, w=2.75, h=3.0, OUT_DIR){
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
        mutate(Cancer_type=gsub("_.*","",Cancer_type)) %>%
        mutate(Regulation=gsub(".*_","",Regulation)) %>%
        mutate(Count = ifelse(Regulation=="up",Count,-Count)) %>% arrange(.,Count)
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


#' Boxplot of limma differential effect (t statistic) by system (Depth1)
#'
#' Build a boxplot summarizing the limma test statistic (t) for tasks
#' grouped at the system level (Depth1). For each cohort in
#' \code{rlst_list} the function expects a data.frame \code{rslt_of_interest}
#' containing at least the limma columns \code{t}, \code{logFC}, and
#' \code{adj.P.Val}. The function merges these results with the
#' task hierarchy in \code{cf} to map tasks (Depth3) to systems (Depth1),
#' orders systems by median effect, and writes a PDF plot file to
#' \code{OUT_DIR}.
#'
#' @param rlst_list named list. Typical input is the output from
#'   \code{runLimmaCellFie} or \code{limmaTest_CellFie}; each element must
#'   contain a data.frame named \code{rslt_of_interest} with at least
#'   the columns \code{t}, \code{adj.P.Val}, and \code{logFC}.
#' @param cf data.frame. CellFie output containing at minimum the columns
#'   \code{Depth1} and \code{Depth3} so tasks can be mapped to systems.
#' @param OUT_DIR character. Directory where the PDF output file will be written.
#' @param w numeric. Figure width in inches (default 2).
#' @param h numeric. Figure height in inches (default 3.5).
#' @param title character. Plot title (default: "Pan-cancer tasks").
#'
#' @return Invisibly returns the ggplot object used to draw the figure.
#'   The primary side effect is the written PDF file
#'   \code{diffEffectBoxplot_system.pdf} inside \code{OUT_DIR}.
#'
#' @examples
#' \dontrun{
#' # rlst_list <- list(OV = list(rslt_of_interest = ov_df), BRCA = list(rslt_of_interest = brca_df))
#' diffEffectBoxplot_bySystem(rlst_list, cf, OUT_DIR = "./plots/", w = 3, h = 4, title = "Study Systems")
#' }
#' @export

diffEffectBoxplot_bySystem = function(rlst_list, cf, OUT_DIR, w=2, h=3.5, title="Pan-cancer tasks") {
    B_Z_N0 = lapply(names(rlst_list), function(cancer){
        temp = rlst_list[[cancer]]$rslt_of_interest
        temp$Cancer = cancer
        temp %>% tibble::rownames_to_column(var="Depth3") 
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>% 
    dplyr::left_join(cf %>% dplyr::select(Depth1, Depth3) %>% unique())

    System_orderZ = B_Z_N0 %>% unique() %>% group_by(Depth1) %>%
        dplyr::summarise(meanEffect = median(t, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(Depth1) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(System=factor(Depth1, levels = System_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=System, y=t)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="System", y="Diff effect")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_system.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}

#' Boxplot of limma differential effect (t statistic) by cancer type
#'
#' Create a boxplot summarizing the limma test statistic (t) for each
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

diffEffectBoxplot_byCancer = function(rlst_list, OUT_DIR, w=2.75, h=3.0, title="Pan-cancer tasks") {
    B_Z_N0 = lapply(names(rlst_list), function(cancer){
        temp = rlst_list[[cancer]]$rslt_of_interest
        temp$Cancer_type = cancer
        temp %>% tibble::rownames_to_column(var="Depth3") 
    }) %>% Reduce(rbind, .) %>% as.data.frame() %>% unique()

    Cancer_orderZ = B_Z_N0 %>% unique() %>% group_by(Cancer_type) %>%
        dplyr::summarise(meanEffect = median(t, na.rm=T), .groups = "drop") %>% 
        arrange(meanEffect) %>% pull(Cancer_type) %>% unique()

    B_Z_N = B_Z_N0 %>%
        dplyr::mutate(Cancer=factor(Cancer_type, levels = Cancer_orderZ)) %>%
        ggplot2::ggplot(., ggplot2::aes(x=Cancer, y=t)) +
        ggplot2::geom_boxplot(fill="green4") +
        ggplot2::labs(title = title,
            x="Cancer", y="Diff effect")+
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust=1, vjust=0))

    pdf(paste0(OUT_DIR, "/diffEffectBoxplot_cancer.pdf"), w=w, h=h)
        print(B_Z_N)
    dev.off()
}


