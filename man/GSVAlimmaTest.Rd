\name{GSVAlimmaTest}
\alias{GSVAlimmaTest}
\title{
Differential Analysis of KEGG Metabolic Pathway Activity Using GSVA and limma
}
\description{
This function for now only supports raw RNAseq data, it normalizes (logCPM) the raw count data using edgeR and limma, calcualtes GSVA enrichment scores based on the normalized counts for each KEGG metabolic pathway, runs limma trend on the GSVA metabolic activity scores.
}
\usage{
GSVAlimmaTest(dat, paired = TRUE, Feature.filtering = TRUE, currentCovariate = NULL, meta = meta, paired_variable = "Patient", OUTDIR, kegg_metab_db)
}
\arguments{
  \item{dat}{A numeric matrix or data frame of raw transcriptomics counts with
    features (genes) as row names and samples as column names.}
  \item{paired}{Logical. Set to \code{TRUE} to perform a paired analysis when
    paired samples exist (at least two pairs). Paired samples must share the
    same identifier provided in \code{paired_variable}.}
  \item{Feature.filtering}{Logical. If \code{TRUE}, features are filtered to
    retain those with counts greater than 5 in at least 15% of samples in each
    group.}
  \item{currentCovariate}{Optional character. Name of an additional covariate
    (column in \code{meta}) to include in the model alongside \code{Group}.}
  \item{meta}{A data frame of sample metadata that contains at minimum the
    columns \code{Sample} (sample IDs matching \code{colnames(dat)}) and
    \code{Group} (a factor indicating group membership).}
  \item{paired_variable}{Character. Name of the metadata column used to
    identify paired samples (for example, patient or subject ID).}
  \item{OUTDIR}{Character. Directory where output files (plots and tables) will
    be written, e.g. normalized counts, GSVA scores, and the MDS plot.}
  \item{kegg_metab_db}{A list mapping KEGG pathway names or IDs to vectors of
    member gene symbols (one of the outputs produced by \code{getdb_metabolism}).}
}
\details{
This function for now only supports raw RNAseq data, it normalizes (logCPM) the raw count data using edgeR and limma, calcualtes GSVA enrichment scores based on the normalized counts for each KEGG metabolic pathway, runs limma trend on the GSVA metabolic activity scores.
}
\value{
 \item{rslt_of_interest}{Differential result table}
 \item{gsva}{sample-wise GSVA score for each KEGG metabolic pathway}
 \item{fitted.model}{fitted limma model}
}
\references{
}
\author{
Yuan Li {yuan.li@nbis.se}
}
\note{
}


\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (dat, paired = TRUE, Feature.filtering = TRUE, currentCovariate = NULL, 
    meta = meta, paired_variable = "Patient", OUTDIR, kegg_metab_db) 
{
    meta = filter(meta, Sample \%in\% colnames(dat))
    dat = dat[, meta$Sample]
    meta = meta \%>\% mutate(across(where(is.factor), droplevels))
    mapping = setNames(meta$Group, meta$Sample)
    print("Only keep features that have at least min.nr number of values that are >5 in each group...")
    min.nr = round(max((ncol(dat)/2) * 0.15, 2), 0)
    if (Feature.filtering) {
        B = mapping[colnames(dat)] \%>\% as.character()
        featureskeep = apply(dat > 5, 1, function(row) {
            A = row
            sum(aggregate(A ~ B, data = data.frame(A = A, B = B), 
                sum, na.rm = T)[, "A"] >= min.nr) > 1
        })
        dat = dat[featureskeep, ]
    }
    Group = factor(meta$Group)
    dge = edgeR::DGEList(counts = dat)
    dge = edgeR::calcNormFactors(dge)
    pdf(paste0(OUTDIR, "/MDS_plot.pdf"))
    print(limma::plotMDS(dge, col = as.numeric(Group)))
    dev.off()
    baseFormula = ~Group
    if (is.null(currentCovariate)) {
        currentFormula = baseFormula
    }
    else {
        currentFormula = as.formula(paste0("~ ", currentCovariate, 
            " + Group"))
    }
    design = model.matrix(currentFormula, data = meta)
    colnames(design) = gsub("Group", "", colnames(design))
    dat2 = limma::voom(dge, design, plot = T)
    normalizedCount = dat2$E
    write.table(normalizedCount \%>\% as.data.frame() \%>\% tibble::rownames_to_column(var = "Features"), 
        paste0(OUTDIR, "/NormalizedCount.csv"), quote = F, row.names = F, 
        col.names = T, sep = "\t")
    normalizedCount = as.matrix(normalizedCount)
    normalizedCount = matrix(as.numeric(normalizedCount), nrow = nrow(normalizedCount), 
        ncol = ncol(normalizedCount), dimnames = dimnames(normalizedCount))
    param = GSVA::gsvaParam(normalizedCount, kegg_metab_db, minSize = 5, 
        maxSize = 500)
    gsva_scores = GSVA::gsva(param, verbose = TRUE)
    write.table(gsva_scores \%>\% as.data.frame() \%>\% tibble::rownames_to_column(var = "Features"), 
        paste0(OUTDIR, "/GSVAscores.csv"), quote = F, row.names = F, 
        col.names = T, sep = "\t")
    if (paired) {
        dupcor = limma::duplicateCorrelation(gsva_scores, design, 
            block = meta[, paired_variable, drop = T])
        fit = limma::lmFit(gsva_scores, design, block = meta[, 
            paired_variable, drop = T], correlation = dupcor$consensus)
    }
    else {
        fit = limma::lmFit(gsva_scores, design)
    }
    gssizes = sapply(kegg_metab_db[rownames(fit$coefficients)], 
        function(gs) {
            length(gs)
        })
    fit.con = limma::eBayes(fit, robust = TRUE, trend = gssizes)
    rlst_interest = limma::topTable(fit.con, n = Inf, coef = levels(meta$Group)[2]) \%>\% 
        arrange(-t)
    return(list(rslt_of_interest = rlst_interest, gsva = gsva_scores, 
        fitted.model = fit))
  }
}

