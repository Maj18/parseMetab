\name{visualizeGSVSscoresGroup}
\alias{visualizeGSVSscoresGroup}
\title{Visualize GSVA pathway means across groups and matched normals}
\description{
Create side-by-side boxplots that summarize mean GSVA pathway activities
across provided groups and, separately, for matched normal samples.
This supports cohort-level comparisons (for example across cancer types).
}
\usage{
visualizeGSVSscoresGroup(GSVA_limma_rslt_gsva, paired_data, OUTDIR)
}
\arguments{
  \item{GSVA_limma_rslt_gsva}{A named list of GSVA result matrices or data frames.
    Each element should correspond to a group (list names indicate groups,
    e.g. cancer types) and contain pathway identifiers as row names and
    sample submitter IDs as column names.}
  \item{paired_data}{A named list of metadata data frames (one per group) that
    corresponds to `GSVA_limma_rslt_gsva`. Each metadata table must contain at
    least the columns \code{sample_type} (e.g. "Solid Tissue Normal" or
    "Primary tumor") and \code{sample.submitter_id} used to select normal samples.}
  \item{OUTDIR}{Character scalar. Directory where output plots (PDF) will be written.
    If \code{NULL}, plots may be returned instead of being saved.}
}
\value{
A list with the ggplot2 objects for (1) all samples and (2) normals-only, and
the path to the saved PDF file (written as \file{GSVA_samples.pdf} in
\code{OUTDIR}) when \code{OUTDIR} is provided. The primary side-effect is the
creation of the PDF file in \code{OUTDIR}.
}
\author{Yuan Li \email{yuan.li@nbis.se}}
\examples{
}
