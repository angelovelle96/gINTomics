#' Example data for a standard workflow.
#' This is an example dataset containing a MultiAssayExperiment of 20 ovarian
#' cancer (OVC) patients extracted from the Cancer Genome Atlas (TCGA) database.
#' The object contains all the available input data types: Gene expression data,
#' miRNA expression data, gene methylation data, gene Copy Number Variations
#' and miRNA Copy Number Variations.
#' @docType data
#' @name mmultiassay_ov
#' @return An object of class \linkS4class{MultiAssayExperiment}.
#' @format MultiAssayExperiment
#' @usage data(mmultiassay_ov)
#' "mmultiomics_ov"
#' @examples
#' # example code
#' data(mmultiassay_ov)
#' mmultiassay_ov
NULL

#' miRNA IDs.
#' Dataset containing lastly definition of miRNAs (Names, Accessions,
#' Sequences, Families and others) from different miRBase versions
#' (From miRBase version 6 to version 22).
#' @docType data
#' @name mirna_hsa
#' @return An object of class \linkS4class{data.frame}.
#' @format data.frame
#' @usage data(mirna_hsa)
#' "mirna_hsa"
#' @examples
#' # example code
#' data(mirna_hsa)
#' head(mirna_hsa)
NULL
