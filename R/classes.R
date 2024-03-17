#' MultiOmics Class
#' @description
#' S4 class containing the output of a multiomics integration. It's a list in
#' which each element represents the result of an integration. If all the
#' available omics are provided, it will be a list of integrations:
#' **gene_genomic_res**, **mirna_cnv_res**, **tf_res**, **tf_mirna_res** and
#' **mirna_target_res**
#' @return MultiOmics Class
setClass("MultiOmics",
    representation("list")
)

#' MultiClass Class
#' @description
#' S4 class containing the output of a single integration integration, for
#' which classes has been provided. It's a list in which each element
#' represents the result of the integration for a given class. The length will
#' be equal to the number of classes defined.
#'@return MultiOmics Class
setClass("MultiClass",
    representation("list")
)
