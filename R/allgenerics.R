#' Setting method for extracting results
#' @importFrom reshape2 melt
#' @importFrom dplyr left_join
#' @importFrom stats IQR quantile
#' @importFrom plyr rbind.fill
#' @param model_results The model results object from which to extract results.
#' @param ... Additional arguments to be passed to specific methods.
#' @return A dataframe containing the results of all the integration models
#' provided
#' @examples
#' # example code
#' library(MultiAssayExperiment)
#' data("mmultiassay_ov")
#' tmp <- lapply(experiments(mmultiassay_ov), function(x) x[1:20,])
#' mmultiassay_ov <- MultiAssayExperiment(experiments = tmp)
#' gene_cnv_matrix <- t(as.matrix(assay(mmultiassay_ov[["cnv_data"]])))
#' gene_exp_matrix <- t(as.matrix(assay(mmultiassay_ov[["gene_exp"]])))
#' cnv_integration <- run_cnv_integration(
#'     expression = gene_exp_matrix,
#'     cnv_data = gene_cnv_matrix
#' )
#' data_table <- extract_model_res(cnv_integration)
#' head(data_table)
#' @rdname extract_model_res
#' @export
setGeneric("extract_model_res", function(model_results, ...)
    standardGeneric("extract_model_res"))



