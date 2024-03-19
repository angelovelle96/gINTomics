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
#' data("ov_test_tcga_omics")
#' tmp <- lapply(mmultiassay_ov@ExperimentList, function(x) x[1:200,])
#' mmultiassay_ov <- MultiAssayExperiment(experiments = tmp)
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' head(data_table)
#' @rdname extract_model_res
#' @export
setGeneric("extract_model_res", function(model_results, ...)
    standardGeneric("extract_model_res"))



