#'  Setting method for extracting results
#' @importFrom reshape2 melt
#' @importFrom dplyr left_join
#' @importFrom stats IQR quantile
#' @importFrom plyr rbind.fill

#' @param model_results The model results object from which to extract results.
#' @param ... Additional arguments to be passed to specific methods.
#'
#' @rdname extract_model_res
#'
#' @export
#'
setGeneric("extract_model_res", function(model_results, ...)
  standardGeneric("extract_model_res"))

