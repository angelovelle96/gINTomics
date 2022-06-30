#' Single gene omics integration
#'
#' @description Compute single gene edgeR's models for each of the genes
#' provided in **response_var**
#' @param response_var Matrix or data.frame containing the expression values
#' for each model. Rows represent samples, while each column represents
#' the different response variables of the models. The number of rows should
#' be identical to that of **covariates** or of design matrices provided in
#' **design_mat_singlegene** and the order of the samples should be the same.
#' In case you need to run a single model you can provide the response
#' variable as an atomic vector, anyway the response variables should always
#' be numeric values.
#' @param design_mat A list of design matrices for
#' each edgeR single gene model, the number of design matrices should be equal
#' to the number of rows of **response_var**.
#' @param offset_singlegene A list of vector containing the offsets for each
#' single gene edgeR model.
#' The number of vectors should be equal to the number of rows
#' of **response_var**. As default, offsets are taken from the all gene fitting.
#' @param fit_all An all gene edgeR fitted model needed to extract the common
#' and tagwise dispersion and the offsets
#' @param threads Number of threads to use for parallelization (Default is 1)
#' @import parallel edgeR stringr

singlegene_edgeR_model <- function( response_var,
                                    offset_singlegene = NULL,
                                    design_mat,
                                    fit_all,
                                    threads = 1) {

    response_var <- t(response_var+1)
    if(!identical(rownames(response_var), names(design_mat))) warning(
        str_wrap("response_var and design_mat have different gene names,
                    assuming that genes in response_var and design_mat are
                    in the same order"))
    fit_list <- mclapply(1:nrow(response_var), function(x) {
        if(nrow(design_mat[[x]])!=ncol(response_var)) stop(str_wrap("Number
                of samples in response_var and design_mat_singlegene
                should not differ"))
        y_gene <- DGEList(counts = t(response_var[x,]))
        if(!is.null(offset_singlegene)) {
            if(is.list(offset_singlegene)) {
                if(length(offset_singlegene)!=nrow(response_var)) stop(str_wrap(
                "If provided as list, the length of offset_singlegene
                    should be equal to the number of genes in response_var"))
                y_gene$offset <- offset_singlegene[[x]]
            } else {
                y_gene$offset <- offset_singlegene
            }
        } else {
            y_gene$offset <- fit_all$offset[1,]
        }
        y_gene$samples$norm.factors <- fit_all$samples$norm.factors
        y_gene$common.dispersion <- fit_all$common.dispersion
        y_gene$tagwise.dispersion <- fit_all$tagwise.dispersion[x]
        fit <- glmFit(y_gene, design_mat[[x]])
        return(fit)
        }, mc.cores = threads)
    names(fit_list) <- rownames(response_var)
    return(fit_list)
}
