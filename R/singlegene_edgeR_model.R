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
#' @param covariates Matrix or data.frame containing the covariates of the
#' models. Rows represent samples, while columns represent the different
#' covariates. This argument is ignored if **design_mat_singlegene**
#' is provided.
#' @param interactions A list of characters containing the covariates to use
#' for each single gene model. Each element of the list should be a character
#' vector corresponding to colnames in **covariates**, for each genes the
#' specified columns will be used to create the design matrix.
#' This argument is ignored if **design_mat_singlegene** is provided.
#' @param design_mat A design matrix (if the model to run is one or if the
#' matrix doesn't change across models) or a list of design matrices for each
#' edgeR single gene model.
#' If provided as list, the number of design matrices should be equal to the
#' number of rows of **response_var**.
#' @param offset_singlegene A vector (if the model to run is one or if the
#' offsets don't change across models) or a list of vectors containing the
#' offsets for each single gene edgeR model.
#' If provided as list, the number of vectors should be equal to the number
#' of rows of **response_var**.
#' @param fit_all An all gene edgeR model needed to extract the common and
#' tagwise dispersion
#' @param threads Number of threads to use for parallelization (Default is 1)
#' @import parallel edgeR stringr

singlegene_edgeR_model <- function( response_var,
                                    interactions = NULL,
                                    covariates = NULL,
                                    offset_singlegene = NULL,
                                    design_mat,
                                    fit_all,
                                    threads = 1) {

    response_var <- t(response_var+1)

    if(is.null(design_mat) & (is.null(interactions)+is.null(covariates))>0) {
        stop(str_wrap("You should provide design_mat if interactions and
                        covariates are not provided"))
    }
    if(!is.null(covariates)){
        if(ncol(response_var)!=nrow(covariates)) stop(str_wrap("response_var
            and covariates should have the same number of samples"))
        if(!identical(colnames(response_var), rownames(covariates))) warning(
            str_wrap("response_var and covariates have different sample names,
                    assuming that samples in response_var and covariates are
                    in the same order"))
    }


    fit_list <- mclapply(1:nrow(response_var), function(x) {

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
            y_gene$samples$norm.factors <- fit_all$samples$norm.factors
            y_gene$offset <- fit_all$offset[1,]
        }
        y_gene$common.dispersion <- fit_all$common.dispersion
        y_gene$tagwise.dispersion <- fit_all$tagwise.dispersion[x]
        fit <- glmFit(y_gene, design_mat[[x]])
        return(fit)
        }, mc.cores = threads)
    names(fit_list) <- rownames(response_var)
    return(fit_list)
}
