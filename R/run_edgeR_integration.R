#' Integration of omics data running single gene edgeR models
#' @description This function is useful for users which need to run single gene
#' integration models, in particular if for each gene the covariates change.
#' For example if you need to integrate gene expression with copy number
#' variations, for each gene you need to integrate its expression levels with
#' CNV status, updating the covariate at each step. In this case you can
#' provide the expression of all genes in **response_var** and the relative
#' CNV values (numeric) as a list of design matrices in
#' **design_mat_singlegene**.
#' The n gene contained in **response_var** is fitted in the negative binomial
#' model of the edgeR package using as design matrix the n element of the
#' **design_mat_singlegene** list. Alternatively, if the covariates are the
#' same for each gene, the user can provide them as a dataframe in
#' **covariates**, in this case all the columns will be included in the design
#' matrix. Since the function runs a model for each of the genes contained in
#' **response_var**, in order to compute the common and tagwise dispersion it
#' starts running an all gene model without covariates, then the dispersion
#' values are passed to the single gene model.
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
#' @param design_mat_singlegene A design matrix (if the model to run is one or if the
#' matrix doesn't change across models) or a list of design matrices for each
#' edgeR single gene model.
#' If provided as list, the number of design matrices should be equal to the
#' number of rows of **response_var**.
#' @param design_mat_allgene A design matrix for the all gene edgeR model
#' @param offset_allgene A matrix containing the offsets of the all gene edgeR
#' model
#' @param offset_singlegene A vector (if the model to run is one) or
#' a list of vector containing the offsets for each single gene edgeR model.
#' The number of vectors should be equal to the number of rows
#' of **response_var**.
#' @param norm_method Normalization method to use for the all gene edgeR
#' model (Default is 'TMM'). The normalization factors will be passed to single
#' gene models for normalization. Ignored if the offsets are provided.
#' @param threads Number of threads to use for parallelization (Default is 1)
#'
#' @return A list containing the results of all the edger single gene models,
#' the pvalues and the coefficients for each gene
#' @import parallel edgeR stringr


run_edgeR_integration <-  function( response_var,
                                    covariates = NULL,
                                    interactions = NULL,
                                    design_mat_allgene = NULL,
                                    design_mat_singlegene = NULL,
                                    offset_allgene = NULL,
                                    offset_singlegene = NULL,
                                    norm_method = "TMM",
                                    threads = 1) {

    if (is.atomic(response_var) & is.vector(response_var)) {
        message("response_var is an atomic vector, converting to matrix")
        tmp <- as.matrix(response_var)
        rownames(tmp) <- names(response_var)
        response_var <- tmp
    }

    response_var <- as.matrix(response_var)

    if (!is.matrix(response_var)) {
        stop(str_wrap("response_var should be a data.frame,
            a matrix or an atomic vector"))
    }

    if(!is.null(covariates)){
        if (is.atomic(covariates)& is.vector(covariates)) {
            message("covariates is an atomic vector, converting to data.frame")
            tmp <- as.data.frame(covariates)
            rownames(tmp) <- names(covariates)
            covariates <- tmp
        }
        if (!is.data.frame(covariates)) {
            stop("covariates should be a data.frame or an atomic vector")
        }
    }

    y_all <- allgene_edgeR_model(
        response_var = response_var,
        covariates = covariates,
        design_mat_allgene = design_mat_allgene,
        offset_allgene = offset_allgene,
        norm_method = norm_method
    )
    y_gene <- singlegene_edgeR_model(
        response_var = response_var,
        covariates = covariates,
        interactions=interactions,
        y_all = y_all,
        design_mat = design_mat_singlegene,
        offset_singlegene = offset_singlegene,
        threads = threads
    )
    model_res <- edger_coef_test(y_gene,
        threads = threads
    )
    coef_pval_mat <- building_edger_result_matrices(model_results = model_res)
    results <- list(
        model_results = model_res,
        coef_data = coef_pval_mat$coef,
        pval_data = coef_pval_mat$pval
    )

    return(results)
}
