#' Integration of omics data running single gene negative binomial edgeR models
#' @description This function is useful for users which need to run single gene
#' integration models, in particular if for each gene the covariates change.
#' For example if you need to integrate gene expression with copy number
#' variations, for each gene you need to integrate its expression levels with
#' CNV status, updating the covariate at each step. In this case you can set
#' **cnv_mode** to TRUE, providing the expression of all genes in
#' **response_var** and the relative CNV values (numeric) in **covariates**.
#' The function will automatically use CNV values of each gene as covariates
#' for the single gene models.
#' If **cnv_mode** is FALSE (default) the user should specify which covariates
#' should be used for each gene in the **interaction** argument.
#' The use can manually provide a list of the design matrices for single gene
#' models in **design_mat_singlegene**, in this case the n gene contained in
#' **response_var** is fitted in the  model  using as design matrix the n
#' element of the**design_mat_singlegene** list.
#' Since the function runs a model for each of the genes contained in
#' **response_var**, in order to compute the common and tagwise dispersion it
#' starts running an all gene model without covariates, then the dispersion
#' values and the offsets are passed to the single gene model.
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
#' specified columns will be used to create the design matrix. The list should
#' contain the covariates for each of the genes of **response_var** and they
#' should be in the same order of **response_var** rownames.
#' This argument is ignored if **design_mat_singlegene** is provided or if
#' **cnv_mode** is set to TRUE.
#' @param design_mat_singlegene A list of design matrices for
#' each edgeR single gene model, the number of design matrices should be equal
#' to the number of rows of **response_var**.
#' @param design_mat_allgene A design matrix for the all gene edgeR model
#' @param offset_allgene A matrix containing the offsets of the all gene edgeR
#' model
#' @param offset_singlegene A list of vector containing the offsets for each
#' single gene edgeR model.
#' The number of vectors should be equal to the number of rows
#' of **response_var**. As default, offsets are taken from the all gene fitting.
#' @param norm_method Normalization method to use for the all gene edgeR
#' model (Default is 'TMM'). The normalization factors will be passed to single
#' gene models for normalization.
#' @param threads Number of threads to use for parallelization (Default is 1)
#' @param cnv_mode logical indicating if the model should perform copy number
#' integration (default=FALSE). If set to TRUE the function will generate the
#' design matrix of a given gene by searching that gene among the columns of
#' **covariates**, so each gene in **response_var** should have a
#' corresponding column in **covariates** with the same name.
#' @param steady_covariates Character vector containing column names of
#' **covariates** corresponding to covariates that should be included in all
#' single gene models (default=NULL).
#'
#' @return A list containing the results of all the edger single gene models,
#' the pvalues and the coefficients for each gene
#' @import parallel edgeR stringr
#' @importMethodsFrom RUVSeq

run_edgeR_integration <-  function( response_var,
                                    covariates = NULL,
                                    interactions = NULL,
                                    design_mat_allgene = NULL,
                                    design_mat_singlegene = NULL,
                                    offset_allgene = NULL,
                                    offset_singlegene = NULL,
                                    norm_method = "TMM",
                                    steady_covariates = NULL,
                                    cnv_mode = F,
                                    threads = 1,
                                    reference=NULL) {


    if(is.null(design_mat_singlegene)){
      tmp <- data_check(response_var = response_var,
                        covariates = covariates,
                        cnv_mode = cnv_mode,
                        interactions = interactions,
                        steady_covariates=steady_covariates)
      response_var <- tmp$response_var
      covariates <- tmp$covariates
      interactions <- tmp$interactions
    }

    fit_all <- allgene_edgeR_model(
        response_var = response_var,
        design_mat_allgene = design_mat_allgene,
        offset_allgene = offset_allgene,
        norm_method = norm_method
    )
    if(is.null(design_mat_singlegene)){
        design_mat_singlegene <- generate_design(
            response_var = response_var,
            covariates = covariates,
            interactions = interactions,
            steady_covariates = steady_covariates,
            cnv_mode = cnv_mode,
            reference = reference,
            threads=threads)
    }
    fit_gene <- singlegene_edgeR_model(
        response_var = response_var,
        fit_all = fit_all,
        design_mat = design_mat_singlegene,
        offset_singlegene = offset_singlegene,
        threads = threads
    )
    model_res <- edger_coef_test(fit_gene,
        threads = threads
    )

    coef_pval_mat <- building_result_matrices(model_results = model_res,
                                              type = "edgeR")

    tmp <- mclapply(fit_gene, function(x) as.data.frame(residuals(x)),
                    mc.cores = threads)
    rresiduals <- rbind.fill(tmp)
    colnames(rresiduals) <- rownames(response_var)
    rownames(rresiduals) <- names(tmp)
    results <- list(
      model_results = model_res,
      coef_data = coef_pval_mat$coef,
      pval_data = coef_pval_mat$pval,
      residuals = rresiduals)
    return(results)
}
