#' Integration of omics data running single gene negative binomial edgeR models
#' @description This function is useful for users which need to run single gene
#' integration models, in particular if for each gene the covariates change.
#' For example if you need to integrate gene expression with copy number
#' variations, for each gene you need to integrate its expression levels with
#' CNV status, updating the covariate at each step. In this case you can set
#' **interactions** to auto (default), providing the expression of all genes in
#' **response_var** and the relative CNV values (numeric) in **covariates**.
#' The function will automatically use CNV values of each gene as covariates
#' for the single gene models.
#' The user can specify which covariates should be used for each gene in the
#' **interaction** argument.
#' The user can manually provide a list of the design matrices for single gene
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
#' This argument is ignored if **design_mat_singlegene** is provided.
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
#' @param steady_covariates Character vector containing column names of
#' **covariates** corresponding to covariates that should be included in all
#' single gene models (default=NULL).
#'
#' @return A list containing the results of all the edger single gene models,
#' the pvalues and the coefficients for each gene
#' @import parallel edgeR BiocParallel stringr
#' @export

run_edgeR_integration <-  function( response_var,
                                    covariates,
                                    interactions = "auto",
                                    design_mat_allgene = NULL,
                                    offset_allgene = NULL,
                                    offset_singlegene = NULL,
                                    normalize=T,
                                    norm_method = "TMM",
                                    steady_covariates = NULL,
                                    reference = NULL,
                                    BPPARAM = BiocParallel::SerialParam()) {

    tmp <- data_check(response_var = response_var,
                      covariates = covariates,
                      interactions = interactions)
    response_var <- tmp$response_var
    covariates <- tmp$covariates
    interactions <- tmp$interactions
    tmp <- unlist(lapply(interactions, length))
    single_cov=F
    if(sum(tmp==1)==length(tmp)){
      single_cov=T
    }

    if(normalize==F) offset_allgene <- matrix(1, ncol(response_var),
                                              nrow(response_var))
    fit_all <- allgene_edgeR_model(
        response_var = response_var,
        design_mat_allgene = design_mat_allgene,
        offset_allgene = offset_allgene,
        norm_method = norm_method
    )

    tmp <- covariates_check(response_var=response_var,
                            covariates=covariates,
                            interactions = interactions,
                            steady_covariates=steady_covariates,
                            reference=reference)
    covariates <- tmp$covariates
    response_var <- tmp$response_var
    interactions <- tmp$interactions
    original_id <- tmp$original_id
    fformula <- generate_formula(interactions = interactions)

    fit_gene <- singlegene_edgeR_model(
        response_var = response_var,
        covariates=covariates,
        fit_all = fit_all,
        fformula = fformula,
        offset_singlegene = offset_singlegene,
        BPPARAM = BPPARAM
    )

    model_res <- edger_coef_test(fit_gene,
                                 BPPARAM = BPPARAM)
    coef_pval_mat <- building_result_matrices(model_results = model_res,
                                              type = "edgeR",
                                              single_cov = single_cov)
    tmp <- lapply(fit_gene, function(x) as.data.frame(residuals(x)))
    rresiduals <- rbind.fill(tmp)
    colnames(rresiduals) <- rownames(response_var)
    rownames(rresiduals) <- names(tmp)
    results <- list(
      model_results = model_res,
      coef_data = coef_pval_mat$coef,
      pval_data = coef_pval_mat$pval,
      residuals = rresiduals)
    if(nrow(original_id)>0){
      results <- id_conversion(dictionary = original_id,
                               results = results)
    }
    return(results)
}


################################################################
#' @export
allgene_edgeR_model <- function(response_var,
                                design_mat_allgene = NULL,
                                offset_allgene = NULL,
                                norm_method = "TMM") {

  response_var <- t(response_var+1)
  if(!norm_method %in% c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
    stop(str_wrap("norm_method should be one of TMM, TMMwsp, RLE,
                        upperquartile, none"))
  }

  if(is.null(design_mat_allgene)) {
    design_mat_allgene <- model.matrix(~1,
                                       data = as.data.frame(t(response_var)))
  }
  if(ncol(response_var) != nrow(design_mat_allgene)) {
    stop(str_wrap("Number of samples differ between respone_var and
                        design_mat_allgene"))
  }
  y_all <- DGEList(counts = response_var)
  y_all <- calcNormFactors(y_all, method = norm_method)
  if (!is.null(offset_allgene)) y_all$offset <- offset_allgene
  y_all <- estimateGLMCommonDisp(y_all, design_mat_allgene)
  y_all <- estimateGLMTagwiseDisp(y_all, design_mat_allgene)
  fit_all <- glmFit(y_all, design_mat_allgene)
  fit_all$common.dispersion <- y_all$common.dispersion
  fit_all$tagwise.dispersion <- y_all$tagwise.dispersion
  return(fit_all)
}


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
#' @import parallel edgeR

singlegene_edgeR_model <- function( response_var,
                                    covariates,
                                    offset_singlegene = NULL,
                                    fformula,
                                    fit_all,
                                    BPPARAM) {

  response_var <- t(response_var+1)
  fformula <- lapply(seq_along(fformula), function(x){
    list(formula=fformula[[x]], var_name=names(fformula)[x])
  })
  names(fformula) <- sapply(fformula, function(x) x$var_name)
  if(!is.null(offset_singlegene)) {
    if(is.list(offset_singlegene)) {
      if(length(offset_singlegene)!=nrow(response_var)) stop(str_wrap(
        "If provided as list, the length of offset_singlegene
            should be equal to the number of genes in response_var"))
    }
  }
  fit_list <- bplapply(fformula, def_edger,
                       response_var=response_var,
                       covariates=covariates,
                       fit_all=fit_all,
                       offset_singlegene=offset_singlegene,
                       BPPARAM = BPPARAM)
  names(fit_list) <- rownames(response_var)
  return(fit_list)
}


##########################################################
#' Define edgeR model

def_edger <- function(formula,
                      response_var,
                      covariates,
                      offset_singlegene = NULL,
                      fit_all){
  design <- model.matrix(formula$formula, data=covariates)
  dge <- DGEList(counts = t(response_var[formula$var_name,]))
  if(!is.null(offset_singlegene)) {
    if(is.list(offset_singlegene)) {
      dge$offset <- offset_singlegene[[formula$var_name]]
    } else {
      dge$offset <- offset_singlegene
    }
  } else {
    dge$offset <- fit_all$offset[1,]
  }
  dge$samples$norm.factors <- fit_all$samples$norm.factors
  dge$common.dispersion <- fit_all$common.dispersion
  tmp <- which(rownames(fit_all$coefficients)==formula$var_name)
  dge$tagwise.dispersion <- fit_all$tagwise.dispersion[tmp]
  fit <- glmFit(dge, design)
  return(fit)
}


#' Perform likelihood ratio test for all the coefficients
#' contained in the edgeR fitted model
#'
#' @param fit A fitted edgeR model
#' @param threads Number of threads to use
#'
#' @return A data frame containing the results of the likelihood ratio test
#' for all the coefficients present in the model
#' @import parallel edgeR
#'
#' @examples
#' edger_coef_test(fitted_model, threads = 2)
#'
edger_coef_test <- function(fit_list,
                            BPPARAM) {
  top_list <- bplapply(fit_list, def_coef_test,
                       BPPARAM = BPPARAM)
  return(top_list)
}

#########################################################
#' Define coef test function

def_coef_test <- function(fit){

  coef <- colnames(fit$coefficients)
  lrt <- lapply(coef, function(z) {
    glmLRT(fit, coef = z)
  })
  names(lrt) <- coef
  top <- as.matrix(sapply(names(lrt),
                          function(z) topTags(lrt[[z]])$table))
  coef <- sapply(names(lrt), function(z) {
    lrt[[z]]$coefficients[, z]
  })
  top <- rbind(top, coef)
  top <- as.data.frame(top)
  return(top)
}


