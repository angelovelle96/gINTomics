#' Running linear model integration
#' @importFrom BiocParallel bplapply
#' @importFrom stats residuals
.run_lm_integration <- function(response_var,
                               covariates,
                               interactions="auto",
                               steady_covariates=NULL,
                               reference=NULL,
                               normalize=TRUE,
                               norm_method="TMM",
                               BPPARAM=SerialParam()){

    tmp <- unlist(lapply(interactions, length))
    single_cov=FALSE
    if(sum(tmp==1)==length(tmp)){
      single_cov=TRUE
    }
  tmp <- .data_check(response_var = response_var,
                      interactions = interactions,
                      covariates = covariates)
    if(is.null(tmp)) return(NULL)
    response_var <- tmp$response_var
    interactions <- tmp$interactions
    covariates <- tmp$covariates

    if(normalize==TRUE) response_var <- .data_norm(response_var,
                                               method = norm_method,
                                               RNAseq = FALSE)
    tmp <- .covariates_check(response_var = response_var,
                            covariates = covariates,
                            interactions = interactions,
                            steady_covariates = steady_covariates,
                            linear = TRUE)
    if(is.null(tmp)) return(NULL)
    covariates <- tmp$covariates
    response_var <- tmp$response_var
    interactions <- tmp$interactions
    steady_covariates <- tmp$steady_covariates
    original_id <- tmp$original_id
    fformula <- .generate_formula(interactions = interactions,
                                 linear=TRUE)
    data <- cbind(response_var, covariates)
    lm_results <-  bplapply(fformula, .def_lm,
                            data=data,
                            BPPARAM = BPPARAM)
    names(lm_results) <- names(interactions)

    coef_pval_mat <- .building_result_matrices(model_results = lm_results,
                                              type = "lm",
                                              single_cov = single_cov)
    tmp <- bplapply(lm_results, function(x) as.data.frame(t(residuals(x))),
                    BPPARAM = BPPARAM)
    rresiduals <- rbind.fill(tmp)
    colnames(rresiduals) <- rownames(response_var)
    rownames(rresiduals) <- names(tmp)

    results <- list(
      coef_data = coef_pval_mat$coef,
      pval_data = coef_pval_mat$pval,
      fdr_data = fdr(coef_pval_mat$pval),
      residuals = rresiduals,
      data = list(response_var = response_var,
                  covariates = covariates,
                  formula = fformula))
    if(nrow(original_id)>0){
      results <- .id_conversion(dictionary = original_id,
                               results = results)
    }
    return(results)
}


##################################
#' Def linear model
#' @importFrom stats lm

.def_lm <- function(formula,
                    data){

  if(length(attr(terms(formula), "term.labels"))>
     as.integer(nrow(data)*0.4)){

    tmp <- as.character(formula[2])
    tmp2 <- .rf_selection(formula = formula,
                          data = data)
    formula <- as.formula(paste(tmp, "~", as.character(tmp2[2])))
  }

  lm_results <- summary(lm(formula, data = data))
  return(lm_results)
}

