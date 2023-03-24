run_lm_integration <- function(response_var,
                               covariates,
                               interactions="auto",
                               threads=1,
                               step=F,
                               steady_covariates=NULL,
                               reference=NULL,
                               BPPARAM=NULL){

  if(is.null(BPPARAM)) BPPARAM <- MulticoreParam(workers = threads,
                                                 exportglobals = F,
                                                 force.GC = T,
                                                 exportvariables=F)
  tmp <- data_check(response_var = response_var,
                      interactions = interactions,
                      covariates = covariates,
                      linear=T,
                      steady_covariates = steady_covariates)
    response_var <- tmp$response_var
    interactions <- tmp$interactions
    covariates <- tmp$covariates

    lm_results <- lm_singlegene(response_var = response_var,
                                covariates = covariates,
                                interactions = interactions,
                                threads = threads,
                                step = step,
                                steady_covariates = steady_covariates,
                                reference = reference,
                                BPPARAM = BPPARAM)

    coef_pval_mat <- building_result_matrices(model_results = lm_results,
                                              type = "lm")
    tmp <- lapply(lm_results, function(x) as.data.frame(t(residuals(x))))
    rresiduals <- rbind.fill(tmp)
    colnames(rresiduals) <- rownames(response_var)
    rownames(rresiduals) <- names(tmp)
    results <- list(
      model_results = lm_results,
      coef_data = coef_pval_mat$coef,
      pval_data = coef_pval_mat$pval,
      residuals = rresiduals)
    #rm(list=ls()[ls()!="results"])
    #gc(full = T)
    return(results)
}
