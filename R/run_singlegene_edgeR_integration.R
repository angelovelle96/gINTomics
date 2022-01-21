run_singlegene_edgeR_integration <- function(response_var, covariates, design_mat_allgene=NULL, design_mat_singlegene=NULL, offset_allgene=NULL, offset_singlegene=NULL, threads=1){

  if(is.atomic(response_var)){
    message("response_var is an atomic vector, converting to matrix")
    tmp <- as.matrix(response_var)
    rownames(tmp) <- names(response_var)
    check_names <-  names(response_var)
    response_var <- tmp
  }else{
    check_names <-  rownames(response_var)
  }

  if(is.atomic(covariates)){
    message("covariates is an atomic vector, converting to data.frame")
    tmp <- as.data.frame(covariates)
    rownames(tmp) <- names(covariates)
    check_names2 <- names(covariates)
    covariates <- tmp
  }else{
    check_names2 <-  rownames(covariates)
  }
  if(!is.data.frame(covariates)) stop("covariates should be a data.frame or an atomic vector")


  if(!identical(check_names, check_names2)){
    message("Sample names in response_var e covariates seems to differ, assuming that response_var and covariates have the same sample order")
  }

  response_var <- as.matrix(t(response_var+1))#######lasciare+1?

  if(!is.matrix(response_var)) stop("response_var should be a data.frame, a matrix or an atomic vector")
  if(ncol(response_var)!=nrow(covariates)) stop("response_var and covariates have a different number of samples")


  y_all <- allgene_edgeR_model(response_var = response_var,  covariates = covariates, design_mat_allgene = design_mat_allgene, offset_allgene = offset_allgene)
  y_gene <- singlegene_edgeR_model(response_var = response_var, covariates = covariates, y_all = y_all, design_mat_singlegene =  design_mat_singlegene, offset_singlegene = offset_singlegene, threads = threads)
  model_results <- edger_coef_test(y_gene, threads = threads)
  coef_pval_mat <- building_edger_result_matrices(model_results = model_results)
  results <- list(model_results=model_results, coef_data=coef_pval_mat$coef, pval_data=coef_pval_mat$pval)

  return(results)
}



##### Offsets dei modelli all e single gene quando non abbiamo RUV calcnormfactors per allgene
