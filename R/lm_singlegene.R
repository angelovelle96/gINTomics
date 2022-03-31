
lasso_model_results_can <- function(data, covariates){

  pbmclapply(1:length(covariates), function(y){


  covariates <- expanded_tf_mirna_couples_can[[y]][,"tf"]
  tmp <- gsub("-", "_", tf)
  tmp2 <- gsub("-", "_", names(expanded_tf_mirna_couples_can)[y])
  tmp2 <- gsub(";", "_", tmp2)
  tmp2 <- gsub(":", "_", tmp2)
  tmp <- formula(paste0(tmp2,"~", paste0(tmp, collapse  = "+")))
  tmp2 <- data
  colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
  colnames(tmp2) <- gsub(";", "_", colnames(tmp2))
  colnames(tmp2) <- gsub(":", "_", colnames(tmp2))

  #lasso2_model <- summary(l1ce(tmp, data = tmp2))

  lm_results <- summary(step(lm(tmp, data = tmp2)))
  #rownames(lm_results$coefficients) <- tf

  return(lm_results)

}, mc.cores = 22)

}



