run_singlegene_edgeR_integration <- function(response_var, covariates, design_mat_allgene=NULL, design_mat_singlegene=NULL, threads=1){


  tmp <- intersect(rownames(response_var), rownames(covariates))
  covariates <- drop.levels(covariates[tmp, ])
  response_var <- as.matrix(t(response_var[tmp, ]+1))

  y_all <- allgene_edgeR_model(response_var = response_var,  covariates = covariates, design_mat_allgene = design_mat_allgene)
  y_gene <- singlegene_edgeR_model(response_var = response_var, covariates = covariates, y_all = y_all, design_mat_singlegene =  design_mat_singlegene, threads = threads)
  top_list <- extracting_edgeR_coef(y_gene, threads = threads)

  return(top_list)
}



##### Offsets dei modelli all e single gene quando non abbiamo RUV calcnormfactors per allgene
