generate_design <- function(response_var,
                            covariates,
                            interactions = NULL,
                            steady_covariates = NULL,
                            cnv_mode = F) {

  response_var <- t(response_var+1)
  if(cnv_mode==F){

    if(nrow(response_var)!=length(interactions)) stop(str_wrap("response_var
            and interactions should have the same number of genes"))
    if(!identical(rownames(response_var), names(interactions))) warning(
      str_wrap("response_var and interactions have different gene names,
                    assuming that genes in response_var and interactions are
                    in the same order"))

    interactions <- lapply(interactions, function(x)
        intersect(x, colnames(covariates)))
    design_matrix <- lapply(1:length(interactions), function(x) {

      cov <- interactions[[x]]
      if(!is.null(steady_covariates)) cov <- c(cov, steady_covariates)
      tmp <- gsub("-", "_", cov)
      tmp <- formula(paste0("~", paste0(tmp, collapse = "+")))
      tmp2 <- covariates
      colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
      design <- model.matrix(tmp, data = tmp2)
      colnames(design)[2:ncol(design)] <- cov
      return(design)
    })
    names(design_matrix) <- names(interactions)
  } else {
    message("Using cnv_mode")
    design_matrix <- lapply(1:nrow(response_var), function(x) {

      cov <- rownames(response_var)[x]
      if(!is.null(steady_covariates)) cov <- c(cov, steady_covariates)
      tmp <- gsub("-", "_", cov)
      tmp <- gsub(";", "_", tmp)
      tmp <- formula(paste0("~", paste0(tmp, collapse = "+")))
      tmp2 <- covariates
      colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
      colnames(tmp2) <- gsub(";", "_", colnames(tmp2))
      design <-  model.matrix(tmp, data = tmp2)
      colnames(design)[2:ncol(design)] <- cov
      colnames(design)[2] <- "cnv"
      return(design)
    })
    names(design_matrix) <- rownames(response_var)
  }
    return(design_matrix)
}



