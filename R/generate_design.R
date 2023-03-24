generate_design <- function(response_var,
                            covariates,
                            interactions,
                            steady_covariates = NULL,
                            reference=NULL,
                            BPPARAM=BPPARAM) {

    response_var <- t(response_var+1)
    tmp <- as.data.frame(covariates)
    tmp <- tmp[sapply(tmp, function(x) !is.numeric(x))]
    if(ncol(tmp)>0){
        tmp <- sapply(tmp, is.factor)
        if(sum(tmp)!=length(tmp)) stop(str_wrap("covariates should be provided
                                                  as numeric or factors"))
    }

    tmp <- unlist(lapply(interactions, length))
    single_cov=F
    if(sum(tmp==1)==length(tmp)){
      single_cov=T
    }

    design_matrix <- bplapply(1:length(interactions), function(x) {
      des_cov <- covariates_check(x=x,
                                  response_var=response_var,
                                  covariates=covariates,
                                  steady_covariates=steady_covariates,
                                  reference=reference,
                                  interactions = interactions,
                                  single_cov=single_cov)

      design <- model.matrix(des_cov$formula, data=des_cov$des_data)
      tmp <- unlist(lapply(interactions, length))
      if(!"cov"%in%colnames(design)){
        colnames(design)[2:(ncol(design)-length(steady_covariates))
                         ] <- interactions[[x]]
      }
      return(design)
    }, BPPARAM = BPPARAM)
    names(design_matrix) <- names(interactions)

    return(design_matrix)

    }




