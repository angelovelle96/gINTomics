generate_design <- function(response_var,
                            covariates,
                            interactions,
                            steady_covariates = NULL,
                            cnv_mode = F,
                            reference=NULL,
                            threads=1) {

    response_var <- t(response_var+1)
    tmp <- as.data.frame(covariates)
    tmp <- tmp[sapply(tmp, function(x) !is.numeric(x))]
    if(ncol(tmp)>0){
        tmp <- sapply(tmp, is.factor)
        if(sum(tmp)!=length(tmp)) stop(str_wrap("covariates should be provided
                                                  as numeric or factors"))
    }


    design_matrix <- mclapply(1:length(interactions), function(x) {
      des_cov <- covariates_check(x=x,
                                  response_var=response_var,
                                  covariates=covariates,
                                  steady_covariates=steady_covariates,
                                  reference=reference,
                                  interactions = interactions,
                                  cnv_mode = cnv_mode)

      design <- model.matrix(des_cov$formula, data=des_cov$des_data)
      if(cnv_mode==F){
        colnames(design)[2:(ncol(design)-length(steady_covariates))
                         ] <- interactions[[x]]
      }
      return(design)
    }, mc.cores = threads)
    names(design_matrix) <- names(interactions)

    return(design_matrix)

    }




