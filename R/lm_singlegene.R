
lm_singlegene <- function(  response_var,
                            covariates,
                            interactions,
                            steady_covariates=NULL,
                            reference=NULL,
                            cnv_mode=F,
                            threads=1,
                            step=F){



    if(length(intersect(colnames(response_var), colnames(covariates)))>0){
        message(str_wrap("response_var and covariates have common colnames,
                         adding '_cov' to covariates colnames"))
        colnames(covariates) <- paste0(colnames(covariates), "_cov")
        interactions <- lapply(interactions, function(x)
          paste0(x, "_cov"))
    }

    lm_results <-  mclapply(1:length(interactions), function(x){

        tmp <- covariates_check(x = x,
                                response_var = t(response_var),
                                covariates = covariates,
                                interactions = interactions,
                                steady_covariates = steady_covariates,
                                reference = reference,
                                cnv_mode = cnv_mode,
                                linear=T)
        if(step==T){
            lm_results <- summary(step(lm(tmp$formula, data = tmp$des_data),
                                       trace = 0))
        } else {
            lm_results <- summary(lm(tmp$formula, data = tmp$des_data))
        }
        return(lm_results)

    }, mc.cores = threads)

    names(lm_results) <- names(interactions)
    return(lm_results)

}



