run_lm_integration <- function(response_var,
                               covariates,
                               interactions,
                               threads=1,
                               step=F,
                               cnv_mode=F){

    if(cnv_mode==T){
      interactions <- as.list(intersect(colnames(covariates),
                                        colnames(response_var)))
      names(interactions) <- unlist(interactions)
    }
    tmp <- data_check(response_var = response_var,
                      interactions = interactions)
    response_var <- tmp$response_var
    interactions <- tmp$interactions

    lm_results <- lm_singlegene(response_var = response_var,
                                covariates = covariates,
                                interactions = interactions,
                                threads = threads,
                                step = step)

    return(lm_results)
}
