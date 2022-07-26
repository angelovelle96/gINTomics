
lm_singlegene <- function( response_var,
                            covariates,
                            interactions,
                            threads=1,
                            step=F){


    if(nrow(response_var)!=nrow(covariates)) stop(str_wrap("response_var
                and covariates should have the same number of samples"))
    if(!identical(rownames(response_var), rownames(covariates))) warning(
    str_wrap("response_var and covariates have different sample names,
                    assuming that samples in response_var and covariates are
                    in the same order"))
    if(length(intersect(colnames(response_var), colnames(covariates)))>0){
        warning(str_wrap("response_var and covariates have common colnames,
                         adding '_cov' to covariates colnames"))
        colnames(covariates) <- paste0(colnames(covariates), "_cov")
        data <- cbind(response_var, covariates)
        interactions <- lapply(interactions, function(x)
          paste0(x, "_cov"))
    } else {
        data <- cbind(response_var, covariates)
    }

    interactions <- lapply(interactions, function(x)
      intersect(x, colnames(covariates)))


    lm_results <-  mclapply(1:length(interactions), function(y){

        tmp <- interactions[[y]]
        tmp <- gsub("-", "_", tmp)
        tmp <- gsub(";", "_", tmp)
        tmp <- gsub(":", "_", tmp)
        tmp <- gsub("\\*", "_", tmp)
        tmp <- gsub("%in%", "_", tmp)
        tmp <- gsub("\\^", "_", tmp)
        tmp2 <- gsub("-", "_", names(interactions)[y])
        tmp2 <- gsub(";", "_", tmp2)
        tmp2 <- gsub(":", "_", tmp2)
        tmp2 <- gsub("\\*", "_", tmp2)
        tmp2 <- gsub("%in%", "_", tmp2)
        tmp2 <- gsub("\\^", "_", tmp2)
        tmp <- formula(paste0(tmp2,"~", paste0(tmp, collapse  = "+")))
        tmp2 <- data
        colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
        colnames(tmp2) <- gsub(";", "_", colnames(tmp2))
        colnames(tmp2) <- gsub(":", "_", colnames(tmp2))
        colnames(tmp2) <- gsub("\\*", "_", colnames(tmp2))
        colnames(tmp2) <- gsub("%in%", "_", colnames(tmp2))
        colnames(tmp2) <- gsub("\\^", "_", colnames(tmp2))

        if(step==T){
            lm_results <- summary(step(lm(tmp, data = tmp2)))
        } else {
            lm_results <- summary(lm(tmp, data = tmp2))
        }
        return(lm_results)

    }, mc.cores = threads)

    names(lm_results) <- names(interactions)
    return(lm_results)

}



