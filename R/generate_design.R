generate_design <- function(response_var,
                            covariates,
                            interactions = NULL,
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


    if(cnv_mode==F){


        design_matrix <- mclapply(1:length(interactions), function(x) {
            des_cov <- covariates_check(x=x,
                                       response_var=response_var,
                                       covariates=covariates,
                                       steady_covariates=steady_covariates,
                                       reference=reference,
                                       interactions = interactions,
                                       cnv_mode = cnv_mode)

            design <- model.matrix(des_cov$formula, data=des_cov$des_data)
            colnames(design)[2:(ncol(design)-length(steady_covariates))
                             ] <- interactions[[x]]
            return(design)
        }, mc.cores = threads)
        names(design_matrix) <- rownames(response_var)
    } else {
        message("Using cnv_mode")
        design_matrix <- mclapply(1:nrow(response_var), function(x) {

            des_cov <- covariates_check(x=x,
                             response_var=response_var,
                             covariates=covariates,
                             steady_covariates=steady_covariates,
                             reference=reference,
                             cnv_mode=cnv_mode,
                             interactions=interactions)
            design <- model.matrix(des_cov$formula, data=des_cov$des_data)
            cov <- des_cov$cov
            tmp <- des_cov$formula

            if(is.numeric(covariates[,rownames(response_var)[x]])){
              colnames(design)[2:ncol(design)] <- cov
              colnames(design)[2] <- "cnv"
            }else{
              tmp <- strsplit(split = " ", as.character(tmp))[[2]][1]
              colnames(design) <- gsub(tmp, "cnv_", colnames(design))
            }
            return(design)
        }, mc.cores = threads)
        names(design_matrix) <- rownames(response_var)
    }
    return(design_matrix)
}



