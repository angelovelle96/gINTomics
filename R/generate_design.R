generate_design <- function(response_var,
                            covariates,
                            interactions = NULL,
                            steady_covariates = NULL,
                            cnv_mode = F,
                            reference=NULL) {

    response_var <- t(response_var+1)
    if(is.null(covariates)) stop(str_wrap("covariates should be provided if
                                design_mat_singlegene is not provided"))
    if(ncol(response_var)!=nrow(covariates)) stop(str_wrap("response_var
            and covariates should have the same number of samples"))
    if(!identical(colnames(response_var), rownames(covariates))) warning(
        str_wrap("response_var and covariates have different sample names,
                    assuming that samples in response_var and covariates are
                    in the same order"))
    tmp <- as.data.frame(covariates)
    tmp <- tmp[sapply(tmp, function(x) !is.numeric(x))]
    if(ncol(tmp)>0){
        tmp <- sapply(tmp, is.factor)
        if(sum(tmp)!=length(tmp)) stop(str_wrap("covariates should be provided
                                                  as numeric or factors"))
    }


    if(cnv_mode==F){

        if(is.null(interactions)) stop(str_wrap("You should provide interactions
                if design_mat_singlegene is not provided and cnv_mode=F"))
        if(nrow(response_var)!=length(interactions)) stop(str_wrap("response_var
            and interactions should have the same number of genes"))
        if(!identical(rownames(response_var), names(interactions))) warning(
            str_wrap("response_var and interactions have different gene names,
                    assuming that genes in response_var and interactions are
                    in the same order"))

        interactions <- lapply(interactions, function(x)
            intersect(x, colnames(covariates)))
        design_matrix <- lapply(1:length(interactions), function(x) {
            des_cov <- covariates_check(x=x,
                                       response_var=response_var,
                                       covariates=covariates,
                                       steady_covariates=steady_covariates,
                                       reference=reference,
                                       interactions = interactions,
                                       cnv_mode = cnv_mode)
            design <- des_cov$design
            colnames(design)[2:(ncol(design)-length(steady_covariates))
                             ] <- interactions[[x]]
            return(design)
        })
        names(design_matrix) <- rownames(response_var)
    } else {
        message("Using cnv_mode")
        design_matrix <- lapply(1:nrow(response_var), function(x) {

            des_cov <- covariates_check(x=x,
                             response_var=response_var,
                             covariates=covariates,
                             steady_covariates=steady_covariates,
                             reference=reference,
                             cnv_mode=cnv_mode,
                             interactions=interactions)
            design <- des_cov$design
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
        })
        names(design_matrix) <- rownames(response_var)
    }
    return(design_matrix)
}



