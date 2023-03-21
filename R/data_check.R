data_check <- function( response_var,
                        covariates,
                        interactions=NULL,
                        cnv_mode=F,
                        steady_covariates=NULL,
                        linear=F){


    if(cnv_mode==F & is.null(interactions)){
      stop(str_wrap("You should provide interactions if you are
                       not using cnv_mode"))
    }

    if(is.atomic(response_var) & is.vector(response_var)) {
      message("response_var is an atomic vector, converting to matrix")
      tmp <- as.matrix(response_var)
      rownames(tmp) <- names(response_var)
      response_var <- tmp
    }

    response_var <- as.matrix(response_var)

    if(!is.matrix(response_var)) {
      stop(str_wrap("response_var should be a data.frame,
              a matrix or an atomic vector"))
    }

    if (is.atomic(covariates)& is.vector(covariates)) {
      message("covariates is an atomic vector, converting to data.frame")
      tmp <- as.data.frame(covariates)
      rownames(tmp) <- names(covariates)
      covariates <- tmp
    }

    covariates <- as.data.frame(covariates)
    if(!is.data.frame(covariates)) {
      stop(str_wrap("covariates should be a data.frame,
              a matrix or an atomic vector"))
    }


    check_sample <- intersect(rownames(response_var), rownames(covariates))
    if(!identical(rownames(response_var), check_sample)|
       !identical(rownames(covariates), check_sample)){
      if(length(check_sample)==0){
        stop(str_wrap("No samples in common between
                      response_var and covariates"))
      }
      message(str_wrap("Aligning samples in response_var and covariates"))
      response_var <- response_var[check_sample,]
      covariates <- covariates[check_sample,]
    }
    check_sd <- apply(response_var, 2, sd)
    if(sum(check_sd==0)>0) {
      message(str_wrap("removing response variables with
                                zero standard deviation"))
      response_var <- response_var[, check_sd>0]
    }


    if(cnv_mode==F){
        tmp <- intersect(names(interactions), colnames(response_var))
        if(length(tmp)==0) stop(str_wrap("No genes in common between response_var
                                                and interactions"))
        response_var <- response_var[,tmp]
        interactions <- lapply(tmp, function(x) interactions[[x]])
        interactions <- lapply(interactions, function(x)
          intersect(x, colnames(covariates)))
        names(interactions) <- tmp
    }else{
      tmp <- intersect(colnames(response_var), colnames(covariates))
      if(length(tmp)==0) stop(str_wrap("No copy number values available
                                       for genes in response_var"))
      if(length(tmp)!=ncol(response_var)){
        message(str_wrap("removing response variables
                         without copy number values"))
      }
      response_var <- response_var[,tmp]
      covariates <- covariates[, c(tmp, steady_covariates)]
      if(linear==T){
        interactions <- lapply(tmp, function(x) interactions[[x]])
        names(interactions) <- tmp
        }

    }

    rresult <- list(response_var=response_var,
                    covariates=covariates,
                    interactions=interactions)


    return(rresult)
}
