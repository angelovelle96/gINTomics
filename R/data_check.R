data_check <- function( response_var,
                        interactions){

    check_sd <- apply(response_var, 2, sd)
    if(sum(check_sd==0)>0) {
      warning(str_wrap("removing response variables with
                                zero standard deviation"))
      response_var <- response_var[, check_sd>0]
    }
    tmp <- intersect(names(interactions), colnames(response_var))
    if(length(tmp)==0) stop(str_wrap("No genes in common between response_var
                                            and interactions"))
    response_var <- response_var[,tmp]
    interactions <- lapply(tmp, function(x) interactions[[x]])
    names(interactions) <- tmp
    rresult <- list(response_var=response_var,
                    interactions=interactions)
    return(rresult)
}
