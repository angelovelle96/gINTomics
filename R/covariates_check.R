covariates_check <- function(x,
                             response_var,
                             covariates,
                             interactions=NULL,
                             steady_covariates=NULL,
                             reference=NULL,
                             cnv_mode=F,
                             linear=F){

    if(cnv_mode==F | linear==T){
      cov <- interactions[[x]]
    }else{
      cov <- rownames(response_var)[x]
    }
    if(!is.null(steady_covariates)) cov <- c(cov, steady_covariates)
    bad_str <- paste0(c("-",";",":","\\*","%in%","\\^"), collapse = "|")
    tmp <- gsub(bad_str, "_", cov)
    tmp2 <- formula(paste0("~", paste0(tmp, collapse = "+")))
    tmp3 <- covariates
    colnames(tmp3) <- gsub(bad_str, "_", colnames(tmp3))
    tmp3 <- tmp3[, intersect(tmp, colnames(tmp3))]
    if(!is.null(reference)){
      tmp3[sapply(tmp3, is.factor)] <- lapply(tmp3[sapply(tmp3, is.factor)],
                                            function(x) relevel(x,
                                            ref = reference))
    }

    tmp3 <- as.data.frame(tmp3)
    colnames(tmp3) <- tmp

    if(linear==T){
      tmp2 <- formula(paste0(gsub(bad_str, "", rownames(response_var)[x]),
                     "~", as.character(tmp2)[2]))
      tmp <- response_var[x,]
      tmp3 <- cbind(tmp, tmp3)
      colnames(tmp3)[1] <- gsub(bad_str, "", rownames(response_var)[x])
    }
    return(list(des_data=tmp3, cov=cov, formula=tmp2))

}
