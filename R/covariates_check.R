covariates_check <- function(x,
                             response_var,
                             covariates,
                             interactions=NULL,
                             steady_covariates=NULL,
                             reference=NULL,
                             cnv_mode=F){

    if(cnv_mode==F){
      cov <- interactions[[x]]
    }else{
      cov <- rownames(response_var)[x]
    }
    if(!is.null(steady_covariates)) cov <- c(cov, steady_covariates)
    tmp <- gsub("-", "_", cov)
    tmp <- gsub(";", "_", tmp)
    tmp <- gsub(":", "_", tmp)
    tmp <- gsub("\\*", "_", tmp)
    tmp <- gsub("%in%", "_", tmp)
    tmp <- gsub("\\^", "_", tmp)
    tmp2 <- formula(paste0("~", paste0(tmp, collapse = "+")))
    tmp3 <- covariates
    colnames(tmp3) <- gsub("-", "_", colnames(tmp3))
    colnames(tmp3) <- gsub(";", "_", colnames(tmp3))
    colnames(tmp3) <- gsub(":", "_", colnames(tmp3))
    colnames(tmp3) <- gsub("\\*", "_", colnames(tmp3))
    colnames(tmp3) <- gsub("%in%", "_", colnames(tmp3))
    colnames(tmp3) <- gsub("\\^", "_", colnames(tmp3))
    tmp3 <- tmp3[, intersect(tmp, colnames(tmp3))]
    if(!is.null(reference)){
      tmp3[sapply(tmp3, is.factor)] <- lapply(tmp3[sapply(tmp3, is.factor)],
                                            function(x) relevel(x,
                                            ref = reference))
    }

    tmp3 <- as.data.frame(tmp3)
    colnames(tmp3) <- tmp
    design <-  model.matrix(tmp2, data = tmp3)
    return(list(design=design, cov=cov, formula=tmp2))

}
