#' @import plyr
building_result_matrices <- function(model_results,
                                           type){

    if(type=="edgeR"){
      tmp <- lapply(model_results, function(x) as.data.frame( x["coef",],
          check.names = F))
    }
    if(type=="lm"){
      tmp <- lapply(model_results, function(x){
        data.frame(t(x[["coefficients"]][, "Estimate"]),check.names = F)
      })
    }

    coef_matrix <- rbind.fill(tmp)
    rownames(coef_matrix) <- names(model_results)
    for(i in 1:ncol(coef_matrix)) {
        tmp <- coef_matrix[, i]
        tmp[sapply(tmp, is.null)] <- NA
        tmp <- unlist(tmp)
        coef_matrix[, i] <- tmp
    }



    if(type=="edgeR"){
      tmp <- lapply(model_results, function(x) as.data.frame( x["PValue",],
        check.names = F))
    }
    if(type=="lm"){
      tmp <- lapply(model_results, function(x){
        data.frame(t(x[["coefficients"]][, "Pr(>|t|)"]),check.names = F)
      })
    }

    pval_matrix <- rbind.fill(tmp)
    rownames(pval_matrix) <- names(model_results)
    for(i in 1:ncol(pval_matrix)) {
        tmp <- pval_matrix[, i]
        tmp[sapply(tmp, is.null)] <- NA
        tmp <- unlist(tmp)
        pval_matrix[, i] <- tmp
    }

    coef_pval_mat <- list()
    coef_pval_mat[["coef"]] <- coef_matrix
    coef_pval_mat[["pval"]] <- pval_matrix
    return(coef_pval_mat)
}
