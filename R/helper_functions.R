#' @export
data_check <- function( response_var,
                        covariates,
                        interactions){

    if(identical(interactions,"auto")){
      message("Generating interactions")
      interactions <- as.list(intersect(colnames(covariates),
                                        colnames(response_var)))
      names(interactions) <- unlist(interactions)
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


    interactions <- lapply(interactions, function(x)
      intersect(x, colnames(covariates)))
    tmp <- lapply(interactions, length)
    tmp <- which(tmp!=0)
    interactions <- lapply(tmp, function(x) interactions[[x]])
    tmp <- intersect(names(interactions), colnames(response_var))
    if(length(tmp)==0) stop(str_wrap("No genes left in common between
                                      response_var and interactions"))
    response_var <- response_var[,tmp]
    interactions <- lapply(tmp, function(x) interactions[[x]])
    names(interactions) <- tmp

    rresult <- list(response_var=response_var,
                    covariates=covariates,
                    interactions=interactions)
    return(rresult)
}

##################################################################
#' @export
covariates_check <- function(response_var,
                             covariates,
                             interactions,
                             steady_covariates=NULL,
                             linear=F,
                             reference=NULL){

  if(length(intersect(colnames(response_var), colnames(covariates)))>0){
    message(str_wrap("response_var and covariates have common colnames,
                           adding '_cov' to covariates colnames"))
    colnames(covariates) <- paste0(colnames(covariates), "_cov")
    interactions <- lapply(interactions, function(x)
      paste0(x, "_cov"))
    if(!is.null(steady_covariates)){
      steady_covariates <- paste0(steady_covariates, "_cov")
    }
  }
  if(!is.null(steady_covariates)){
    interactions <- lapply(interactions, function(x)
      unique(c(x,steady_covariates)))
  }
  bad_str <- paste0(c("-",";",":","\\*","%in%","\\^"), collapse = "|")
  original_id <- unique(c(unlist(interactions), names(interactions)))
  original_id <- data.frame(original=original_id,
                            tranformed=gsub(bad_str, "_", original_id))
  tmp <- apply(original_id, 1, function(x) x[1]!=x[2])
  original_id <- original_id[tmp,]
  interactions <- lapply(interactions, function(x) gsub(bad_str, "_", x))
  colnames(covariates) <- gsub(bad_str, "_", colnames(covariates))
  covariates <- covariates[, colnames(covariates)%in%unlist(interactions)]
  if(linear==T){
    colnames(response_var) <- gsub(bad_str, "_", colnames(response_var))
    names(interactions) <- gsub(bad_str, "_", names(interactions))
  }
  if(!is.null(reference)){
    covariates[sapply(covariates, is.factor)] <- lapply(
      covariates[sapply(covariates, is.factor)],
      function(x) relevel(x, ref = reference))
  }
  return(list(covariates=covariates,
              response_var=response_var,
              interactions=interactions,
              original_id=original_id,
              steady_covariates=steady_covariates))
}


######################################################
#' Function for formula generation

generate_formula <- function(interactions,
                             linear=F){

  fformula <- lapply(seq_along(interactions), function(x){
    cov <- interactions[[x]]
    ans <- formula(paste0("~", paste0(cov, collapse = "+")))
    if(linear==T){
      ans <- formula(paste0(names(interactions)[x],
                            "~", as.character(ans)[2]))
    }
    return(ans)
  })
  names(fformula) <- names(interactions)
  return(fformula)
}

#####################################################
#' @import plyr
#' @export

building_result_matrices <- function(model_results,
                                     type,
                                     single_cov=F){
  ####coef extraction
  if(type=="edgeR"){
    tmp <- lapply(model_results, function(x) as.data.frame( x["coef",],
                                                            check.names = F))
  }
  if(type=="lm"){
    tmp <- lapply(model_results, function(x){
      tmp <- data.frame(t(x[["coefficients"]][, "Estimate"]),check.names = F)
      colnames(tmp) <- rownames(x[["coefficients"]])
      return(tmp)
    })
  }
  if(single_cov==T){
    tmp <- lapply(tmp, function(x){
      ans=x
      pos <- length(grep("Intercept", names(ans)))+1
      names(ans)[pos] <- "cov"
      return(ans)
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

  ######pvalue extraction
  if(type=="edgeR"){
    tmp <- lapply(model_results, function(x) as.data.frame( x["PValue",],
                                                            check.names = F))
  }
  if(type=="lm"){
    tmp <- lapply(model_results, function(x){
      tmp <- data.frame(t(x[["coefficients"]][, "Pr(>|t|)"]),check.names = F)
      colnames(tmp) <- rownames(x[["coefficients"]])
      return(tmp)
    })
  }
  if(single_cov==T){
    tmp <- lapply(tmp, function(x){
      ans=x
      pos <- length(grep("Intercept", names(ans)))+1
      names(ans)[pos] <- "cov"
      return(ans)
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


################################################
#' @export

create_multiassay <- function(methylation=NULL,
                              cnv_data=NULL,
                              gene_exp=NULL,
                              miRNA_exp=NULL,
                              miRNA_cnv_data=NULL,
                              ...){

  if((is.null(methylation)+is.null(cnv_data)+
      is.null(gene_exp)+is.null(miRNA_exp)+is.null(miRNA_cnv_data))>=4){
    stop(str_wrap("You need to provide at least two assays or
                  SummarizedExperiment"))
  }
  eexperiments <- Filter(Negate(is.null),
                         list(gene_exp=gene_exp,
                              methylation=methylation,
                              cnv_data=cnv_data,
                              miRNA_exp=miRNA_exp,
                              miRNA_cnv_data=miRNA_cnv_data))
  mmultiassay <- MultiAssayExperiment(experiments = eexperiments)
  mmultiassay <- intersectColumns(mmultiassay)
  return(mmultiassay)
}


######################################################
#' Data normalization for linear models


    data_norm <- function(data,
                          method="TMM"){

      data <- edgeR::DGEList(t(data))
      data <- edgeR::calcNormFactors(data, method=method)
      data <- t(edgeR::cpm(data, log = F))
      return(data)
    }

#########################################################
#' conversion back to the orginal IDs
#' @importFrom stringi stri_replace_all_regex

    id_conversion <- function(dictionary,
                              results){

      dictionary$tranformed <- paste0("^", dictionary$tranformed, "$")

      names(results$model_results) <- stri_replace_all_regex(
        names(results$model_results),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      rownames(results$residuals) <- stri_replace_all_regex(
        rownames(results$residuals),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      colnames(results$coef_data) <- stri_replace_all_regex(
        colnames(results$coef_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      rownames(results$coef_data) <- stri_replace_all_regex(
        rownames(results$coef_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      colnames(results$pval_data) <- stri_replace_all_regex(
        colnames(results$pval_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      rownames(results$pval_data) <- stri_replace_all_regex(
        rownames(results$pval_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

       return(results)
    }




