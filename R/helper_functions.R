
#' @importFrom stats sd

.data_check <- function( response_var,
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
    check_sd <- apply(response_var, 2, function(x) sd(x, na.rm = T))
    check_sd[is.na(check_sd)] <- 0
    if(sum(check_sd==0)>0) {
      message(str_wrap("removing response variables with
                                zero standard deviation"))
      response_var <- response_var[, check_sd>0]
    }

    if(is.numeric(covariates[,1])){
      check_sd <- apply(covariates, 2, function(x) sd(x, na.rm = T))
      check_sd[is.na(check_sd)] <- 0
      if(sum(check_sd==0)>0) {
        message(str_wrap("removing covariates with
                                  zero standard deviation"))
        covariates <- covariates[, check_sd>0]
      }
    }else{
      check_sd <- apply(covariates, 2, function(x){
        ans <- length(unique(x))
        return(ans)
      })
      if(sum(check_sd==1)>0) {
        message(str_wrap("removing covariates with
                                  zero standard deviation"))
        covariates <- covariates[, check_sd>1]
      }
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


#' @importFrom stats relevel

.covariates_check <- function(response_var,
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


#' @importFrom stats formula

.generate_formula <- function(interactions,
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
#' @importFrom plyr rbind.fill

.building_result_matrices <- function(model_results,
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
      if(length(ans)>=pos){
        names(ans)[pos] <- "cov"
      }
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
      if(length(ans)>=pos){
        names(ans)[pos] <- "cov"
      }
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
#' MultiAssayExperiment generation
#' @description
#' This function will generate a proper MultiAssayExperiment suitable for the
#' **run_multiomics** function.
#' @param methylation Matrix or SummarizedExperiment for Methylation data
#' @param cnv_data Matrix or SummarizedExperiment for genes' Copy Number
#' Variation data
#' @param gene_exp Matrix or SummarizedExperiment for Gene expression data
#' @param miRNA_exp Matrix or SummarizedExperiment for miRNA expression data
#' @param miRNA_cnv_data Matrix or SummarizedExperiment for miRNA's Copy
#' Number Variations data
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
  return(mmultiassay)
}


######################################################

    .data_norm <- function(data,
                          method="TMM"){

      data <- edgeR::DGEList(t(data))
      data <- edgeR::calcNormFactors(data, method=method)
      data <- t(edgeR::cpm(data, log = F))
      return(data)
    }

#########################################################
#' @importFrom stringi stri_replace_all_regex

    .id_conversion <- function(dictionary,
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

      colnames(results$data$response_var) <- stri_replace_all_regex(
        colnames(results$data$response_var),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      colnames(results$data$covariates) <- stri_replace_all_regex(
        colnames(results$data$covariates),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)
       return(results)
    }



#########################################################
#' @importFrom randomForest randomForest importance

.rf_selection <- function(data,
                          formula){

  if(length(attr(terms(formula), "term.labels"))>
     as.integer(nrow(data)*0.4)){
    ans <- randomForest(formula, data)
    ans <- importance(ans)
    ans <- ans[order(ans[,1], decreasing = T),]
    tmp <- as.integer(nrow(data)*0.4)
    ans <- ans[1:tmp]
    tmp <- paste(names(ans), collapse = "+")
    tmp <- paste(as.character(formula[2]), "~", tmp)
    ans <- as.formula(tmp)
  }else{
      ans <- formula
  }

  ans <- as.formula(paste(as.character(ans[1]), as.character(ans[3])))
  return(ans)
}



####################################################
#' @importFrom reshape2 melt
#' @importFrom stats IQR quantile

setMethod("extract_model_res", "list",
          function(model_results,
                   outliers=F,
                   species="hsa",
                   filters="hgnc_symbol",
                   genes_info=NULL){

            if(names(model_results)[length(model_results)]!="data"){
              tmp <- names(model_results)
              data <- lapply(tmp, function(x){
               ans <- extract_model_res(model_results = model_results[[x]],
                                        outliers=outliers,
                                        species=species,
                                        filters=filters,
                                        genes_info=genes_info)
               ans <- cbind(ans, class=x)
               return(ans)
              })
              names(data) <- tmp
              data <- rbind.fill(data)
              return(data)
            }

            data <- cbind(response=rownames(model_results$coef_data),
                          model_results$coef_data)
            data <- melt(data,
                         id.vars="response",
                         variable.name = "cov")
            rownames(data) <- paste0(data$response, "_", data$cov)
            data <- data[!is.na(data$value),]
            pval <- cbind(response=rownames(model_results$pval_data),
                          model_results$pval_data)
            pval <- reshape2::melt(pval,
                                   id.vars="response",
                                   variable.name = "cov")
            rownames(pval) <- paste0(pval$response, "_", pval$cov)
            pval <- pval[rownames(data),]
            tmp <- rep("not_significant", nrow(pval))
            tmp[pval$value<=0.05] <- "significant"
            names(tmp) <- rownames(pval)
            data$pval <- pval[rownames(data), "value"]
            data$significativity <- tmp[rownames(data)]
            data$sign <- rep("negative", nrow(data))
            data$sign[data$value>0]="positive"
            data$cov <- as.character(data$cov)
            data$cov <- gsub("_cov", "", data$cov)
            data$significativity[data$significativity=="significant" &
                                   data$cov=="cnv"] <- "significant_cnv"
            data$significativity[data$significativity=="significant" &
                                   data$cov=="met"] <- "significant_met"
            data$cov[data$cov%in%c("cov", "cnv", "met")] <-
              data$response[data$cov%in%c("cov", "cnv", "met")]

            if(is.null(genes_info)){
              genes_info <- .download_gene_info(c(data$cov, data$response),
                                               filters=filters,
                                               species = species)
              }
            tmp <- intersect(data$cov, rownames(genes_info))
            tmp <- genes_info[tmp,]
            tmp2 <- intersect(data$cov, genes_info$ensembl_gene_id)
            tmp2 <- genes_info[genes_info$ensembl_gene_id%in%tmp2,]
            tmp2 <- tmp2[!duplicated(tmp2$ensembl_gene_id),]
            rownames(tmp2) <- tmp2$ensembl_gene_id
            tmp <- rbind(tmp, tmp2)
            data$chr_cov <- tmp[as.character(data$cov),
                                       "chromosome_name"]
            data$cytoband_cov <- tmp[as.character(data$cov), "band"]
            data$start_cov <- tmp[as.character(data$cov), "start_position"]
            data$end_cov <- tmp[as.character(data$cov), "end_position"]

            tmp <- intersect(data$response, rownames(genes_info))
            tmp <- genes_info[tmp,]
            tmp2 <- intersect(data$response, genes_info$ensembl_gene_id)
            tmp2 <- genes_info[genes_info$ensembl_gene_id%in%tmp2,]
            tmp2 <- tmp2[!duplicated(tmp2$ensembl_gene_id),]
            rownames(tmp2) <- tmp2$ensembl_gene_id
            tmp <- rbind(tmp, tmp2)
            data$chr_response <- tmp[as.character(data$response),
                                     "chromosome_name"]
            data$cytoband_response <- tmp[as.character(data$response),
                                          "band"]
            data$start_response <- tmp[as.character(data$response),
                                       "start_position"]
            data$end_response <- tmp[as.character(data$response),
                                     "end_position"]
            colnames(data)[colnames(data)=="value"] <- "coef"
            tmp <-  apply(model_results$data$response_var, 2, mean)
            data$response_value <- tmp[data$response]
            tmp <-  apply(model_results$data$covariates, 2, mean)
            names(tmp) <- gsub("_cov", "", names(tmp))
            names(tmp) <- gsub("_cnv", "", names(tmp))
            names(tmp) <- gsub("_met", "", names(tmp))
            data$cov_value <- tmp[data$cov]
            mmin <- quantile(data$coef, 0.25) - 1.5*IQR(data$coef)
            mmax <- quantile(data$coef, 0.75) + 1.5*IQR(data$coef)

            if(outliers==F){
              data <- data[data$coef>mmin,]
              data <- data[data$coef<mmax,]
            }

            return(data)

          }
)


####################################################
#' @importFrom plyr rbind.fill

setMethod("extract_model_res", "MultiOmics",
          function(model_results,
                   outliers=F,
                   species="hsa",
                   filters="hgnc_symbol",
                   genes_info=NULL){

            tmp <- model_results
            if(is.null(genes_info)){
                  if(names(model_results[[1]])[
                    length(model_results[[1]])]!="data"){
                    tmp <- unlist(tmp, recursive = F)
                  }
                  tmp <- lapply(tmp, function(x){
                      ans <- cbind(response=rownames(x$coef_data),
                                    x$coef_data)
                      ans <- melt(ans,
                                  id.vars="response",
                                  variable.name = "cov")
                      ans <- ans[!is.na(ans$value),]
                      return(ans)
                  })

              tmp <- rbind.fill(tmp)
              tmp$cov <- as.character(tmp$cov)
              tmp$cov <- gsub("_cov", "", tmp$cov)
              tmp$cov[tmp$cov%in%c("cov", "cnv", "met")] <-
                tmp$response[tmp$cov%in%c("cov", "cnv", "met")]

              genes_info <- .download_gene_info(c(tmp$cov, tmp$response),
                                                filters=filters,
                                                species = species)
            }

            tmp <- lapply(model_results, function(x)
              extract_model_res(x,
                                outliers=outliers,
                                species=species,
                                genes_info=genes_info))
            data <- lapply(names(tmp), function(x) cbind(tmp[[x]], omics=x))
            names(data) <- names(tmp)
            data <- rbind.fill(data)
            return(data)
          }
)


#######################################################
#DA ELIMINARE
setMethod("extract_data", "list",
          function(model_results,
                   species="hsa",
                   filters="hgnc_symbol"){

    res_layer <- model_results$data$response_var
    res_layer <- apply(res_layer, 2, mean)
    cov_layer <- model_results$data$covariates[
                    rownames(model_results$data$response_var),]
    colnames(cov_layer) <- gsub("_cov", "", colnames(cov_layer))
    cov_layer <- apply(cov_layer, 2, mean)
    coef_layer <- NULL
    tmp <- as.data.frame(model_results$coef_data)
    tmp <- tmp[, colnames(tmp)!="(Intercept)", drop=F]
    if(identical(colnames(tmp),"cov")){
      coef_layer <- tmp
      coef_layer <- as.data.frame(cbind(coef_layer,
                                        pval=model_results$pval_data$cov))
    }
    genes_info <- .download_gene_info(unique(c(names(res_layer),
                                              names(cov_layer))),
                                     filters=filters,
                                     species = species)
    tmp <- genes_info$ensembl_gene_id%in%c(names(res_layer), names(cov_layer))
    tmp <- genes_info[tmp,]
    tmp <- tmp[!duplicated(tmp$ensembl_gene_id),]
    rownames(tmp) <- tmp$ensembl_gene_id
    genes_info <- rbind(genes_info, tmp)
    res_layer <- data.frame(mean=res_layer, genes_info[names(res_layer),])
    tmp <- is.na(res_layer$chromosome_name)+
      is.na(res_layer$start_position)+
      is.na(res_layer$end_position)
    res_layer <- res_layer[tmp==0,]
    cov_layer <- data.frame(mean=cov_layer, genes_info[names(cov_layer),])
    tmp <- is.na(cov_layer$chromosome_name)+
      is.na(cov_layer$start_position)+
      is.na(cov_layer$end_position)
    cov_layer <- cov_layer[tmp==0,]
    if(!is.null(coef_layer)){
      coef_layer <- as.data.frame(cbind(coef_layer,
                                        genes_info[rownames(coef_layer),]))
      tmp <- is.na(coef_layer$chromosome_name)+
        is.na(coef_layer$start_position)+
        is.na(coef_layer$end_position)
      coef_layer <- coef_layer[tmp==0,]
    }

    ans <- Filter(Negate(is.null),
                  list(res_layer=res_layer,
                       cov_layer=cov_layer,
                       coef_layer=coef_layer))

    return(ans)
})



###############################################
#DA ELIMINARE
setMethod("extract_data", "MultiOmics",
          function(model_results,
                   species="hsa"){

    ans <- lapply(model_results, function(x)
      extract_data(x,species=species))
    return(ans)

})




