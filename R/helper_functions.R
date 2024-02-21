
#' @importFrom stats sd

.data_check <- function(response_var,
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
      pos <- length(grep("Intercept", names(x)))+1
      if(length(x)>=pos){
        names(x)[pos] <- "cov"
      }
      return(x)
    })
  }else{
    tmp <- lapply(tmp, function(x){
      names(x) <- gsub("^.*_cnv", "cnv", names(x))
      names(x) <- gsub("^.*_met", "met", names(x))
      return(x)
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
      pos <- length(grep("Intercept", names(x)))+1
      if(length(x)>=pos){
        names(x)[pos] <- "cov"
      }
      return(x)
    })
  }else{
    tmp <- lapply(tmp, function(x){
      names(x) <- gsub("^.*_cnv", "cnv", names(x))
      names(x) <- gsub("^.*_met", "met", names(x))
      return(x)
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
#' @importFrom limma normalizeBetweenArrays
#' @importFrom edgeR DGEList calcNormFactors cpm
    .data_norm <- function(data,
                          method="TMM",
                          RNAseq=T){
      if(RNAseq){
          data <- DGEList(t(data))
          data <- calcNormFactors(data, method=method)
          data <- t(cpm(data, log = F))
      }else{
          data <- normalizeBetweenArrays(data)
          }
      return(data)
    }

#########################################################
#' @importFrom stringi stri_replace_all_regex

    .id_conversion <- function(dictionary,
                              results){

      dictionary$tranformed <- paste0("^", dictionary$tranformed, "$")

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

      colnames(results$fdr_data) <- stri_replace_all_regex(
        colnames(results$fdr_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = F)

      rownames(results$fdr_data) <- stri_replace_all_regex(
        rownames(results$fdr_data),
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
                   outliers=T,
                   species="hsa",
                   filters=c("hgnc_symbol",
                             "ensembl_gene_id",
                             "entrezgene_id"),
                   genes_info=NULL,
                   ...){

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
            pval <- melt(pval,
                         id.vars="response",
                         variable.name = "cov")
            rownames(pval) <- paste0(pval$response, "_", pval$cov)
            pval <- pval[rownames(data),]
            fdr <- cbind(response=rownames(model_results$fdr_data),
                          model_results$fdr_data)
            fdr <- melt(fdr,
                        id.vars="response",
                        variable.name = "cov")
            rownames(fdr) <- paste0(fdr$response, "_", fdr$cov)
            fdr <- fdr[rownames(data),]
            tmp <- rep("not_significant", nrow(pval))
            tmp[pval$value<=0.05] <- "significant"
            names(tmp) <- rownames(pval)
            data$pval <- pval[rownames(data), "value"]
            data$fdr <- fdr[rownames(data), "value"]
            data$significativity <- tmp[rownames(data)]
            data$sign <- rep("negative", nrow(data))
            data$sign[data$value>0]="positive"
            data$cov <- as.character(data$cov)
            data$cov <- gsub("_cov", "", data$cov)
            data$significativity[data$significativity=="significant" &
                                   data$cov=="cnv"] <- "significant_cnv"
            data$significativity[data$significativity=="significant" &
                                   data$cov=="met"] <- "significant_met"
            data$cnv_met <- rep(NA, nrow(data))
            data$cnv_met[data$cov=="cnv"] <- "cnv"
            data$cnv_met[data$cov=="met"] <- "met"
            data$cov[data$cov%in%c("cov", "cnv", "met")] <-
              data$response[data$cov%in%c("cov", "cnv", "met")]

            if(is.null(genes_info)){
              genes_info <- .download_gene_info(c(data$cov, data$response),
                                               filters=filters,
                                               species = species,
                                               ...)
              }
            tmp2 <- data$cov
            tmp <- intersect(tmp2, rownames(genes_info))
            tmp <- genes_info[tmp,]
            data$chr_cov <- tmp[as.character(tmp2),"chromosome_name"]
            data$cytoband_cov <- tmp[as.character(tmp2), "band"]
            data$start_cov <- tmp[as.character(tmp2), "start_position"]
            data$end_cov <- tmp[as.character(tmp2), "end_position"]
            data$ensg_cov <- tmp[as.character(tmp2), "ensembl_gene_id"]
            data$entrez_cov <- tmp[as.character(tmp2), "entrezgene_id"]
            data$hgnc_cov <- tmp[as.character(tmp2), "hgnc_symbol"]

            tmp2 <- data$response
            tmp <- intersect(tmp2, rownames(genes_info))
            tmp <- genes_info[tmp,]
            data$chr_response <- tmp[as.character(tmp2),"chromosome_name"]
            data$cytoband_response <- tmp[as.character(tmp2),"band"]
            data$start_response <- tmp[as.character(tmp2),"start_position"]
            data$end_response <- tmp[as.character(tmp2),"end_position"]
            data$ensg_response <- tmp[as.character(tmp2),"ensembl_gene_id"]
            data$entrez_response <- tmp[as.character(tmp2),"entrezgene_id"]
            data$hgnc_response <- tmp[as.character(tmp2), "hgnc_symbol"]

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

###################################################
setMethod("extract_model_res", "MultiClass",
          function(model_results,
                   outliers=T,
                   species="hsa",
                   filters=c("hgnc_symbol",
                             "ensembl_gene_id",
                             "entrezgene_id"),
                   genes_info=NULL,
                   ...){

      deg <- model_results$deg
      model_results <- lapply(model_results, function(x) x)
      data <- extract_model_res(model_results,
                                outliers=outliers,
                                species=species,
                                genes_info=genes_info,
                                filters=filters,
                                ...)
      data$deg <- rep(F, nrow(data))
      data$deg[data$response%in%deg] <- T
      return(data)
          })


####################################################
#' @importFrom plyr rbind.fill

setMethod("extract_model_res", "MultiOmics",
          function(model_results,
                   outliers=T,
                   species="hsa",
                   filters=c("hgnc_symbol",
                             "ensembl_gene_id",
                             "entrezgene_id"),
                   genes_info=NULL,
                   ...){

            tmp <- model_results
            if(is.null(genes_info)){
                  if(names(model_results[[1]])[
                    length(model_results[[1]])]!="data"){
                    tmp <- unlist(tmp, recursive = F)
                    tmp2 <- names(tmp)[-grep("\\.deg$", names(tmp))]
                    tmp <- lapply(tmp2, function(x) tmp[[x]])
                    names(tmp) <- tmp2
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
                                                species = species,
                                                ...)
            }

            tmp <- lapply(model_results, function(x)
              extract_model_res(x,
                                outliers=outliers,
                                species=species,
                                genes_info=genes_info,
                                filters=filters))
            data <- lapply(names(tmp), function(x) cbind(tmp[[x]], omics=x))
            names(data) <- names(tmp)
            data <- rbind.fill(data)
            if("deg"%in%colnames(data)){
              deg <- unique(data$response[data$deg])
              data$deg[data$response%in%deg] <- T
            }
            return(data)
          }
)



##################################
search_gene <- function(genes,
                        model_res=NULL,
                        data_frame=NULL){
  if(sum(is.null(model_res), is.null(data_frame))>1){
    stop(str_wrap("One of model_res amd data_frame should be provided"))
  }
  genes <- genes[!is.na(genes)]
  genes <- as.character(genes)
  genes <- unique(genes)
  if(is.null(data_frame)) data_frame <- extract_model_res(model_res)
  ans <- rbind(data_frame[data_frame$ensg_cov%in%genes,],
               data_frame[data_frame$entrez_cov%in%genes,],
               data_frame[data_frame$hgnc_cov%in%genes,],
               data_frame[data_frame$ensg_response%in%genes,],
               data_frame[data_frame$entrez_response%in%genes,],
               data_frame[data_frame$hgnc_response%in%genes,])
  return(ans)

}

#########################################
#' @import edgeR limma

.find_deg <- function(eexpression,
                      class,
                      RNAseq=T,
                      norm_method="TMM",
                      normalize=T){

  tmp <- apply(eexpression, 1, function(x) sum(is.na(x)))
  eexpression <- eexpression[tmp==0,]
  design_mat <- model.matrix(~class, data = as.data.frame(t(eexpression)))
  if(RNAseq){
    y <- DGEList(counts = eexpression)
    y <- calcNormFactors(y, method = norm_method)
    if(!normalize) y$offset <- matrix(1, ncol(eexpression),
                                      nrow(eexpression))
    y <- estimateGLMCommonDisp(y, design_mat)
    y <- estimateGLMTagwiseDisp(y, design_mat)
    fit <- glmFit(y, design_mat)
    lrt <- glmLRT(fit)
    top <- topTags(lrt, n = nrow(eexpression))$table
  }else{
    if(normalize) eexpression <- t(.data_norm(t(eexpression), norm_method))
    fit <- lmFit(eexpression, design_mat)
    fit <- eBayes(fit)
    top <- topTable(fit, number = nrow(eexpression))
    colnames(top) <- gsub("adj.P.Val", "FDR", colnames(top))
  }
    return(top)

}

############################################

#' @importFrom BiocGenerics lapply
setMethod("lapply", "MultiClass",
          function(X, FUN){
            ans <- X
            attributes(ans)$class <- "list"
            tmp <- names(ans)[names(ans)!="deg"]
            ans <- lapply(tmp, function(y) ans[[y]])
            names(ans) <- tmp
            return(lapply(ans, FUN))
          })


#'
#' @importFrom stats p.adjust

fdr <- function(pval_mat){

  fdr_data <- apply(pval_mat, 2, function(x){
    ans <- p.adjust(x[!is.na(x)], method = "fdr", )
    x[!is.na(x)] <- ans
    return(x)
  })
  fdr_data <- as.data.frame(fdr_data)
  return(fdr_data)
}

##############################################

.shiny_preprocess <- function(data){

  data_table <- data
  data_table <- filter(data_table, cov != '(Intercept)')   # elimino intercetta
  rownames(data_table) <- 1:nrow(data_table)
  data_table <- data_table[, !(colnames(data_table) %in% c('significativity', 'sign'))]
  colnames(data)[colnames(data) == "response"] <- "gene"
  colnames(data)[colnames(data) == "start_response"] <- "start"
  colnames(data)[colnames(data) == "end_response"] <-"end"
  colnames(data)[colnames(data) == "chr_response"] <-"chr"
  data$direction_cov <- ifelse(data$cov_value < 0, 'negative', 'positive')
  data$direction_coef <- ifelse(data$coef < 0, 'negative', 'positive')
  data$cov_value <- abs(data$cov_value)
  data$coef <- abs(data$coef)
  ans <- list(data=data, data_table=data_table)
  return(ans)
}

#############################################################################

.circos_preprocess <- function(data){

  library(GenomicRanges)
  dataframes <- lapply(unique(data$omics), function(x) {
    single_omic_df <- data[data$omics==x,]
    return(single_omic_df)
  })
  names(dataframes) <- paste0("df_", unique(data$omics))
  if("df_gene_genomic_res"%in%names(dataframes)){
    dataframes$df_cnv <- filter(dataframes$df_gene_genomic_res,
                                cnv_met == 'cnv')
    dataframes$df_met <- filter(dataframes$df_gene_genomic_res,
                                cnv_met == 'met')
    tmp <- lapply(which(names(dataframes)!="df_gene_genomic_res"), function(x){
      ans <- dataframes[[x]]
      ans <- ans[ans$cov!="(Intercept)",]
      return(ans)
    })
    names(tmp) <- names(dataframes)[
      which(names(dataframes)!="df_gene_genomic_res")]
    dataframes <- tmp
  }
  gr <- lapply(dataframes, function(x){
    ans <- makeGRangesFromDataFrame(x,
                                    seqnames.field = 'chr',
                                    start.field = 'start',
                                    end.field = 'end',
                                    keep.extra.columns = TRUE,
                                    na.rm = TRUE)
    return(ans)
  })
  return(gr)
}

# funzione per creare una singola track. questa viene usata nella funzione .create_tracks che in base ai dati disponibili crea tutte le track possibili che poi andranno unite (composed_view) e usate in arrange_views
.create_single_track <- function(data,
                                 dataValue,
                                 x_axis,
                                 xe_axis,
                                 y_axis,
                                 colorField,
                                 colorDomain,
                                 colorRange,
                                 tooltipField1,
                                 tooltipTitle,
                                 tooltipAlt1,
                                 tooltipField2,
                                 tooltipAlt2,
                                 tooltipField3,
                                 tooltipAlt3,
                                 tooltipField4,
                                 tooltipAlt4,
                                 legend) {
  return(
    add_single_track(
      data = track_data_gr(data, chromosomeField = 'seqnames', genomicFields = c('start','end'), value = dataValue),
      mark = 'bar',
      x = visual_channel_x(field = 'start', type = 'genomic', axis = x_axis),
      xe = visual_channel_x(field = 'end', type = 'genomic', axis = xe_axis),
      y = visual_channel_y(field = yValue, type = 'quantitative', axis = y_axis),
      color = visual_channel_color(field = colorField, type = 'nominal', domain = colorDomain, range = colorRange),
      tooltip = visual_channel_tooltips(
        visual_channel_tooltip(field = "start", type = "genomic", alt = 'Start Position:'),
        visual_channel_tooltip(field = "end", type = "genomic", alt = "End Position:"),
        visual_channel_tooltip(field = tooltipField1, title = tooltipTitle, type = "quantitative", alt = paste(tooltipTitle, "Value:"), format = "0.2"),
        visual_channel_tooltip(field = tooltipField2, type = 'nominal', alt = tooltipAlt2),
        visual_channel_tooltip(field = tooltipField3, type = 'nominal', alt = tooltipAlt3),
        visual_channel_tooltip(field = tooltipField4, type = 'nominal', alt = tooltipAlt4)
      ),
      size = list(value = 1),
      legend = legend
    )
  )
}

#######################################################################
########################################################################

.create_tracks <- function(data_table, gr){

tracks <- list()

track_cyto <- add_single_track(
  id = "track2",
  data = track_data(
    url = "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
    type = "csv",
    chromosomeField = "Chromosome",
    genomicFields = c("chromStart",
                      "chromEnd")
  ),
  mark = "rect",
  x = visual_channel_x(field = "chromStart",
                       type = "genomic"),
  xe = visual_channel_x(field = "chromEnd",
                        type = "genomic"),
  color = visual_channel_color(
    field = "Stain",
    type = "nominal",
    domain = c(
      "gneg",
      "gpos25",
      "gpos50",
      "gpos75",
      "gpos100",
      "gvar",
      "acen"           # acen: centromeric region (UCSC band files)
    ),
    range = c(
      "white",
      "#D9D9D9",
      "#979797",
      "#636363",
      "black",
      "#F0F0F0",
      "red"
    )
  ),
  stroke = visual_channel_stroke(
    value = "lightgray"
  ),
  strokeWidth = visual_channel_stroke_width(
    value = 0.5
  ),
)

  if ("gene_genomic_res" %in% unique(data_table$omics)){
    #####da mettere in un'altra funzione
    track_cnv <- create_single_track(data=cnv_gr,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="cnv",
                                     tooltipAlt1="CNV Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_met <- create_single_track(data=met_gr,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="met",
                                     tooltipAlt1="MET Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_expr <- create_single_track(data=cnv_gr,
                                     dataValue='response_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="response_value",
                                     tooltipTitle="expr",
                                     tooltipAlt1="Expression Value (log2):",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)
   tracks <- c(tracks, track_cnv=track_cnv, track_met=track_met, track_expr=track_expr, track_cyto=track_cyto)
    }


  if ("cnv_gene_res" %in% unique(data_table$omics)){
    track_cnv <- create_single_track(data=cnv_gr,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="cnv",
                                     tooltipAlt1="CNV Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_expr <- create_single_track(data=cnv_gr,
                                      dataValue='response_value',
                                      x_axis="none",
                                      xe_axis="none",
                                      y_axis="none",
                                      colorField="direction_cov",
                                      colorDomain=c("positive","negative"),
                                      colorRange=c("red","blue"),
                                      tooltipField1="response_value",
                                      tooltipTitle="expr",
                                      tooltipAlt1="Expression Value (log2):",
                                      tooltipField2="gene",
                                      tooltipAlt2="Gene Name:",
                                      tooltipField3="class",
                                      tooltipAlt3="Class:",
                                      tooltipField4="cnv_met",
                                      tooltipAlt4="Integration Type:",
                                      legend=FALSE)

    tracks <- c(tracks,track_cnv=track_cnv, track_expr=track_expr, track_cyto=track_cyto)
    }

  if("met_gene_res"%in%unique(data_table$omics)){

    track_met <- create_single_track(data=met_gr,
                                     dataValue='cov_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="direction_cov",
                                     colorDomain=c("positive","negative"),
                                     colorRange=c("red","blue"),
                                     tooltipField1="cov_value",
                                     tooltipTitle="met",
                                     tooltipAlt1="MET Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     legend=FALSE)

    track_expr <- create_single_track(data=cnv_gr,
                                      dataValue='response_value',
                                      x_axis="none",
                                      xe_axis="none",
                                      y_axis="none",
                                      colorField="direction_cov",
                                      colorDomain=c("positive","negative"),
                                      colorRange=c("red","blue"),
                                      tooltipField1="response_value",
                                      tooltipTitle="expr",
                                      tooltipAlt1="Expression Value (log2):",
                                      tooltipField2="gene",
                                      tooltipAlt2="Gene Name:",
                                      tooltipField3="class",
                                      tooltipAlt3="Class:",
                                      tooltipField4="cnv_met",
                                      tooltipAlt4="Integration Type:",
                                      legend=FALSE)

    tracks <- c(tracks, track_met=track_met, track_expr=track_expr, track_cyto=track_cyto)
    }

return(tracks)
}
