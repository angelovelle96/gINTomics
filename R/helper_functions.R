#' Data check
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
    check_sd <- apply(response_var, 2, function(x) sd(x, na.rm = TRUE))
    check_sd[is.na(check_sd)] <- 0
    if(sum(check_sd==0)>0) {
      message(str_wrap("removing response variables with
                                zero standard deviation"))
      response_var <- response_var[, check_sd>0]
    }

    if(is.numeric(covariates[,1])){
      check_sd <- apply(covariates, 2, function(x) sd(x, na.rm = TRUE))
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
    if(length(tmp)==0){
      warning(str_wrap("No genes left in common between
                        response_var and interactions"))
      return(NULL)
      }
    response_var <- as.matrix(response_var[,tmp])
    colnames(response_var) <- tmp
    interactions <- lapply(tmp, function(x) interactions[[x]])
    names(interactions) <- tmp

    rresult <- list(response_var=response_var,
                    covariates=covariates,
                    interactions=interactions)
    return(rresult)
}

#' Covariates check
#' @importFrom stats relevel

.covariates_check <- function(response_var,
                             covariates,
                             interactions,
                             steady_covariates=NULL,
                             linear=FALSE,
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
  tmp <- intersect(colnames(covariates), unlist(interactions))
  if(length(tmp)==0) return(NULL)
  covariates <- as.data.frame(covariates[,tmp])
  colnames(covariates) <- tmp
  if(linear==TRUE){
    colnames(response_var) <- gsub(bad_str, "_", colnames(response_var))
    names(interactions) <- gsub(bad_str, "_", names(interactions))
  }
  if(!is.null(reference)){
    covariates[vapply(covariates, is.factor, FUN.VALUE = TRUE)] <- lapply(
      covariates[vapply(covariates, is.factor, FUN.VALUE = TRUE)],
      function(x) relevel(x, ref = reference))
  }
  return(list(covariates=covariates,
              response_var=response_var,
              interactions=interactions,
              original_id=original_id,
              steady_covariates=steady_covariates))
}

#' Generating formula for models
#' @importFrom stats formula

.generate_formula <- function(interactions,
                             linear=FALSE){

  fformula <- lapply(seq_along(interactions), function(x){
    cov <- interactions[[x]]
    ans <- formula(paste0("~", paste0(cov, collapse = "+")))
    if(linear==TRUE){
      ans <- formula(paste0(names(interactions)[x],
                            "~", as.character(ans)[2]))
    }
    return(ans)
  })
  names(fformula) <- names(interactions)
  return(fformula)
}

#' Building result matrices
#' @importFrom plyr rbind.fill

.building_result_matrices <- function(model_results,
                                     type,
                                     single_cov=FALSE){
  ####coef extraction
  if(type=="edgeR"){
    tmp <- lapply(model_results, function(x) as.data.frame( x["coef",],
                                                            check.names = FALSE))
  }
  if(type=="lm"){
    tmp <- lapply(model_results, function(x){
      tmp <- data.frame(t(x[["coefficients"]][, "Estimate"]),check.names = FALSE)
      colnames(tmp) <- rownames(x[["coefficients"]])
      return(tmp)
    })
  }
  if(single_cov==TRUE){
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
  for(i in seq_len(ncol(coef_matrix))) {
    tmp <- coef_matrix[, i]
    tmp[vapply(tmp, is.null, FUN.VALUE = TRUE)] <- NA
    tmp <- unlist(tmp)
    coef_matrix[, i] <- tmp
  }

  ######pvalue extraction
  if(type=="edgeR"){
    tmp <- lapply(model_results, function(x) as.data.frame( x["PValue",],
                                                            check.names = FALSE))
  }
  if(type=="lm"){
    tmp <- lapply(model_results, function(x){
      tmp <- data.frame(t(x[["coefficients"]][, "Pr(>|t|)"]),check.names = FALSE)
      colnames(tmp) <- rownames(x[["coefficients"]])
      return(tmp)
    })
  }
  if(single_cov==TRUE){
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
  for(i in seq_along(pval_matrix)) {
    tmp <- pval_matrix[, i]
    tmp[vapply(tmp, is.null, FUN.VALUE = TRUE)] <- NA
    tmp <- unlist(tmp)
    pval_matrix[, i] <- tmp
  }
  coef_pval_mat <- list()
  coef_pval_mat[["coef"]] <- coef_matrix
  coef_pval_mat[["pval"]] <- pval_matrix
  return(coef_pval_mat)
}


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
#' @param ... Additional arguments to be passed to the function
#' @return A MultiAssayExperiment object containing the provided assays.
#' @examples
#' # Example usage:
#' create_multiassay(methylation, cnv_data, gene_exp, miRNA_exp, miRNA_cnv_data)
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


#' Data normalization
#' @importFrom limma normalizeBetweenArrays
#' @importFrom edgeR DGEList calcNormFactors cpm
    .data_norm <- function(data,
                          method="TMM",
                          RNAseq=TRUE){
      if(RNAseq){
          data <- DGEList(t(data))
          data <- calcNormFactors(data, method=method)
          data <- t(cpm(data, log = FALSE))
      }else{
          data <- normalizeBetweenArrays(data)
          }
      return(data)
    }

#' ID conversion
#' @importFrom stringi stri_replace_all_regex

    .id_conversion <- function(dictionary,
                              results){

      dictionary$tranformed <- paste0("^", dictionary$tranformed, "$")

      rownames(results$residuals) <- stri_replace_all_regex(
        rownames(results$residuals),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      colnames(results$coef_data) <- stri_replace_all_regex(
        colnames(results$coef_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      rownames(results$coef_data) <- stri_replace_all_regex(
        rownames(results$coef_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      colnames(results$pval_data) <- stri_replace_all_regex(
        colnames(results$pval_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      rownames(results$pval_data) <- stri_replace_all_regex(
        rownames(results$pval_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      colnames(results$fdr_data) <- stri_replace_all_regex(
        colnames(results$fdr_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      rownames(results$fdr_data) <- stri_replace_all_regex(
        rownames(results$fdr_data),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      colnames(results$data$response_var) <- stri_replace_all_regex(
        colnames(results$data$response_var),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)

      colnames(results$data$covariates) <- stri_replace_all_regex(
        colnames(results$data$covariates),
        pattern = dictionary$tranformed,
        replacement = dictionary$original,
        vectorize_all = FALSE)
       return(results)
    }


#' RandomForest selection
#' @importFrom randomForest randomForest importance

.rf_selection <- function(data,
                          formula){

  if(length(attr(terms(formula), "term.labels"))>
     as.integer(nrow(data)*0.4)){
    ans <- randomForest(formula, data)
    ans <- importance(ans)
    ans <- ans[order(ans[,1], decreasing = TRUE),]
    tmp <- as.integer(nrow(data)*0.4)
    ans <- ans[seq_len(tmp)]
    tmp <- paste(names(ans), collapse = "+")
    tmp <- paste(as.character(formula[2]), "~", tmp)
    ans <- as.formula(tmp)
  }else{
      ans <- formula
  }

  ans <- as.formula(paste(as.character(ans[1]), as.character(ans[3])))
  return(ans)
}


#' @rdname extract_model_res
setMethod("extract_model_res", "list",
          function(model_results,
                   outliers=TRUE,
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
            data$sign[data$value>0] <- "positive"
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
            tmp2 <- data.frame(gene = as.character(data$cov))
            tmp <- intersect(tmp2$gene, rownames(genes_info))
            tmp <- genes_info[tmp,]
            tmp$gene <- rownames(tmp)
            tmp <- left_join(tmp2, tmp, by = "gene")
            data$chr_cov <- tmp[,"chromosome_name"]
            data$cytoband_cov <- tmp[, "band"]
            data$start_cov <- as.numeric(tmp[, "start_position"])
            data$end_cov <- as.numeric(tmp[, "end_position"])
            data$ensg_cov <- tmp[, "ensembl_gene_id"]
            data$entrez_cov <- tmp[, "entrezgene_id"]
            data$hgnc_cov <- tmp[, "hgnc_symbol"]

            tmp2 <- data.frame(gene = as.character(data$response))
            tmp <- intersect(tmp2$gene, rownames(genes_info))
            tmp <- genes_info[tmp,]
            tmp$gene <- rownames(tmp)
            tmp <- left_join(tmp2, tmp, by = "gene")
            data$chr_response <- tmp[,"chromosome_name"]
            data$cytoband_response <- tmp[,"band"]
            data$start_response <- as.numeric(tmp[,"start_position"])
            data$end_response <- as.numeric(tmp[,"end_position"])
            data$ensg_response <- tmp[,"ensembl_gene_id"]
            data$entrez_response <- tmp[,"entrezgene_id"]
            data$hgnc_response <- tmp[, "hgnc_symbol"]

            colnames(data)[colnames(data)=="value"] <- "coef"
            tmp <-  apply(model_results$data$response_var, 2, mean)
            data$response_value <- tmp[data$response]
            tmp <-  apply(model_results$data$covariates, 2, mean)
            names(tmp) <- gsub("_cov", "", names(tmp))
            tmp2 <- c(length(grep("_cnv$", names(tmp)))>0,
                      length(grep("_met$", names(tmp)))>0)
            if(sum(tmp2)==2){
              data$cov_value <- tmp[rownames(data)]
              }else{data$cov_value <- tmp[data$cov]}
            mmin <- quantile(data$coef, 0.25) - 1.5*IQR(data$coef)
            mmax <- quantile(data$coef, 0.75) + 1.5*IQR(data$coef)

            if(outliers==FALSE){
              data <- data[data$coef>mmin,]
              data <- data[data$coef<mmax,]
            }

            return(data)

          }
)

#' @rdname extract_model_res
setMethod("extract_model_res", "MultiClass",
          function(model_results,
                   outliers=TRUE,
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
      data$deg <- rep(FALSE, nrow(data))
      data$deg[data$response%in%deg] <- TRUE
      return(data)
          })


#' @rdname extract_model_res
setMethod("extract_model_res", "MultiOmics",
          function(model_results,
                   outliers=TRUE,
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
                    tmp <- unlist(tmp, recursive = FALSE)
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
              data$deg[data$response%in%deg] <- TRUE
            }
            return(data)
          }
)


#' Search genes
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

#' Finding DEGs
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#' estimateGLMTagwiseDisp glmFit glmLRT topTags
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix

.find_deg <- function(eexpression,
                      class,
                      RNAseq=TRUE,
                      norm_method="TMM",
                      normalize=TRUE){

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

#' Setting method for lapply
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


#' FDR calculation
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


#' Shiny data preprocessing
#' @importFrom dplyr filter
.shiny_preprocess <- function(data){

  data_table <- data
  data_table <- filter(data_table, cov != '(Intercept)')
  rownames(data_table) <- seq_len(nrow(data_table))
  data_table$cnv_met[data_table$omics=="gene_cnv_res"] <- "cnv"
  data_table$cnv_met[data_table$omics=="gene_met_res"] <- "met"
  data_table <- data_table[, !(colnames(data_table) %in% c('significativity', 'sign'))]
  colnames(data)[colnames(data) == "response"] <- "gene"
  colnames(data)[colnames(data) == "start_response"] <- "start"
  colnames(data)[colnames(data) == "end_response"] <-"end"
  colnames(data)[colnames(data) == "chr_response"] <-"chr"
  ans <- list(data=data, data_table=data_table)
  return(ans)
}

#' Change integration names
.change_int_names <- function(nnames){

  tmp <- gsub("gene_genomic_res", "Gene CNV-Met", nnames)
  tmp <- gsub("gene_cnv_res", "Gene CNV", tmp)
  tmp <- gsub("gene_met_res", " Gene Met", tmp)
  tmp <- gsub("mirna_cnv_res", "miRNA CNV", tmp)
  tmp <- gsub("tf_res", "Gene's TFs", tmp)
  tmp <- gsub("tf_mirna_res", "miRNA's TFs", tmp)
  tmp <- gsub("mirna_target_res", "miRNA's targets", tmp)
  names(nnames) <- tmp
  return(nnames)
}


#' @importClassesFrom edgeR DGEGLM
#' @importFrom stats poisson
#' @importFrom MASS negative.binomial
#' @importFrom stats residuals
#' @exportS3Method gINTomics::residuals
residuals.DGEGLM <- function(object, type=c("deviance", "pearson"), ...) {
  y <- as.matrix(object$counts)
  mu <- as.matrix(object$fitted.values)
  theta <- 1/object$dispersion
  if(is.null(object$weights)) {
    wts <- rep(1, ncol(object$counts))
  } else {
    wts <- as.matrix(object$weights)
  }
  type <- match.arg(type)
  ymut <- cbind(y, mu, theta)

  res <- t(apply(ymut, 1, function(x) {
    yy <- as.vector(x[seq_len(ncol(y))])
    mm <- as.vector(x[seq((ncol(y)+1),(ncol(y)+ncol(mu)))])
    t <- x[length(x)]
    if(type=="deviance") {
      if(t==Inf) {
        d.res <- sqrt(pmax((poisson()$dev.resids)(yy, mm, wts), 0))
      } else {
        d.res <- sqrt(pmax((negative.binomial(theta=t)$dev.resids)(yy, pmax(mm, 1e-8), wts), 0))
      }
      return(ifelse(yy > mm, d.res, -d.res))
    } else if(type=="pearson") {
      if(t==Inf) {
        return((yy - mm) * sqrt(wts) / pmax(sqrt(poisson()$variance(mm)), 1))
      } else {
        return((yy - mm) * sqrt(wts) / pmax(sqrt(negative.binomial(theta=t)$variance(mm)), 1))
      }
    }
  }))
  return(res)
}


