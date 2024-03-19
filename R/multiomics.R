
#' Complete Multi-Omics integration
#' @description
#' This function will perform a complete Multi-Omics integration on a
#' MultiAssayExperiment
#' @param data A MultiAssayExperiment. It can be generated exploiting the
#' **generate_multiassay** function.
#' @param interactions_met **interactions** as for **run_met_integration**
#' @param interactions_miRNA_target miRNA-target interactions as requested by
#'  **run_tf_integration**
#' @param interactions_tf TF-target interactions as requested by
#'  **run_tf_integration**
#' @param interactions_tf_miRNA TF-target interactions as
#' requested by **run_tf_integration**
#' @param RNAseq logical. Are gene expression data obtained from RNA
#' sequencing ? Default is set to TRUE
#' @param miRNAseq logical. Are miRNA expression data obtained from miRNA
#' sequencing ? Default is set to TRUE
#' @param normalize_miRNA_expr logical.Should miRNA expression data be
#' normalized ? Default is set to TRUE
#' @param normalize_gene_expr logical.Should gene expression data be
#' normalized ? Default is set to TRUE
#' @param norm_method_gene_expr Normalization method to be used for gene
#' expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
#' "none".
#' @param norm_method_miRNA_expr Normalization method to be used for miRNA
#' expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
#' "none".
#' @param class Character vector specifying the classes for differential
#' expression analysis.
#' @param BPPARAM A BiocParallelParam object specifying the parallel backend to
#' be used.
#' @import MultiAssayExperiment SummarizedExperiment
#' @importFrom methods new
#' @importFrom plyr rbind.fill
#' @importFrom methods is
#' @import BiocParallel
#' @return A \linkS4class{MultiOmics} object containing the results of all the
#' possible integration models
#' @examples
#' # Example usage_multiomics:
#' data("ov_test_tcga_omics")
#' multiomics_integration <-run_multiomics(data = mmultiassay_ov)
#' @export
run_multiomics <- function(data,
                           interactions_met=NULL,
                           interactions_miRNA_target=NULL,
                           interactions_tf=NULL,
                           interactions_tf_miRNA=NULL,
                           RNAseq=TRUE,
                           miRNAseq=TRUE,
                           normalize_miRNA_expr=TRUE,
                           normalize_gene_expr=TRUE,
                           norm_method_gene_expr="TMM",
                           norm_method_miRNA_expr="TMM",
                           class=NULL,
                           BPPARAM=SerialParam()){



    if(!is(data, "MultiAssayExperiment")){
      stop(str_wrap("You need to provide a MultiassayExperiment as input
                    data. The names of the experiments included in the
                    MultiAssay should follow the rule provided in the
                    help section. You can use the function create_multiassay
                    to create proper input data"))
    }


    normalize_gene_expr2 <- normalize_gene_expr
    normalize_miRNA_expr2 <- normalize_miRNA_expr
    if(!is.null(data@ExperimentList$gene_exp))
      data <- c(data, gene_exp_original=data@ExperimentList$gene_exp)
    if(!is.null(data@ExperimentList$miRNA_exp))
      data <- c(data, miRNA_exp_original=data@ExperimentList$miRNA_exp)

    deg_gene <- TRUE
    deg_mirna <- TRUE
    gene_genomic_res <- NULL
    geno <- FALSE
    if(!is.null(data@ExperimentList$cnv_data) &
       !is.null(data@ExperimentList$gene_exp)&
       !is.null(data@ExperimentList$methylation)){
      message("--------------Running gene genomic integration--------------")
      gene_genomic_res <- run_genomic_integration(
        expression = t(assay(data, i = "gene_exp")),
        cnv_data = t(assay(data, i = "cnv_data")),
        methylation = t(assay(data, i = "methylation")),
        sequencing_data = RNAseq,
        normalize=normalize_gene_expr,
        norm_method=norm_method_gene_expr,
        class=class,
        run_deg=deg_gene,
        BPPARAM=BPPARAM)
      if(!is.null(gene_genomic_res)){
          if(!is.null(class)){
            tmp <- lapply(gene_genomic_res, function(x)
              as.data.frame(t(x$residuals)))
            tmp2 <- rbind.fill(tmp)
            rownames(tmp2) <- unlist(lapply(tmp, rownames))
            if(!is.null(tmp2)) rresiduals <- t(tmp2)
          }else{
            rresiduals <- gene_genomic_res$residuals
          }
          data@ExperimentList$gene_exp <- data@ExperimentList$gene_exp[
            rownames(rresiduals),]
          assay(data@ExperimentList$gene_exp) <- as.matrix(rresiduals[,
            colnames(assay(data@ExperimentList$gene_exp))])
          RNAseq <- FALSE
          normalize_gene_expr2 <- FALSE
          geno <- TRUE
          deg_gene <- FALSE
        }
      }



    gene_cnv_res <- NULL
    if(!is.null(data@ExperimentList$cnv_data) &
       !is.null(data@ExperimentList$gene_exp) &
       geno==FALSE){
      message("----------------Running gene CNV integration----------------")
      gene_cnv_res <- run_cnv_integration(
        expression = t(assay(data, i = "gene_exp")),
        cnv_data = t(assay(data, i = "cnv_data")),
        sequencing_data = RNAseq,
        normalize=normalize_gene_expr,
        norm_method=norm_method_gene_expr,
        class=class,
        run_deg=deg_gene,
        BPPARAM=BPPARAM)
      if(!is.null(gene_cnv_res)){
          if(!is.null(class)){
          tmp <- lapply(gene_cnv_res, function(x)
            as.data.frame(t(x$residuals)))
          tmp2 <- rbind.fill(tmp)
          rownames(tmp2) <- unlist(lapply(tmp, rownames))
          if(!is.null(tmp2)) rresiduals <- t(tmp2)
        }else{
          rresiduals <- gene_cnv_res$residuals
        }
          data@ExperimentList$gene_exp <- data@ExperimentList$gene_exp[
            rownames(rresiduals),]
          assay(data@ExperimentList$gene_exp) <- as.matrix(rresiduals[,
            colnames(assay(data@ExperimentList$gene_exp))])
        RNAseq <- FALSE
        normalize_gene_expr2 <- FALSE
        deg_gene <- FALSE
        }
      }

    gene_met_res <- NULL
    if(!is.null(data@ExperimentList$methylation) &
       !is.null(data@ExperimentList$gene_exp) &
       geno==FALSE){
      message("------------Running gene methylation integration------------")
      gene_met_res <- run_met_integration(
        expression = t(assay(data, i = "gene_exp")),
        methylation = t(assay(data, i = "methylation")),
        sequencing_data = RNAseq,
        class=class,
        run_deg=deg_gene,
        normalize = normalize_gene_expr2,
        norm_method=norm_method_gene_expr,
        BPPARAM=BPPARAM)
      if(!is.null(gene_met_res)) deg_gene <- FALSE
    }

    mirna_cnv_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_cnv_data) &
       !is.null(data@ExperimentList$miRNA_exp)){
      message("---------------Running miRNA CNV integration----------------")
      mirna_cnv_res <- run_cnv_integration(
        expression = t(assay(data, i = "miRNA_exp")),
        cnv_data = t(assay(data, i = "miRNA_cnv_data")),
        sequencing_data = miRNAseq,
        normalize=normalize_miRNA_expr,
        norm_method=norm_method_miRNA_expr,
        class=class,
        run_deg=deg_mirna,
        BPPARAM=BPPARAM)
      if(!is.null(mirna_cnv_res)){
          if(!is.null(class)){
            tmp <- lapply(mirna_cnv_res, function(x)
              as.data.frame(t(x$residuals)))
            tmp2 <- rbind.fill(tmp)
            rownames(tmp2) <- unlist(lapply(tmp, rownames))
            if(!is.null(tmp2)) rresiduals <- t(tmp2)
          }else{
            rresiduals <- mirna_cnv_res$residuals
          }
          data@ExperimentList$miRNA_exp <- data@ExperimentList$miRNA_exp[
            rownames(rresiduals),]
          assay(data@ExperimentList$miRNA_exp) <- as.matrix(rresiduals[,
            colnames(assay(data@ExperimentList$miRNA_exp))])
          miRNAseq <- FALSE
          normalize_miRNA_expr2 <- FALSE
          deg_mirna <- FALSE
        }
      }

    tf_res <- NULL
    if(!is.null(data@ExperimentList$gene_exp)){
      message("-------------------Running TF integration-------------------")
      tf_res <- run_tf_integration(
        expression = t(assay(data,i = "gene_exp")),
        tf_expression = t(assay(data,i = "gene_exp_original")),
        interactions = interactions_tf,
        sequencing_data = RNAseq,
        normalize=normalize_gene_expr2,
        norm_method=norm_method_gene_expr,
        normalize_cov = normalize_gene_expr,
        norm_method_cov = norm_method_gene_expr,
        class=class,
        run_deg=deg_gene,
        type="tf",
        BPPARAM=BPPARAM)
      if(!is.null(tf_res)) deg_gene <- FALSE
    }


    tf_mirna_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_exp) &
       !is.null(data@ExperimentList$gene_exp)){
      message("----------------Running TF miRNA integration----------------")
      tf_mirna_res <- run_tf_integration(
        expression = t(assay(data,i = "miRNA_exp")),
        tf_expression = t(assay(data,i = "gene_exp_original")),
        interactions = interactions_tf_miRNA,
        sequencing_data = miRNAseq,
        normalize=normalize_miRNA_expr2,
        norm_method=norm_method_miRNA_expr,
        normalize_cov = normalize_gene_expr,
        norm_method_cov = norm_method_gene_expr,
        class=class,
        run_deg=deg_mirna,
        type="tf_miRNA",
        BPPARAM=BPPARAM)
      if(!is.null(tf_mirna_res)) deg_mirna <- FALSE
    }


    mirna_target_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_exp) &
       !is.null(data@ExperimentList$gene_exp)){
      message("--------------Running miRNA target integration--------------")
      mirna_target_res <- run_tf_integration(
        expression = t(assay(data,i = "gene_exp")),
        tf_expression = t(assay(data,i = "miRNA_exp_original")),
        interactions = interactions_miRNA_target,
        sequencing_data = RNAseq,
        normalize=normalize_gene_expr2,
        norm_method=norm_method_gene_expr,
        normalize_cov = normalize_miRNA_expr,
        norm_method_cov = norm_method_miRNA_expr,
        class=class,
        run_deg=deg_gene,
        type="miRNA_target",
        BPPARAM=BPPARAM)
      if(!is.null(mirna_target_res)) deg_gene <- FALSE
    }

  ans <- new("MultiOmics", Filter(Negate(is.null),
                                  list(gene_cnv_res=gene_cnv_res,
                                       gene_genomic_res=gene_genomic_res,
                                       mirna_cnv_res=mirna_cnv_res,
                                       gene_met_res=gene_met_res,
                                       tf_res=tf_res,
                                       tf_mirna_res=tf_mirna_res,
                                       mirna_target_res=mirna_target_res)))
  return(ans)
}


# cnv integration
.def_cnv_integration <- function(expression,
                                cnv_data,
                                sequencing_data,
                                normalize,
                                norm_method,
                                BPPARAM,
                                ...){

  if(sequencing_data==TRUE){
    cnv_res <- .run_edgeR_integration(response_var = expression,
                                      covariates = cnv_data,
                                      normalize = normalize,
                                      norm_method = norm_method,
                                      BPPARAM = BPPARAM,
                                      ...)
    if(normalize & !is.null(cnv_res)){
      cnv_res$data$response_var <- .data_norm(cnv_res$data$response_var,
                                              method = norm_method)
    }
  }else{
    cnv_res <- .run_lm_integration(response_var = expression,
                                   covariates = cnv_data,
                                   normalize = normalize,
                                   norm_method = norm_method,
                                   BPPARAM = BPPARAM,
                                   ...)
  }
  return(cnv_res)
}


#' Integration of expression and Copy Number Variations
#' @description
#' This function will perform an integration of expression data and Copy Number
#' Variations data
#' @param expression Matrix or data.frame containing the expression values
#' for each model. Rows represent samples, while each column represents
#' the different response variables of the models.
#' @param cnv_data Matrix or data.frame containing the Copy Number variation
#' status for the models. Rows represent samples, while columns represent
#' the different covariates. If **interactions** are not provided, they will be
#' automatically generated and for each gene contained in **expression**
#' the model will look for the same gene in **cnv_data**
#' @param sequencing_data logical. Are expression data obtained from RNA
#' sequencing ? Default is set to TRUE
#' @param normalize logical.Should expression data be
#' normalized ? Default is set to TRUE
#' @param norm_method Normalization method to be used for
#' expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
#' "none".
#' @param class Character vector specifying the classes for differential
#' expression analysis.
#' @param run_deg Logical. Should differential expression analysis be performed?
#'  Default is set to TRUE.
#' @param BPPARAM A BiocParallelParam object specifying the parallel backend to
#' be used.
#' @param ... Additional arguments to be passed to internal functions.
#' @return A list or a \linkS4class{MultiClass} object if **class** is provided
#' containing the results of the CNV integration
#' @examples
#' # Example usage_multi:
#' data("ov_test_tcga_omics")
#' gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
#' gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
#' cnv_integration <- run_cnv_integration(expression=gene_exp_matrix,
#'  cnv_data=gene_cnv_matrix)
#' @export
#' @importFrom BiocParallel bpparam SerialParam
run_cnv_integration <- function(expression,
                                cnv_data,
                                sequencing_data=TRUE,
                                normalize=TRUE,
                                norm_method="TMM",
                                class=NULL,
                                run_deg=TRUE,
                                BPPARAM=SerialParam(),
                                ...){

  if(is.null(class)){
    cnv_res <- .def_cnv_integration(expression = expression,
                                   cnv_data = cnv_data,
                                   sequencing_data = sequencing_data,
                                   normalize = normalize,
                                   norm_method = norm_method,
                                   BPPARAM = BPPARAM,
                                   ...)
  }else{
    if(!is.character(class) |
       length(class)!=nrow(expression) |
       !identical(names(class),rownames(expression))){
        stop(str_wrap("class should be a named character vector, names should
                      match sample names contained in expression data"))
    }

    tmp <- unique(class)
    tmp <- lapply(tmp, function(x) names(class)[class==x])
    names(tmp) <- unique(class)
    cnv_res <- lapply(tmp, function(x){
      ans <- .def_cnv_integration(expression = expression[x,],
                                 cnv_data = cnv_data[x,],
                                 sequencing_data = sequencing_data,
                                 normalize = normalize,
                                 norm_method = norm_method,
                                 BPPARAM = BPPARAM,
                                 ...)
      return(ans)
    })
    deg <- NULL
    if(run_deg){
      deg <- .find_deg(eexpression = t(expression),
                       class = class,
                       normalize = normalize,
                       norm_method = norm_method)
      deg <- list(rownames(deg)[deg$FDR<=0.1])
      }
    cnv_res <- new("MultiClass", c(cnv_res, deg=deg))
  }
  return(cnv_res)
}


# defining met integration
.def_met_integration <- function(expression,
                                 methylation,
                                 sequencing_data,
                                 normalize,
                                 norm_method,
                                 BPPARAM,
                                 ...){

  if(sequencing_data==TRUE){
    met_res <- .run_edgeR_integration(response_var = expression,
                                      covariates = methylation,
                                      normalize = normalize,
                                      norm_method = norm_method,
                                      BPPARAM = BPPARAM,
                                      ...)
    if(normalize & !is.null(met_res)){
      met_res$data$response_var <- .data_norm(met_res$data$response_var,
                                              method = norm_method)
    }
  }else{
    met_res <- .run_lm_integration(response_var = expression,
                                   covariates = methylation,
                                   normalize = normalize,
                                   norm_method = norm_method,
                                   BPPARAM = BPPARAM,
                                   ...)
  }
  return(met_res)
}




#' Integration of expression and methylation
#' @description
#' This function will perform an integration of expression data and methylation
#'  data
#' @param expression Matrix or data.frame containing the expression values
#' for each model. Rows represent samples, while each column represents
#' the different response variables of the models.
#' @param methylation Matrix or data.frame containing the methylation
#' values for the models. Rows represent samples, while columns represent
#' the different covariates. If **interactions** are not provided, they will be
#' automatically generated and for each gene contained in **expression**
#' the model will look for the same gene in **methylation**
#' @param sequencing_data logical. Are expression data obtained from RNA
#' sequencing ? Default is set to TRUE
#' @param normalize logical.Should expression data be
#' normalized ? Default is set to TRUE
#' @param norm_method Normalization method to be used for
#' expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
#' "none".
#' @param class Character vector specifying the classes for differential
#' expression analysis.
#' @param run_deg Logical. Should differential expression analysis be performed?
#'  Default is set to TRUE.
#' @param BPPARAM A BiocParallelParam object specifying the parallel backend to
#'  be used.
#' @param ... Additional arguments to be passed to internal functions.
#' @return A list or a \linkS4class{MultiClass} object if **class** is provided
#' containing the results of the Methylation integration
#' @examples
#' # Example usage_multi:
#' data("ov_test_tcga_omics")
#' meth_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
#' gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
#' met_integration <- run_met_integration(expression=gene_exp_matrix,
#' methylation=meth_matrix)
#' @export
#' @importFrom BiocParallel bpparam SerialParam

run_met_integration <- function(expression,
                                 methylation,
                                 sequencing_data=TRUE,
                                 normalize=TRUE,
                                 norm_method="TMM",
                                 class=NULL,
                                 run_deg=TRUE,
                                 BPPARAM=SerialParam(),
                                 ...){

  if(is.null(class)){
    met_res <- .def_met_integration(expression = expression,
                                    methylation = methylation,
                                    sequencing_data = sequencing_data,
                                    normalize = normalize,
                                    norm_method = norm_method,
                                    BPPARAM = BPPARAM,
                                    ...)
  }else{
    if(!is.character(class) |
       length(class)!=nrow(expression) |
       !identical(names(class),rownames(expression))){
      stop(str_wrap("class should be a named character vector, names should
                      match sample names contained in expression data"))
    }

    tmp <- unique(class)
    tmp <- lapply(tmp, function(x) names(class)[class==x])
    names(tmp) <- unique(class)
    met_res <- lapply(tmp, function(x){
      ans <- .def_met_integration(expression = expression[x,],
                                  methylation = methylation[x,],
                                  sequencing_data = sequencing_data,
                                  normalize = normalize,
                                  norm_method = norm_method,
                                  BPPARAM = BPPARAM,
                                  ...)
      return(ans)
    })
    deg <- NULL
    if(run_deg){
      deg <- .find_deg(eexpression = t(expression),
                       class = class,
                       normalize = normalize,
                       norm_method = norm_method)
      deg <- list(rownames(deg)[deg$FDR<=0.1])
    }
    met_res <- new("MultiClass", c(met_res, deg=deg))
  }
  return(met_res)
}

# def genomic integration
.def_genomic_integration <- function(expression,
                                    cnv_data,
                                    methylation,
                                    sequencing_data,
                                    normalize,
                                    norm_method,
                                    interactions,
                                    scale,
                                    BPPARAM,
                                    ...){


  if(scale){
    original_cnv <- cnv_data
    original_met <- methylation
    colnames(original_cnv) <- paste0(colnames(original_cnv), "_cnv")
    colnames(original_met) <- paste0(colnames(original_met), "_met")
    cnv_data <- scale(cnv_data)
    methylation <- scale(methylation)
  }

  if(is.null(interactions)){

    tmp <- Reduce(intersect,
                  list(colnames(expression),
                       colnames(cnv_data),
                       colnames(methylation)))
    expression <- expression[, tmp]
    cnv_data <- cnv_data[, tmp]
    methylation <- methylation[, tmp]
    colnames(cnv_data) <- paste0(colnames(cnv_data), "_cnv")
    colnames(methylation) <- paste0(colnames(methylation), "_met")
    message("Generating interactions")
    interactions <- lapply(seq_along(colnames(expression)), function(x)
      c(colnames(cnv_data)[x], colnames(methylation)[x]))
    names(interactions) <- colnames(expression)
  }



  if(sequencing_data==TRUE){


    gen_res <- .run_edgeR_integration(response_var = expression,
                                      covariates = cbind(cnv_data, methylation),
                                      normalize = normalize,
                                      norm_method = norm_method,
                                      interactions = interactions,
                                      BPPARAM = BPPARAM,
                                      ...)
    if(normalize & !is.null(gen_res)){
      gen_res$data$response_var <- .data_norm(gen_res$data$response_var,
                                              method = norm_method)
    }
  }else{
    gen_res <- .run_lm_integration(response_var = expression,
                                   covariates = cbind(cnv_data, methylation),
                                   normalize = normalize,
                                   norm_method = norm_method,
                                   interactions = interactions,
                                   BPPARAM = BPPARAM,
                                   ...)
  }

  if(scale & !is.null(gen_res)){
    tmp <- cbind(original_cnv, original_met[rownames(original_cnv),])
    gen_res$data$covariates <- tmp[rownames(gen_res$data$covariates),
                                   colnames(gen_res$data$covariates)]
  }
  return(gen_res)
}



#' Integration of expression, Copy Number Variations and methylation data
#' @description
#' This function will perform an integration of expression data and Copy Number
#' Variations data
#' @param expression Matrix or data.frame containing the expression values
#' for each model. Rows represent samples, while each column represents
#' the different response variables of the models.
#' @param cnv_data Matrix or data.frame containing the Copy Number variation
#' status for the models. Rows represent samples, while columns represent
#' the different covariates. If **interactions** are not provided, they will be
#' automatically generated and for each gene contained in **expression**
#' the model will look for the same gene in **cnv_data**
#' @param methylation Matrix or data.frame containing the methylation
#' values for the models. Rows represent samples, while columns represent
#' the different covariates. If **interactions** are not provided, they will be
#' automatically generated and for each gene contained in **expression**
#' the model will look for the same gene in **methylation**
#' @param sequencing_data logical. Are expression data obtained from RNA
#' sequencing ? Default is set to TRUE
#' @param normalize logical.Should expression data be
#' normalized ? Default is set to TRUE
#' @param norm_method Normalization method to be used for
#' expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
#' "none".
#' @param interactions A list of character vectors containing the interactions
#' between response variable and covariates. The names of the list should
#' match the response variables while the character contained in each element
#' of the list should match the covariates. If NULL (default), the interactions
#' will be automatically defined according to response variable's colnames.
#' @param class Character vector specifying the classes for differential
#' expression analysis.
#' @param scale Logical. Should the data be scaled? Default is set to TRUE.
#' @param run_deg Logical. Should differential expression analysis be performed?
#'  Default is set to TRUE.
#' @param BPPARAM A BiocParallelParam object specifying the parallel backend to
#'  be used.
#' @param ... Additional arguments to be passed to internal functions.
#' @importFrom plyr rbind.fill
#' @importFrom BiocParallel bpparam SerialParam
#' @return A list or a \linkS4class{MultiClass} object if **class** is provided
#' containing the results of the Genomic integration
#' @examples
#' # Example usage_multi:
#' data("ov_test_tcga_omics")
#' meth_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
#' gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
#' gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
#' genomic_integration <- run_genomic_integration(expression=gene_exp_matrix,
#' cnv_data=gene_cnv_matrix, methylation=meth_matrix)
#' @export
run_genomic_integration <- function(expression,
                                cnv_data,
                                methylation,
                                sequencing_data=TRUE,
                                normalize=TRUE,
                                norm_method="TMM",
                                interactions=NULL,
                                class=NULL,
                                scale=TRUE,
                                run_deg=TRUE,
                                BPPARAM = SerialParam(),
                                ...){


  if(is.null(class)){
    gen_res <- .def_genomic_integration(expression = expression,
                                        cnv_data = cnv_data,
                                        methylation = methylation,
                                        sequencing_data = sequencing_data,
                                        normalize = normalize,
                                        norm_method = norm_method,
                                        interactions = interactions,
                                        scale = scale,
                                        BPPARAM = BPPARAM,
                                        ...)
  }else{
    if(!is.character(class) |
       length(class)!=nrow(expression) |
       !identical(names(class),rownames(expression))){
      stop(str_wrap("class should be a named character vector, names should
                      match sample names contained in expression data"))
    }

    tmp <- unique(class)
    tmp <- lapply(tmp, function(x) names(class)[class==x])
    names(tmp) <- unique(class)
    gen_res <- lapply(tmp, function(x){
      ans <- .def_genomic_integration(expression = expression[x,],
                                      cnv_data = cnv_data[x,],
                                      methylation = methylation[x,],
                                      sequencing_data = sequencing_data,
                                      normalize = normalize,
                                      norm_method = norm_method,
                                      interactions = interactions,
                                      scale = scale,
                                      BPPARAM = BPPARAM,
                                      ...)
      return(ans)
    })
    deg <- NULL
    if(run_deg){
      deg <- .find_deg(eexpression = t(expression),
                       class = class,
                       normalize = normalize,
                       norm_method = norm_method)
      deg <- list(rownames(deg)[deg$FDR<=0.1])
    }
    gen_res <- new("MultiClass", c(gen_res, deg=deg))
  }
  return(gen_res)
}

# def tf integration
.def_tf_integration <- function(expression,
                                tf_expression,
                                interactions,
                                type,
                                sequencing_data,
                                species,
                                normalize,
                                norm_method,
                                normalize_cov,
                                norm_method_cov,
                                BPPARAM,
                                ...){

  if(is.null(interactions)){
    if(!type%in%c("tf_miRNA", "tf", "miRNA_target"))
      stop(str_wrap("If interactions are not provided, type should be one
                      of tf_miRNA, tf, miRNA_target"))
    if(type=="tf_miRNA"){
      interactions <- .download_tf_mirna(miRNAs = colnames(expression),
                                         species = species)
    }

    if(type=="tf"){
      interactions <- .download_tf(genes = colnames(expression),
                                   species = species)
    }

    if(type=="miRNA_target"){
      interactions <- .download_mirna_target(miRNAs = colnames(tf_expression),
                                             species = species)
    }


  }

  if(normalize_cov) tf_expression <- .data_norm(tf_expression,
                                                method = norm_method_cov)


  if(length(interactions)==0) return(NULL)
  if(sequencing_data==TRUE){
    tf_res <- .run_edgeR_integration(response_var = expression,
                                     covariates = tf_expression,
                                     interactions = interactions,
                                     normalize = normalize,
                                     norm_method = norm_method,
                                     BPPARAM = BPPARAM,
                                     ...)
    if(normalize & !is.null(tf_res)){
      tf_res$data$response_var <- .data_norm(tf_res$data$response_var,
                                             method = norm_method)
    }
  }else{
    tf_res <- .run_lm_integration(response_var = expression,
                                  covariates = tf_expression,
                                  interactions = interactions,
                                  normalize = normalize,
                                  norm_method = norm_method,
                                  BPPARAM = BPPARAM,
                                  ...)
  }
  return(tf_res)
}




#' Integration of expression and Transcription Factors / Generic Regulators
#' @description
#' This function will perform an integration of gene/miRNA expression data
#' and Transcription Factors expression. Moreover, every type of regulator can
#' be provided to the function as covariate through the **tf_expression**
#' argument. Interactions for TF-target, miRNA-target and TF-miRNA integration
#' will be automatically downloaded by the function as defined by the **type**
#' argument. Other types of interactions should be provided through the
#' **interactions** argument.
#' @param expression Matrix or data.frame containing the expression values
#' for each model. Rows represent samples, while each column represents
#' the different response variables of the models.
#' @param tf_expression Matrix or data.frame containing the expression
#' values for the models. Rows represent samples, while columns represent
#' the different covariates. If not provided, it will be set equal to
#' **expression**.
#' @param interactions A list of character vectors containing the interactions
#' between response variable and covariates. The names of the list should
#' match the response variables while the character contained in each element
#' of the list should match the covariates. If NULL (default), the interactions
#' will be automatically downloaded according to the **type** argument.
#' @param type A character defining the type of regulation under analysis.
#' Should be one of "tf_miRNA", "tf", "miRNA_target".
#' @param species species information for interactions download. Fully
#' supported species are "hsa" (default) and "mmu".
#' @param sequencing_data logical. Are expression data obtained from RNA
#' sequencing ? Default is set to TRUE
#' @param normalize logical.Should expression data be
#' normalized ? Default is set to TRUE
#' @param norm_method Normalization method to be used for
#' expression data. One of "TMM" (default), "TMMwsp", "RLE", "upperquartile",
#' "none".
#' @param normalize_cov Same as **normalize** but for covariates.
#' @param norm_method_cov Same as **norm_method** but for covariates.
#' @param class Character vector specifying the classes for differential
#' expression analysis.
#' @param run_deg Logical. Should differential expression analysis be performed?
#'  Default is set to TRUE.
#' @param BPPARAM A BiocParallelParam object specifying the parallel backend to
#'  be used.
#' @param ... Additional arguments to be passed to internal functions.
#' @importFrom stats quantile
#' @importFrom BiocParallel bpparam SerialParam
#' @return A list or a \linkS4class{MultiClass} object if **class** is provided
#' containing the results of the transcriptional integration
#' @examples
#' # Example usage_multi:
#' data("ov_test_tcga_omics")
#' gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
#' tf_integration <- run_tf_integration(expression=gene_exp_matrix)
#' @export

run_tf_integration <- function(expression,
                                tf_expression=expression,
                                interactions=NULL,
                                type="none",
                                sequencing_data=TRUE,
                                species="hsa",
                                normalize=TRUE,
                                norm_method="TMM",
                                normalize_cov=TRUE,
                                norm_method_cov="TMM",
                                class=NULL,
                                run_deg=TRUE,
                                BPPARAM=SerialParam(),
                                ...){


  if(is.null(class)){
    tf_res <- .def_tf_integration(expression=expression,
                                  tf_expression=tf_expression,
                                  interactions=interactions,
                                  type=type,
                                  sequencing_data=sequencing_data,
                                  species=species,
                                  normalize=normalize,
                                  norm_method=norm_method,
                                  normalize_cov=normalize_cov,
                                  norm_method_cov=norm_method_cov,
                                  BPPARAM = BPPARAM,
                                  ...)
  }else{
    if(!is.character(class) |
       length(class)!=nrow(expression) |
       !identical(names(class),rownames(expression))){
      stop(str_wrap("class should be a named character vector, names should
                      match sample names contained in expression data"))
    }

    tmp <- unique(class)
    tmp <- lapply(tmp, function(x) names(class)[class==x])
    names(tmp) <- unique(class)
    tf_res <- lapply(tmp, function(x){
      ans <- .def_tf_integration(expression=expression[x,],
                                 tf_expression=tf_expression[x,],
                                 interactions=interactions,
                                 type=type,
                                 sequencing_data=sequencing_data,
                                 species=species,
                                 normalize=normalize,
                                 norm_method=norm_method,
                                 normalize_cov=normalize_cov,
                                 norm_method_cov=norm_method_cov,
                                 BPPARAM = BPPARAM,
                                 ...)
      return(ans)
    })
    deg <- NULL
    if(run_deg){
      deg <- .find_deg(eexpression = t(expression),
                       class = class,
                       normalize = normalize,
                       norm_method = norm_method)
      deg <- list(rownames(deg)[deg$FDR<=0.1])
    }
    tf_res <- new("MultiClass", c(tf_res, deg=deg))
  }
  return(tf_res)
}


