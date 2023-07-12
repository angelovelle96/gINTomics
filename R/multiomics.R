
#' @import MultiAssayExperiment SummarizedExperiment
#' @export
run_multiomics <- function(data,
                           use_cnv_res=T,
                           interactions_met=NULL,
                           interactions_miRNA_target=NULL,
                           interactions_tf=NULL,
                           interactions_tf_miRNA=NULL,
                           interactions_regulators=NULL,
                           RNAseq=T,
                           miRNAseq=T,
                           normalize_miRNA_expr=T,
                           normalize_gene_expr=T,
                           norm_method_gene_expr="TMM",
                           norm_method_miRNA_expr="TMM",
                           ...){



    if(class(data)!="MultiAssayExperiment"){
      stop(str_wrap("You need to provide a MultiassayExperiment as input
                    data. The names of the experiments included in the
                    MultiAssay should follow the rule provided in the
                    help section. You can use the function create_multiassay
                    to create proper input data"))
    }


    normalize_gene_expr2 <- normalize_gene_expr
    normalize_miRNA_expr2 <- normalize_miRNA_expr
    data <- c(data, gene_exp_original=data@ExperimentList$gene_exp)
    gene_cnv_res <- NULL
    if(!is.null(data@ExperimentList$cnv_data) &
         !is.null(data@ExperimentList$gene_exp)){
          gene_cnv_res <- run_cnv_integration(
            expression = t(assay(data, i = "gene_exp")),
            cnv_data = t(assay(data, i = "cnv_data")),
            sequencing_data = RNAseq,
            normalize=normalize_gene_expr,
            norm_method=norm_method_gene_expr,
            BPPARAM=SerialParam())
        data@ExperimentList$gene_exp <- data@ExperimentList$gene_exp[
          rownames(gene_cnv_res$residuals),
          colnames(gene_cnv_res$residuals)]
        assay(data@ExperimentList$gene_exp) <- gene_cnv_res$residuals[
          rownames(assay(data, i = "gene_exp")),
          colnames(assay(data, i = "gene_exp"))]
        RNAseq <- F
        normalize_gene_expr2 <- F
      }


    data <- c(data, miRNA_exp_original=data@ExperimentList$miRNA_exp)
    mirna_cnv_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_cnv_data) &
       !is.null(data@ExperimentList$miRNA_exp)){
      mirna_cnv_res <- run_cnv_integration(
        expression = t(assay(data, i = "miRNA_exp")),
        cnv_data = t(assay(data, i = "miRNA_cnv_data")),
        sequencing_data = miRNAseq,
        normalize=normalize_miRNA_expr,
        norm_method=norm_method_miRNA_expr,
        BPPARAM=SerialParam())
      data@ExperimentList$miRNA_exp <- data@ExperimentList$miRNA_exp[
        rownames(mirna_cnv_res$residuals),
        colnames(mirna_cnv_res$residuals)]
      assay(data@ExperimentList$miRNA_exp) <- mirna_cnv_res$residuals[
        rownames(assay(data, i = "miRNA_exp")),
        colnames(assay(data, i = "miRNA_exp"))]
      miRNAseq <- F
      normalize_miRNA_expr2 <- F

    }

    met_res <- NULL
    if(!is.null(data@ExperimentList$methylation) &
       !is.null(data@ExperimentList$gene_exp)){
      met_res <- run_met_integration(
        expression = t(assay(data, i = "gene_exp")),
        methylation = t(assay(data, i = "methylation")),
        sequencing_data = RNAseq,
        normalize = normalize_gene_expr2,
        norm_method=norm_method_gene_expr)
      tmp <- t(assay(data, i = "gene_exp_original"))
      tmp <- tmp[rownames(met_res$data$response_var),
                 colnames(met_res$data$response_var)]
      if(normalize_gene_expr) tmp <-  data_norm(tmp,
                                                method = norm_method_gene_expr)
      met_res$data$response_var <- tmp

    }


    tf_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_exp)){

      tf_res <- run_tf_integration(
        expression = t(assay(data,i = "gene_exp")),
        tf_expression = t(assay(data,i = "gene_exp_original")),
        interactions = interactions_tf,
        sequencing_data = RNAseq,
        normalize=normalize_gene_expr2,
        norm_method=norm_method_gene_expr,
        normalize_cov = normalize_gene_expr,
        norm_method_cov = norm_method_gene_expr,
        type="tf")
      tmp <- t(assay(data, i = "gene_exp_original"))
      tmp <- tmp[rownames(tf_res$data$response_var),
                 colnames(tf_res$data$response_var)]
      if(normalize_gene_expr) tmp <-  data_norm(tmp,
                                                method = norm_method_gene_expr)
      tf_res$data$response_var <- tmp
    }


    tf_mirna_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_exp) &
       !is.null(data@ExperimentList$gene_exp)){

        tf_mirna_res <- run_tf_integration(
          expression = t(assay(data,i = "miRNA_exp")),
          tf_expression = t(assay(data,i = "gene_exp_original")),
          interactions = interactions_tf_miRNA,
          sequencing_data = miRNAseq,
          normalize=normalize_miRNA_expr2,
          norm_method=norm_method_miRNA_expr,
          normalize_cov = normalize_gene_expr,
          norm_method_cov = norm_method_gene_expr,
          type="tf_miRNA")
        tmp <- t(assay(data, i = "miRNA_exp_original"))
        tmp <- tmp[rownames(tf_mirna_res$data$response_var),
                   colnames(tf_mirna_res$data$response_var)]
        if(normalize_miRNA_expr) tmp <-  data_norm(tmp,
                                                  method = norm_method_miRNA_expr)
        tf_mirna_res$data$response_var <- tmp
    }


    mirna_target_res <- NULL
    if(!is.null(data@ExperimentList$miRNA_exp) &
       !is.null(data@ExperimentList$gene_exp)){

      mirna_target_res <- run_tf_integration(
        expression = t(assay(data,i = "gene_exp")),
        tf_expression = t(assay(data,i = "miRNA_exp_original")),
        interactions = interactions_miRNA_target,
        sequencing_data = RNAseq,
        normalize=normalize_gene_expr2,
        norm_method=norm_method_gene_expr,
        normalize_cov = normalize_miRNA_expr,
        norm_method_cov = norm_method_miRNA_expr,
        type="miRNA_target")
      tmp <- t(assay(data, i = "gene_exp_original"))
      tmp <- tmp[rownames(mirna_target_res$data$response_var),
                 colnames(mirna_target_res$data$response_var)]
      if(normalize_gene_expr) tmp <-  data_norm(tmp,
                                                 method = norm_method_gene_expr)
      mirna_target_res$data$response_var <- tmp
    }





#proseguo con le altre integrazioni

  ans <- new("MultiOmics", Filter(Negate(is.null),
                                  list(gene_cnv_res=gene_cnv_res,
                                       mirna_cnv_res=mirna_cnv_res,
                                       met_res=met_res,
                                       tf_res=tf_res,
                                       tf_mirna_res=tf_mirna_res,
                                       mirna_target_res=mirna_target_res)))
  return(ans)
}

####################################################################
#' @export
run_cnv_integration <- function(expression,
                                cnv_data,
                                sequencing_data=T,
                                normalize=T,
                                norm_method="TMM",
                                ...){

  if(sequencing_data==T){
    cnv_res <- run_edgeR_integration(response_var = expression,
                                     covariates = cnv_data,
                                     normalize = normalize,
                                     norm_method = norm_method,
                                     ...)
    if(normalize){
      cnv_res$data$response_var <- data_norm(cnv_res$data$response_var,
                                             method = norm_method)
      }
  }else{
    cnv_res <- run_lm_integration(response_var = expression,
                                  covariates = cnv_data,
                                  normalize = normalize,
                                  norm_method = norm_method,
                                  ...)
  }
  return(cnv_res)
}


####################################################################
#' @export

run_met_integration <- function( expression,
                                 methylation,
                                 sequencing_data=T,
                                 normalize=T,
                                 norm_method="TMM",
                                 ...){

  if(sequencing_data==T){
    met_res <- run_edgeR_integration(response_var = expression,
                                     covariates = methylation,
                                     normalize = normalize,
                                     norm_method = norm_method,
                                     ...)
    if(normalize){
      met_res$data$response_var <- data_norm(met_res$data$response_var,
                                             method = norm_method)
    }
  }else{
    met_res <- run_lm_integration(response_var = expression,
                                  covariates = methylation,
                                  normalize = normalize,
                                  norm_method = norm_method,
                                  ...)
  }
  return(met_res)
}


####################################################################
#' @export

run_tf_integration <- function( expression,
                                tf_expression=expression,
                                interactions=NULL,
                                type=NULL,
                                sequencing_data=T,
                                species="hsa",
                                normalize=T,
                                norm_method="TMM",
                                normalize_cov=T,
                                norm_method_cov="TMM",
                                ...){

    if(is.null(interactions)){
      if(type=="tf_miRNA"){
        interactions <- tf_mirna[[species]]
        interactions <- interactions[interactions$level=="literature",]
        tmp <- unique(interactions$miRNA)
        interactions <- lapply(tmp, function(x){
          interactions$TF[grep(x,interactions$miRNA)]
        })
        names(interactions) <- tmp
      }

      if(type=="tf"){
        interactions <- download_tf(genes = colnames(expression),
                                    species = species)
      }

      if(type=="miRNA_target"){
        interactions <- download_mirna_target(miRNAs = colnames(tf_expression),
                                              targets = colnames(expression),
                                              species = species)
        interactions <- lapply(interactions, function(x)
          gsub("-miR-", "-mir-", x))
      }


    }

    if(normalize_cov) tf_expression <- data_norm(tf_expression,
                                                 method = norm_method_cov)

    tmp <- unlist(lapply(interactions, length))
    if(quantile(tmp, 0.75)>=10) sequencing_data <- F

    if(sequencing_data==T){
      tf_res <- run_edgeR_integration(response_var = expression,
                                      covariates = tf_expression,
                                      interactions = interactions,
                                      normalize = normalize,
                                      norm_method = norm_method,
                                      ...)
      if(normalize){
        tf_res$data$response_var <- data_norm(tf_res$data$response_var,
                                               method = norm_method)
      }
    }else{
      tf_res <- run_lm_integration(response_var = expression,
                                   covariates = tf_expression,
                                   interactions = interactions,
                                   normalize = normalize,
                                   norm_method = norm_method,
                                   ...)
    }
    return(tf_res)
}


