setClass("MultiOmics",
         representation("list")
)

#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @export
run_multiomics <- function(data=NULL,
                           use_cnv_res=T,
                           interactions_met=NULL,
                           interactions_miRNA=NULL,
                           interactions_tf=NULL,
                           interactions_regulators=NULL,
                           RNAseq=T,
                           miRNAseq=T,
                           normalize_miRNA_expr=T,
                           normalize_gene_expr=T,
                           norm_method_gene_expr="TMM",
                           norm_method_miRNA_expr="TMM",
                           ...){


#si lavora con multiassay experiment, va creato con la funzione create_multia..
#prima integro il copy number e l'espressione e do la possibilità
#di usare i residui per le analisi successive
#prima scelta, se l'espressione è un dato di conte usiamo edgeR altrimenti
#usiamo modello lineare

    if(class(data)!="MultiAssayExperiment"){
      stop(str_wrap("You need to provide a MultiassayExperiment as input
                    data. The names of the experiments included in the
                    MultiAssay should follow the rule provided in the
                    help section. You can use the function create_multiassay
                    to create proper input data"))
    }


  normalize_gene_expr2 <- normalize_gene_expr
  normalize_miRNA_expr2 <- normalize_miRNA_expr
  if(!is.null(data@ExperimentList$cnv_data) &
       !is.null(data@ExperimentList$gene_exp)){
        gene_cnv_res <- run_cnv_integration(
          expression = t(assay(data, i = "gene_exp")),
          cnv_data = t(assay(data, i = "cnv_data")),
          sequencing_data = RNAseq,
          normalize=normalize_gene_expr,
          norm_method=norm_method_gene_expr,
          BPPARAM=SerialParam())
      data <- c(data, gene_exp_original=data@ExperimentList$gene_exp)
      data@ExperimentList$gene_exp <- data@ExperimentList$gene_exp[
        rownames(gene_cnv_res$residuals),
        colnames(gene_cnv_res$residuals)]
      assay(data@ExperimentList$gene_exp) <- gene_cnv_res$residuals[
        rownames(assay(data, i = "gene_exp")),
        colnames(assay(data, i = "gene_exp"))]
      RNAseq <- F
      normalize_gene_expr2 <- F
    }


    if(!is.null(data@ExperimentList$miRNA_cnv_data) &
       !is.null(data@ExperimentList$miRNA_exp)){
      mirna_cnv_res <- run_cnv_integration(
        expression = t(assay(data, i = "miRNA_exp")),
        cnv_data = t(assay(data, i = "miRNA_cnv_data")),
        sequencing_data = miRNAseq,
        normalize=normalize_miRNA_expr,
        norm_method=norm_method_miRNA_expr,
        BPPARAM=SerialParam())
      data <- c(data, miRNA_exp_original=data@ExperimentList$miRNA_exp)
      data@ExperimentList$miRNA_exp <- data@ExperimentList$miRNA_exp[
        rownames(mirna_cnv_res$residuals),
        colnames(mirna_cnv_res$residuals)]
      assay(data@ExperimentList$miRNA_exp) <- mirna_cnv_res$residuals[
        rownames(assay(data, i = "miRNA_exp")),
        colnames(assay(data, i = "miRNA_exp"))]
      miRNAseq <- F
      normalize_miRNA_expr2 <- F

    }

  if(!is.null(data@ExperimentList$methylation) &
     !is.null(data@ExperimentList$miRNA_exp)){
    met_res <- run_met_integration(
      expression = t(assay(data, i = "gene_exp")),
      methylation = t(assay(data, i = "methylation")),
      sequencing_data = RNAseq,
      normalize = normalize_gene_expr2)

  }



  if(!is.null(data@ExperimentList$miRNA_exp) &
     !is.null(data@ExperimentList$gene_exp)){

    tmp <- t(assay(data,i = "gene_exp_original"))
    if(normalize_gene_expr==T) tmp <- data_norm(tmp)

    if(!is.null(data@ExperimentList$gene_exp_original)){
      tf_mirna_res <- run_tf_integration(
        expression = t(assay(data, i = "miRNA_exp")),
        tf_expression = tmp,
        interactions = interactions_tf,
        sequencing_data = RNAseq, normalize=normalize_miRNA_expr2)
    }else{
      tf_mirna_res <- run_tf_integration(
        expression = t(assay(data, i = "miRNA_exp")),
        tf_expression = tmp,
        interactions = interactions_tf,
        sequencing_data = RNAseq, normalize=normalize_miRNA_expr2)
      }

  }

#proseguo con le altre integrazioni

  ans <- new("MultiOmics",  list(gene_cnv_res=gene_cnv_res,
              mirna_cnv_res=mirna_cnv_res,
              met_res=met_res,
              tf_mirna_res=tf_mirna_res))
  return(ans)
}

####################################################################
#' @export
run_cnv_integration <- function(expression,
                                cnv_data,
                                sequencing_data=T,
                                ...){

  if(sequencing_data==T){
    cnv_res <- run_edgeR_integration(response_var = expression,
                                     covariates = cnv_data,
                                     ...)
  }else{
    cnv_res <- run_lm_integration(response_var = expression,
                                  covariates = cnv_data,
                                  ...)
  }
  return(cnv_res)
}


####################################################################
#' @export

run_met_integration <- function( expression,
                                 methylation,
                                 sequencing_data=F,
                                 interactions="auto",
                                 normalize=F,
                                 ...){


  if(sequencing_data==T){
    met_res <- run_edgeR_integration(response_var = expression,
                                     covariates = methylation,
                                     interactions = interactions,
                                     normalize = normalize,
                                     ...)
  }else{
    met_res <- run_lm_integration(response_var = expression,
                                  covariates = methylation,
                                  interactions = interactions,
                                  normalize = normalize,
                                  ...)
  }
  return(met_res)
}


####################################################################
#' @export

run_tf_integration <- function( expression,
                                tf_expression=expression,
                                interactions=NULL,
                                sequencing_data=T,
                                species="hsa",
                                ...){

    if(is.null(interactions)){
      interactions <- integrazione::tf_mirna[[species]]
      interactions <- interactions[interactions$level=="literature",]
      tmp <- unique(interactions$miRNA)
      interactions <- lapply(tmp, function(x){
        interactions$TF[grep(x,interactions$miRNA)]
      })
      names(interactions) <- tmp

    }


  if(sequencing_data==T){
    tf_res <- run_edgeR_integration(response_var = expression,
                                    covariates = tf_expression,
                                    interactions = interactions,
                                    ...)
  }else{
    tf_res <- run_lm_integration(response_var = expression,
                                 covariates = tf_expression,
                                 interactions = interactions,
                                 ...)
  }
  return(tf_res)
}


