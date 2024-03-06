#########################
###############################
plot_network <- function(data_table,
                         num_interactions=300,
                         class=NULL,
                         pval=0.05){

  if(!is.null(class)) data_table <- data_table[data_table$class==class,]
  nnet <- .prepare_network(data_table)
  nnet$edges <- nnet$edges[nnet$edges$pval<=pval,]
  nnet$edges <- nnet$edges[seq_len(length.out = num_interactions),]
  nodes_with_edges <- unique(c(nnet$edges$from, nnet$edges$to))
  nnet$nodes <- nnet$nodes[nnet$nodes$id %in% nodes_with_edges,]
  .build_network(nodes = nnet$nodes,
                 edges = nnet$edges,
                 legend_nodes = nnet$legend_nodes,
                 legend_edges = nnet$legend_edges)

}

#########################
###############################


plot_venn <- function(data_table,
                      class=NULL){
  if(is.null(class) & "class"%in%colnames(data_table))
     stop(str_wrap("Please, specify class"))
  input <- list()
  input$classSelectVenn <- class
  input$fdrRangeVenn <- c(0, 0.05)
  input$pvalRangeVenn <- c(0,0.05)
  input$significativityCriteriaVenn <- "pval"
  reactive_venn <- gINTomics:::.prepare_reactive_venn(data_table = data_table,
                                                      input = input,
                                                      output = output,
                                                      deg = FALSE)
  tmp <- isolate(reactive_venn())
  ans <- .build_venn(tmp)
  ggplotly(ans)
}


#########################
###############################


plot_volcano <- function(data_table,
                         class=NULL,
                         omics=NULL){
  if(is.null(class) & "class"%in%colnames(data_table))
    stop(str_wrap("Please, specify class"))
  if(is.null(omics) & length(unique(data_table$omics))>0)
    stop(str_wrap("Please, specify omics"))

  if(!is.null(class)) data_table <- data_table[data_table$class==class,]
  if(!is.null(omics)) data_table <- data_table[data_table$omics==omics,]
  data_table["group"] <- "Not Significant"
  data_table[data_table$pval <= 0.05, 'group'] <- "Significant"
  data_table$pval_fdr <- -log10(data_table$fdr)
  .build_volcano(data_table)
}


#########################
###############################


plot_ridge <- function(data_table,
                       class=NULL,
                       omics=NULL,
                       cnv_met=NULL){
  if(is.null(class) & "class"%in%colnames(data_table))
    stop(str_wrap("Please, specify class"))
  if(is.null(omics) & length(unique(data_table$omics))>0)
    stop(str_wrap("Please, specify omics"))

  if(!is.null(class)) data_table <- data_table[data_table$class==class,]
  if(!is.null(omics)) data_table <- data_table[data_table$omics==omics,]
  data_table$significance <- ifelse(data_table$pval <= 0.05,
                            "Significant", "Not Significant")
  if(omics == "gene_genomic_res" & !is.null(cnv_met)){
    if(!cnv_met%in%c("cnv", "met"))
      stop(str_wrap("cnv_met should be one of cnv and met"))
    data_table <- data_table[data_table$cnv_met == cnv_met,]
    data_table <- data_table[!is.na(data_table$cnv_met), ]
  }
  lower_quantile <- quantile(data_table$coef, 0.001)
  upper_quantile <- quantile(data_table$coef, 0.999)
  quantiles = c(lower_quantile, upper_quantile)
  .build_ridge(ridge_data=data_table,quantiles=quantiles)
}


#########################
###############################


plot_heatmap <- function(data_table,
                         num_interactions=300,
                         omics,
                         class=NULL,
                         pval=0.05){

  if (omics == "gene_genomic_res"){
    ans <- .prepare_gen_heatmap(data_table = data_table,
                               df_heatmap = df_heatmap,
                               df_heatmap_t = df_heatmap_t,
                               significativityCriteria=significativityCriteria,
                               pvalRange = pvalRange,
                               fdrRange = fdrRange,
                               numTopCNV = numTopCNV,
                               numTopMET = numTopMET)
  }
  if(omics == "gene_cnv_res"){
    ans <- .prepare_cnv_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria=significativityCriteria,
                              pvalRange = pvalRange,
                              fdrRange = fdrRange,
                              numTopCNVonly = numTopCNVonly)
  }
  if(omics == "gene_met_res"){
    ans <-  .prepare_met_heatmap(data_table = data_table,
                               df_heatmap = df_heatmap,
                               df_heatmap_t = df_heatmap_t,
                               significativityCriteria=significativityCriteria,
                               pvalRange = pvalRange,
                               fdrRange = fdrRange,
                               numTopMETonly = numTopMETonly)
  }
  if(omics == "mirna_cnv_res"){
    ans <- .prepare_mirna_heatmap(data_table = data_table,
                               df_heatmap = df_heatmap,
                               df_heatmap_t = df_heatmap_t,
                               significativityCriteria=significativityCriteria,
                               pvalRange = pvalRange,
                               fdrRange = fdrRange,
                               numTopMiCNV = numTopMiCNV)
  }
  return(ht)

}








