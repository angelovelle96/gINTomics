##############################
###############################
#' Plotting network
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

#' plotting venn
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

#' plotting volcano
plot_volcano <- function(data_table,
                         class=NULL,
                         omics=NULL,
                         cnv_met=NULL){
  if(!is.null(class) & !"class"%in%colnames(data_table))
    stop(str_wrap("Class not available"))
  if(is.null(omics) & length(unique(data_table$omics))>0)
    stop(str_wrap("Please, specify omics"))

  if(!is.null(class)) data_table <- data_table[data_table$class==class,]
  if(!is.null(omics)) data_table <- data_table[data_table$omics==omics,]
  if(omics == "gene_genomic_res" & !is.null(cnv_met)){
    if(!cnv_met%in%c("cnv", "met"))
      stop(str_wrap("cnv_met should be one of cnv and met"))
    data_table <- data_table[data_table$cnv_met == cnv_met,]
    data_table <- data_table[!is.na(data_table$cnv_met), ]
  }
  data_table["group"] <- "Not Significant"
  data_table[data_table$pval <= 0.05, 'group'] <- "Significant"
  data_table$pval_fdr <- -log10(data_table$pval)
  .build_volcano(data_table)
}


#########################
###############################

#' plotting ridge
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
#' plotting heatmap
plot_heatmap <- function(multiomics_integration,
                         data_table,
                         omics,
                         scale="none",
                         genes_number=50,
                         class=NULL,
                         pval=0.05){

  if(!"omics"%in%colnames(data_table)) data_table$omics <- omics
  if(is.null(class) & "class"%in%colnames(data_table))
    stop(str_wrap("Please, specify class"))
  tmp <- c("gene_genomic_res", "gene_cnv_res", "gene_met_res", "mirna_cnv_res")
  if(!omics%in%tmp) stop(paste("Omics should be one of",
                                paste(tmp, collapse = ", ")))
  if(class(multiomics_integration)!="MultiOmics"){
    tmp <- list()
    tmp[[omics]] <- multiomics_integration
    multiomics_integration <- tmp
  }

  if("class"%in%colnames(data_table)){
    df_heatmap <- multiomics_integration[[omics]][[
      class]]$data$response_var
    data_table <- data_table[data_table$omics == omics,]
    data_table <- data_table[data_table$class == class,]
  }else{
    df_heatmap <- multiomics_integration[[
      omics]]$data$response_var
    data_table <- data_table[data_table$omics == omics,]
  }

  df_heatmap_t <- t(as.matrix(df_heatmap))
  if (omics == "gene_genomic_res"){
    ans <- .prepare_gen_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria="pval",
                              pvalRange = pval,
                              fdrRange = pval,
                              numTopCNV = round(genes_number/2),
                              numTopMET = round(genes_number/2),
                              scale = scale)
  }
  if(omics == "gene_cnv_res"){
    ans <- .prepare_cnv_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria="pval",
                              pvalRange = pval,
                              fdrRange = pval,
                              numTopCNVonly = genes_number,
                              scale = scale)
  }
  if(omics == "gene_met_res"){
    ans <-  .prepare_met_heatmap(data_table = data_table,
                               df_heatmap = df_heatmap,
                               df_heatmap_t = df_heatmap_t,
                               significativityCriteria="pval",
                               pvalRange = pval,
                               fdrRange = pval,
                               numTopMETonly = genes_number,
                               scale = scale)
  }
  if(omics == "mirna_cnv_res"){
    ans <- .prepare_mirna_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria="pval",
                              pvalRange = pval,
                              fdrRange = pval,
                              numTopMiCNV = genes_number,
                              scale = scale)
  }
  return(ans)
}


#########################
###############################
#' plotting chr distribution
plot_chr_distribution <- function(data_table,
                                  class=NULL,
                                  omics=NULL,
                                  cnv_met=NULL,
                                  pval=0.05){
  if(is.null(class) & "class"%in%colnames(data_table))
    stop(str_wrap("Please, specify class"))
  if(is.null(omics) & length(unique(data_table$omics))>0)
    stop(str_wrap("Please, specify omics"))
  data_table <- data_table[!is.na(data_table$chr_cov),]
  data_table <- data_table[!is.na(data_table$chr_response),]
  chr_order <- mixedsort(unique(data_table$chr_cov))
  chr_order <- chr_order[!is.na(chr_order)]
  data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)

  if(!is.null(class)) data_table <- data_table[data_table$class==class,]
  if(!is.null(omics)) data_table <- data_table[data_table$omics==omics,]
  if(omics == "gene_genomic_res" & !is.null(cnv_met)){
    if(!cnv_met%in%c("cnv", "met"))
      stop(str_wrap("cnv_met should be one of cnv and met"))
    data_table <- data_table[data_table$cnv_met == cnv_met,]
    data_table <- data_table[!is.na(data_table$cnv_met), ]
  }
  data_table$significance <- ifelse(data_table$pval <= pval,
                                    "Significant",
                                    "Not Significant")
  if (nrow(data_table) == 0) return(NULL)
  .build_histo(data_table)

}

  #########################
  ###############################
#' plotting tf distribution
  plot_tf_distribution <- function(data_table,
                                    class=NULL,
                                    pval=0.05){
    omics="tf_res"
    if(is.null(class) & "class"%in%colnames(data_table))
      stop(str_wrap("Please, specify class"))
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    if("omics"%in%colnames(data_table))
      data_table <- data_table[data_table$omics==omics,]
    data_table <- data_table[,
                             colnames(data_table)%in%c("response", "cov",
                                                       "pval","fdr","chr_cov",
                                                       "deg","class")]
    if(!is.null(class)) data_table <- data_table[data_table$class==class,]
    df_filtered_histo_tf <- data_table
    df_filtered_histo_tf <- df_filtered_histo_tf[df_filtered_histo_tf$pval <= pval, ]
    genes_count <- table(df_filtered_histo_tf$cov, df_filtered_histo_tf$chr_cov)
    genes_count_df <- as.data.frame.table(genes_count)
    genes_count_df <- subset(genes_count_df, Freq != 0)
    colnames(genes_count_df) <- c("TF", "Chromosome", "Count")
    genes_count_df <- genes_count_df[order(-genes_count_df$Count), ]
    .build_histo_TFbyChr(genes_count_df)

}

  ####################################################################
  #########################################################################
#' plotting enrichment
#' @importFrom ggtree fortify
#' @importFrom plotly add_markers subplot plot_ly
  dot_plotly <- function(enrich_result, showCategory=10, width=800, height=700){
    df <- fortify(enrich_result, showCategory = showCategory)
    df$Description <- as.character(df$Description)
    df <- df[order(df$GeneRatio, decreasing = TRUE),]
    df$Description <- unlist(lapply(df$Description, function(label) {
      words <- strsplit(label, " ")[[1]]
      split_words <- lapply(seq(1, length(words), by = 2), function(i) {
        paste(words[i:min(i+1, length(words))], collapse = " ")
      })
      paste(split_words, collapse = "<br>")
    }))
    legend.sizes = seq(min(df$Count),
                       max(df$Count),
                       max(c(1,round(((max(df$Count)-min(df$Count))/4)))))
    lprop <- c(0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55)
    lprop <- lprop[length(legend.sizes)]
    ax = list(zeroline = FALSE,
              showline = FALSE,
              showticklabels = TRUE,
              showgrid = FALSE,
              side="top")
    mk = list(sizeref=0.1, sizemode="area")

    pplot <- plot_ly(df, width=width, height=height) %>%
      add_markers(x=~GeneRatio,
                  y=~Description,
                  name = "DotPlot",
                  text = ~paste("Count:", Count,"<br>",
                                "pValue:", round(pvalue, digits = 4),"<br>",
                                "qValue:", round(qvalue, digits = 4)),
                  size=~Count,
                  color=~pvalue,
                  marker=mk,
                  type="scatter",
                  mode="markers")%>%
      layout(yaxis=list(automargin=TRUE,
                        tickfont = list(size = 7),
                        categoryorder = "array",
                        categoryarray = rev(df$Description)),
             xaxis=list(title="GeneRatio"))

    llegend <- plot_ly() %>%
      add_markers(x = "Count",
                  y = legend.sizes,
                  size = legend.sizes,
                  showlegend = FALSE,
                  fill = ~'',
                  marker = c(mk, color="black"),
                  text=legend.sizes,
                  hoverinfo = "text") %>%
      layout(xaxis = ax,
             yaxis = list(showgrid = FALSE, tickvals = legend.sizes))

    empty_trace <- plot_ly(x = numeric(0),
                           y = numeric(0),
                           type = "scatter",
                           mode = "markers",
                           showlegend = FALSE) %>%
      layout(xaxis = list(showgrid = FALSE,
                          zeroline = FALSE,
                          showticklabels = FALSE),
             yaxis = list(showgrid = FALSE,
                          zeroline = FALSE,
                          showticklabels = FALSE))

    ans <- subplot(empty_trace,llegend,empty_trace,
                   heights = c(0.02,lprop, (0.98-lprop)), nrows = 3)
    ans <- subplot(pplot, ans,
                   widths = c(0.9, 0.1), titleX = TRUE, shareX = FALSE)
    return(ans)

  }


