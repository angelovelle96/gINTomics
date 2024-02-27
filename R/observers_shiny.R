
.render_reactive_network <- function(reactive_network,
                                     input,
                                     output){
  observe({
  output$networkPlot <- renderVisNetwork({
    network_data <- reactive_network()
    network <- .build_network(nodes=network_data$nodes,
                              edges=network_data$edges,
                              legend_nodes=network_data$legend_nodes,
                              legend_edges=network_data$legend_edges)
    network <- network%>%
      visHierarchicalLayout(enabled = input$layoutNetwork)%>%
      visIgraphLayout(physics = input$physics,
                      layout = "layout.fruchterman.reingold",
                      randomSeed = 20)

    return(network)
    })
  })%>%bindEvent(input$layoutNetwork,
                 input$physics)
}

################################################################
#################################################################

.select_deg_network <- function(data_table,
                                network_data,
                                input,
                                output){
  reactive({
    if(input$degNetwork==FALSE){
    ans <- network_data
    data_table <- data_table[data_table$cov != "(Intercept)",]
    if(input$significativityCriteriaNetwork == "pval"){
      data_table <- data_table[data_table$pval <= input$pvalNetwork,]
    } else {
      data_table <- data_table[data_table$fdr <= input$fdrNetwork,]
    }
    ssel <- unique(data_table$response)
    ssel <- head(ssel, input$numNodes)
    ans$nodes <- ans$nodes[ans$nodes$id%in%ssel,]
    ans$edges <- ans$edges[ans$edges$from%in%ssel,]
    ans$edges <- ans$edges[ans$edges$to%in%ssel,]
    nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
    ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges,]
    }

    if(input$degNetwork==TRUE){
      ans <- network_data
      data_table <- data_table[data_table$cov != "(Intercept)",]
      if(input$significativityCriteriaNetwork == "pval"){
        data_table <- data_table[data_table$pval <= input$pvalNetwork,]
      } else {
        data_table <- data_table[data_table$fdr <= input$fdrNetwork,]
      }
      ssel <- unique(data_table$response[data_table$deg])
      ssel <- head(ssel, input$numNodes)
      ans$nodes <- ans$nodes[ans$nodes$id%in%ssel,]
      ans$edges <- ans$edges[ans$edges$from%in%ssel,]
      ans$edges <- ans$edges[ans$edges$to%in%ssel,]
      # nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
      # ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges,]
    }

    return(ans)
  })%>%bindEvent(input$layoutNetwork,
                 input$numNodes,
                 input$significativityCriteriaNetwork,
                 input$pvalNetwork,
                 input$fdrNetwork,
                 input$degNetwork,
                 input$physics)
}

################################################################
#################################################################
.prepare_reactive_venn <- function(data_table,
                                   input,
                                   output,
                                   type = "genomic"){
  reactive_venn <- reactive({
    if(type=="genomic") {
      data_table <- data_table[data_table$omics == 'gene_genomic_res', ]
    }else{
      return(NULL)}
    if('class' %in% names(data_table)){data_table <- data_table[data_table$class == input$classSelectVenn, ]}

    if(input$significativityCriteriaVenn == 'pval'){

      cnv_sign_genes <- subset(data_table,
                               cnv_met == 'cnv' & pval >= input$pvalRangeVenn[1] &
                                 pval <= input$pvalRangeVenn[2],
                               select = c('cov'))
      met_sign_genes <- subset(data_table,
                               cnv_met == 'met' & pval >= input$pvalRangeVenn[1] &
                                 pval <= input$pvalRangeVenn[2],
                               select = c('cov'))
    }else{
      cnv_sign_genes <- subset(data_table,
                               cnv_met == 'cnv' & fdr >= input$FDRRangeVenn[1] &
                                 fdr <= input$FDRRangeVenn[2],
                               select = c('cov'))
      met_sign_genes <- subset(data_table,
                               cnv_met == 'met' & fdr >= input$FDRRangeVenn[1] &
                                 fdr <= input$FDRRangeVenn[2],
                               select = c('cov'))
    }

    data_venn <- list(cnv_sign_genes = cnv_sign_genes,
                met_sign_genes = met_sign_genes)
    return(data_venn)
    })%>%bindEvent(input$classSelectVenn,
                   input$FDRRangeVenn,
                   input$pvalRangeVenn,
                   input$significativityCriteriaVenn)
}
#####################################################################
######################################################################

.prepare_reactive_volcano <- function(data_table,
                                      input,
                                      output,
                                      type = "genomic",
                                      deg = FALSE){
  reactive({


    if(!"class"%in%colnames(data_table)) {data_table$class <- data_table$omics}
    if(type=="genomic"){

      if("gene_genomic_res"%in%unique(data_table$omics)){
        data_volcano <- data_table[data_table$omics == 'gene_genomic_res',]
        data_volcano <- data_volcano[data_volcano$cnv_met == input$genomicTypeSelectVolcano,]
      }
      if("cnv_gene_res"%in%unique(data_table$omics)){
        data_volcano <- data_table[data_table$omics == 'cnv_gene_res',]
      }
      if("met_gene_res"%in%unique(data_table$omics)){
        data_volcano <- data_table[data_table$omics == 'met_gene_res',]
      }
      if(deg==TRUE){data_volcano <- data_volcano[data_volcano$deg,]}
    }

    if(input$significativityCriteriaVolcano == 'pval'){
      data_volcano["group"] <- "Not Significant"
      data_volcano[data_volcano$pval <= input$pvalRangeVolcano, 'group'] <- "Significant"
      top_peaks <- data_volcano[with(data_volcano,
                                     order(pval, coef)),]
      data_volcano <- rbind(top_peaks, data_volcano[with(data_volcano,
                                                         order(-coef, pval)), ][1:10,])

      data_volcano$pval_fdr <- -log10(data_volcano$pval)
    }else{
      data_volcano[data_volcano$fdr <= input$FDRRangeVolcano, 'group'] <- "Significant"
      top_peaks <- data_volcano[with(data_volcano,
                                     order(fdr, coef)),]
      data_volcano <- rbind(top_peaks, data_volcano[with(data_volcano,
                                                         order(-coef, fdr)),][1:10,])
      data_volcano$pval_fdr <- -log10(data_volcano$fdr)
    }
    return(data_volcano)
  })%>%bindEvent(input$genomicTypeSelectVolcano,
                 input$FDRRangeVolcano,
                 input$pvalRangeVolcano,
                 input$significativityCriteriaVolcano)
}
#######################################################################
#######################################################################

.prepare_reactive_heatmap <- function(data_table,
                             multiomics_integration,
                             input,
                             output,
                             session){
observe({

  data_table <- data_table[data_table$cov != '(Intercept)',]

  if ("class" %in% colnames(data_table)) {
    df_heatmap <- multiomics_integration[[input$integrationSelectHeatmap]][[input$selectClassHeatmap]]$data$response_var
    data_table <- data_table[data_table$class == input$selectClassHeatmap,]
  }

  if ("deg" %in% colnames(data_table) & input$degSelectHeatmap == "Only DEGs") {
    df_heatmap <- multiomics_integration[[input$integrationSelectHeatmap]][[input$selectClassHeatmap]]$data$response_var
    data_table <- data_table[data_table$deg,]
    df_heatmap <- df_heatmap[, unique(data_table$response)]
  }

  if (!"class" %in% colnames(data_table)) {
    df_heatmap <- multiomics_integration[[input$integrationSelectHeatmap]]$data$response_var
    data_table <- data_table[data_table$omics == input$integrationSelectHeatmap,]
  }

  df_heatmap_t <- t(as.matrix(df_heatmap))
  if (input$integrationSelectHeatmap == "gene_genomic_res") {
    data_table <- data_table[data_table$omics == 'gene_genomic_res',]
    tmp <- data.frame(cnv = data_table$coef[data_table$cnv_met == 'cnv'],
                      pval_cnv = data_table$pval[data_table$cnv_met == 'cnv'],
                      fdr_cnv = data_table$fdr[data_table$cnv_met == 'cnv'],
                      row.names = data_table$response[data_table$cnv_met == 'cnv'])

    tmp2 <- data.frame(met = data_table$coef[data_table$cnv_met == 'met'],
                       pval_met = data_table$pval[data_table$cnv_met == 'met'],
                       fdr_met = data_table$fdr[data_table$cnv_met == 'met'],
                       row.names = data_table$response[data_table$cnv_met == 'met'])

    tmp3 <- unique(c(rownames(tmp), rownames(tmp2)))

    data_table <- cbind(cnv = tmp[tmp3,], met = tmp2[tmp3,])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp3
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])

    if (input$significativityCriteriaHeatmap == 'pval') {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= input$pvalRangeHeatmap | df_heatmap_t$pval_met <= input$pvalRangeHeatmap,]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= input$FDRRangeHeatmap | df_heatmap_t$fdr_met <= input$FDRRangeHeatmap,]
    }

    if (nrow(df_heatmap_t) == 0) {
      return(NULL)
    }

    top_met <- df_heatmap_t %>% arrange(desc(abs(met))) %>% head(input$numTopGenesHeatmapCNV)
    top_cnv <- df_heatmap_t %>% arrange(desc(abs(cnv))) %>% head(input$numTopGenesHeatmapMET)
    expr_top <- rbind(top_cnv, top_met[!rownames(top_met) %in% rownames(top_cnv),])
    expr_top_subset <- expr_top[, -c((ncol(expr_top) - 5):ncol(expr_top))]

    set.seed(123)
    row_ha <- rowAnnotation(coef_cnv = expr_top$cnv, coef_met = expr_top$met)
    ht <- ComplexHeatmap:::Heatmap(expr_top_subset, right_annotation = row_ha)
    ht <- draw(ht)
    ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')
  }

  if(input$integrationSelectHeatmap == "gene_cnv_res"){

    data_table <- data_table[data_table$omics == 'gene_cnv_res',]
    tmp <-  data.frame(cnv=data_table$coef[data_table$omics == 'cnv'],
                       pval_cnv=data_table$pval[data_table$omics == 'cnv'],
                       fdr_cnv=data_table$fdr[data_table$omics == 'cnv'],
                       row.names = data_table$response[data_table$omics == 'cnv'])
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(cnv=tmp[tmp2,])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
    if (input$significativityCriteriaHeatmap == 'pval'){
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= input$pvalRangeHeatmap,]
    }else{
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= input$FDRRangeHeatmap,]
    }

    top_cnv <- df_heatmap_t %>%
      arrange(desc(abs(cnv))) %>%
      head(input$numTopGenesHeatmapCNVonly)
      expr_top <- top_cnv

    set.seed(123)
    row_ha <- rowAnnotation(coef_cnv=expr_top$cnv)
    ht <- Heatmap(as.matrix(expr_top, right_annotation=row_ha))
    ht = draw(ht)
    ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')

  }

  if(input$integrationSelectHeatmap == "gene_met_res"){

    data_table <- data_table[data_table$omics == 'gene_met_res',]
    tmp <-  data.frame(met=data_table$coef[data_table$omics == 'met'],
                       pval_met=data_table$pval[data_table$omics == 'met'],
                       fdr_met=data_table$fdr[data_table$omics == 'met'],
                       row.names = data_table$response[data_table$omics == 'met'])
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(met=tmp[tmp2,])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
    if (input$significativityCriteriaHeatmap == 'pval'){
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_met <= input$pvalRangeHeatmap,]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_met <= input$FDRRangeHeatmap,]
    }

    top_met <- df_heatmap_t %>%
      arrange(desc(abs(met))) %>%
      head(input$numTopGenesHeatmapMETonly)
      expr_top <- top_met

    set.seed(123)
    row_ha <- rowAnnotation(coef_met=expr_top$met)
    ht <- Heatmap(as.matrix(expr_top, right_annotation=row_ha))
    ht = draw(ht)
    ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')

  }
  if(input$integrationSelectHeatmap == "mirna_cnv_res"){


    data_table <- data_table[data_table$omics == 'mirna_cnv_res',]
    tmp <-  data.frame(mirna_cnv=data_table$coef[data_table$omics == 'mirna_cnv_res'],
                       pval_mirna_cnv=data_table$pval[data_table$omics == 'mirna_cnv_res'],
                       fdr_mirna_cnv=data_table$fdr[data_table$omics == 'mirna_cnv_res'],
                       row.names = data_table$response[data_table$omics == 'mirna_cnv_res'])
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(mirna_cnv=tmp[tmp2,])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
    if (input$significativityCriteriaHeatmap == 'pval'){
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_mirna_cnv <= input$pvalRangeHeatmap,]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_mirna_cnv <= input$FDRRangeHeatmap,]
    }

    df_heatmap_t <- as.data.frame(df_heatmap_t)
    top_mirna_cnv <- df_heatmap_t %>%
    arrange(desc(abs(mirna_cnv)))  %>%
    head(input$numTopGenesHeatmapmirna_cnv)
    expr_top <- top_mirna_cnv
    expr_top_subset <- expr_top[, -c((ncol(expr_top) - 3):ncol(expr_top))]

    set.seed(123)
    row_ha <- rowAnnotation(coef_mirna=expr_top$mirna_cnv)
    ht <- ComplexHeatmap:::Heatmap(expr_top_subset, right_annotation=row_ha)
    ht = draw(ht)
    ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')

  }
  })%>%bindEvent(input$integrationSelectHeatmap,
                 input$numTopGenesHeatmapCNV,
                 input$numTopGenesHeatmapCNVonly,
                 input$numTopGenesHeatmapMET,
                 input$numTopGenesHeatmapMETonly,
                 input$numTopGenesHeatmapmirna_cnv,
                 input$FDRRangeHeatmap,
                 input$pvalRangeHeatmap,
                 input$significativityCriteriaHeatmap,
                 input$selectClassHeatmap,
                 input$degSelectHeatmap
                 )
}

######################################################################
######################################################################

.prepare_reactive_ridge <- function(data_table,
                                    input,
                                    output,
                                    type = "genomic",
                                    deg = FALSE){

  reactive({
    df <- data_table
    if ('class' %in% colnames(df)) {
      df <- df[df$class == input$classSelectRidge,]
    }
    if(type == "genomic"){
      if("gene_genomic_res"%in%unique(df$omics)){
      df <- df[df$omics == "gene_genomic_res",]
      df <- df[df$cnv_met == input$genomicTypeSelectRidge,]
      df <- df[!is.na(df$cnv_met),]
      }
      if("cnv_gene_res"%in%unique(df$omics)){
        df <- df[df$omics == "cnv_gene_res",]
      }
      if("met_gene_res"%in%unique(df$omics)){
        df <- df[df$omics == "met_gene_res",]
      }

      if(deg == TRUE){df <- df[df$deg,]}
    }

    if (input$significativityCriteriaRidge == 'pval') {
      df$significance <- ifelse(df$pval >= input$pvalRangeRidge[1] &
                                  df$pval <= input$pvalRangeRidge[2],
                                "Significant",
                                "Not Significant")
    } else {
      df$significance <- ifelse(df$fdr >= input$FDRRangeRidge[1] &
                                  df$fdr <= input$FDRRangeRidge[2],
                                "Significant",
                                "Not Significant")
    }
    lower_quantile <- quantile(df$coef,
                               0.001)
    upper_quantile <- quantile(df$coef,
                               0.999)
    ans <- list(df=df, quantiles=c(lower_quantile, upper_quantile))
    return(ans)
  })%>%bindEvent(input$classSelectRidge,
                 input$pvalRangeRidge,
                 input$FDRRangeRidge,
                 input$significativityCriteriaRidge,
                 input$genomicTypeSelectRidge)
}

#######################################################################
########################################################################

.prepare_reactive_ridge_table <- function(data_table,
                                    input,
                                    output,
                                    type = "genomic",
                                    deg = FALSE){

  reactive({
    df <- data_table
    if ('class' %in% colnames(data_table)) {
      df <- df[df$class == input$classSelectRidge,]
    }
    if(type=="genomic"){
      if("gene_genomic_res"%in%unique(df$omics)){
        df <- df[df$omics == "gene_genomic_res", ]
        df <- df[df$cnv_met == input$genomicTypeSelectRidge,]
        df <- df[!is.na(df$cnv_met),]
      }
      if("cnv_gene_res"%in%unique(df$omics)){
        df <- df[df$omics == "cnv_gene_res", ]
      }
      if("met_gene_res"%in%unique(df$omics)){
        df <- df[df$omics == "met_gene_res", ]
      }
      if(deg==TRUE){df <- df[df$deg,]}
    }

    if (input$significativityCriteriaRidge == 'pval') {
      df <- df[df$pval >= input$pvalRangeRidge[1] &
                                  df$pval <= input$pvalRangeRidge[2],
      ]
    } else {
      df <- df[df$fdr >= input$FDRRangeRidge[1] &
                                  df$fdr <= input$FDRRangeRidge[2],
      ]
    }

    return(df)
  })%>%bindEvent(input$classSelectRidge,
                 input$pvalRangeRidge,
                 input$FDRRangeRidge,
                 input$significativityCriteriaRidge,
                 input$genomicTypeSelectRidge)
}

#######################################################################
########################################################################

.prepare_reactive_histo <- function(data_table,
                                    input,
                                    output){
  reactive({
    data_table <- data_table[!is.na(data_table$chr_cov),]
    chr_order <- gtools::mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    df_filtered_histo <- data_table[data_table$omics == input$integrationSelectHisto, ]
    if(input$chrSelectHisto != "All"){
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$chr_cov == input$chrSelectHisto, ]
    }
    if('class' %in% colnames(df_filtered_histo)){df_filtered_histo <- df_filtered_histo[df_filtered_histo$class == input$classSelectHisto, ]
    }
    if(input$significativityCriteriaHisto == 'pval'){
      df_filtered_histo$significance <- ifelse(df_filtered_histo$pval >= input$pvalRangeHisto[1] &
                                                 df_filtered_histo$pval <= input$pvalRangeHisto[2], "Significant",
                                               "Not Significant")
    }else{
      df_filtered_histo$significance <- ifelse(df_filtered_histo$fdr >= input$FDRRangeHisto[1] &
                                                 df_filtered_histo$fdr <= input$FDRRangeHisto[2], "Significant",
                                               "Not Significant")
    }
    if(input$degSelectHisto == 'Only DEGs' & 'class' %in% names(df_filtered_histo)){
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$deg,]
    }
    return(df_filtered_histo)
  })%>%bindEvent(input$classSelectHisto,
                 input$integrationSelectHisto,
                 input$chrSelectHisto,
                 input$pvalRangeHisto,
                 input$FDRRangeHisto,
                 input$degSelectHisto,
                 input$significativityCriteriaHisto)
}

#######################################################################
########################################################################
.prepare_reactive_histo_table <- function(data_table,
                                          input,
                                          output){
  reactive({
    data_table <- data_table[!is.na(data_table$chr_cov),]
    chr_order <- gtools::mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    df_filtered_histo <- data_table[data_table$omics == input$integrationSelectHisto,]
    if(input$chrSelectHisto != "All"){
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$chr_cov == input$chrSelectHisto,]
    }
    if('class' %in% colnames(df_filtered_histo)){df_filtered_histo <- df_filtered_histo[df_filtered_histo$class == input$classSelectHisto,]
    }
    if(input$significativityCriteriaHisto == 'pval'){
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$pval >= input$pvalRangeHisto[1] &
                                                 df_filtered_histo$pval <= input$pvalRangeHisto[2],]
    }else{
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$fdr >= input$FDRRangeHisto[1] &
                                                 df_filtered_histo$fdr <= input$FDRRangeHisto[2],]
    }
    if(input$degSelectHisto == 'Only DEGs' & 'class' %in% colnames(df_filtered_histo)){
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$deg,]
    }
    return(df_filtered_histo)
  })%>%bindEvent(input$classSelectHisto,
                 input$integrationSelectHisto,
                 input$chrSelectHisto,
                 input$pvalRangeHisto,
                 input$FDRRangeHisto,
                 input$degSelectHisto,
                 input$significativityCriteriaHisto)
}
#######################################################################
########################################################################

.prepare_reactive_histo_tf <- function(data_table,
                                    input,
                                    output){
  reactive({
    chr_order <- gtools::mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    data_table <- data_table[data_table$omics == 'tf_res', ]
    data_table <- data_table[,
                             colnames(data_table)%in%c("response", "cov",
                                                       "pval","fdr","chr_cov",
                                                       "deg","class")]
    df_filtered_histo_tf <- data_table
    if('class' %in% names(df_filtered_histo_tf)){
      df_filtered_histo_tf <- df_filtered_histo_tf[
        df_filtered_histo_tf$class == input$classSelectHistoTFs, ]
    }
    if(input$degSelectHistoTFs == 'Only DEGs'){
      df_filtered_histo_tf <- df_filtered_histo_tf[df_filtered_histo_tf$deg, ]
    }
    df_filtered_histo_tf <- df_filtered_histo_tf[df_filtered_histo_tf$pval <= 0.05, ]
    genes_count <- table(df_filtered_histo_tf$cov, df_filtered_histo_tf$chr_cov)
    genes_count_df <- as.data.frame.table(genes_count)
    genes_count_df <- subset(genes_count_df, Freq != 0)
    colnames(genes_count_df) <- c("TF", "Chromosome", "Count")
    genes_count_df <- genes_count_df[order(-genes_count_df$Count), ]
    return(genes_count_df)
  })%>%bindEvent(input$classSelectHistoTFs,
                 input$degSelectHistoTFs)
}

#######################################################################
########################################################################

.prepare_reactive_table <- function(data_table,
                                    input,
                                    output){
    reactive({
      filtered_df <- data_table[!is.na(data_table$chr_cov),]
      filtered_df <- filtered_df[filtered_df$omics == input$integrationSelectTable,]

      if('class' %in% colnames(filtered_df)){
        filtered_df <- filtered_df[filtered_df$class == input$classSelectTable,]}

      if(input$significativityCriteriaTable == 'pval'){
        filtered_df <- filtered_df[filtered_df$pval >= input$pvalRangeTable[1] &
                                   filtered_df$pval <= input$pvalRangeTable[2],]
      }else{
        filtered_df <- filtered_df[filtered_df$fdr >= input$FDRRangeTable[1] &
                                   filtered_df$fdr <= input$FDRRangeTable[2],]
      }
      if(input$degSelectTable == 'Only DEGs'){
        filtered_df <- filtered_df[filtered_df$deg,]}

      filtered_df <- filtered_df[filtered_df$chr_cov == input$chrSelectTable,]
      return(filtered_df)
  })%>%bindEvent(input$integrationSelectTable,
                 input$classSelectTable,
                 input$chrSelectTable,
                 input$pvalRangeTable,
                 input$FDRRangeTable,
                 input$degSelectTable,
                 input$significativityCriteriaTable)
}

#######################################################################
########################################################################

.background_srv <- function(input,
                            output,
                            session,
                            data_gen_enrich,
                            data_tf_enrich){

  gen_enr <- gINTomics:::.run_bg(FFUN = run_genomic_enrich,
                                  input = input,
                                  output = output,
                                  args = list(model_results = NULL,
                                              qvalueCutoff = 1,
                                              pvalueCutoff = 1,
                                              extracted_data = data_gen_enrich))
  check <- gINTomics:::.check_reactive_bg_enrich(bg_enrich = gen_enr,
                                                 input = input,
                                                 output = output,
                                                 session = session)

  output$gen_enrichment <- renderText({
    check()
  })

  gen_plot <- gINTomics:::.reactive_gen_dotplot(bg_enrich = gen_enr,
                                                 input = input,
                                                 output = output,
                                                 session = session)
  gINTomics:::.render_gen_dotplot(gen_plot=gen_plot,
                                  data_gen_enrich = data_gen_enrich,
                                  input = input,
                                  output = output,
                                  session = session)


  tf_enr <- gINTomics:::.run_bg(FFUN = run_tf_enrich,
                                 input = input,
                                 output = output,
                                 args = list(model_results = NULL,
                                             qvalueCutoff = 1,
                                             pvalueCutoff = 1,
                                             extracted_data = data_tf_enrich))
  check_tf <- gINTomics:::.check_reactive_bg_enrich(bg_enrich = tf_enr,
                                                 input = input,
                                                 output = output,
                                                 session = session)
  output$tf_enrichment <- renderText({
    check_tf()
  })
}


#######################################################################
########################################################################

.run_bg <- function(FFUN,
                    args,
                    input,
                    output){
  #reactive({
    ans <- callr::r_bg(func = FFUN,
                           args = args,
                           supervise = T)
    return(ans)
  #})
}

#######################################################################
########################################################################

.check_reactive_bg_enrich <- function(bg_enrich,
                                        input,
                                        output,
                                        session){
  reactive({
    invalidateLater(millis = 1000, session = session)
    if (bg_enrich$is_alive()) {
      x <- "Enrichment running in background, this may take several minutes"
    } else {
      x <- "Enrichment completed"
    }
    return(x)
  })
}

#######################################################################
########################################################################

.reactive_gen_dotplot <- function(bg_enrich,
                                  input,
                                  output,
                                  session){
  reactive({
    if (bg_enrich$is_alive()){
      invalidateLater(millis = 1000, session = session)
      }
    if (!bg_enrich$is_alive()) {
      data <- bg_enrich$get_result()
      if(sum(c("cnv", "met")%in%names(data))){
      cnv <- data[["cnv"]][[1]][[1]]
      cnv <- clusterProfiler::dotplot(cnv)
      cnv <- ggplotly(cnv)
      met <- data[["met"]][[1]][[1]]
      met <- clusterProfiler::dotplot(met)
      met <- ggplotly(met, width = 800, height = 700)
      ans <- list(cnv=cnv, met=met)
      }else{
        ans <- data[[1]][[1]]
        ans <- clusterProfiler::dotplot(ans) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))

        cnv <- ggplotly(ans, width = 800 ,height = 700)
      }
      return(ans)
      }
  })
}

#######################################################################
########################################################################

.render_gen_dotplot <- function(gen_plot,
                                data_gen_enrich,
                                input,
                                output,
                                session){
  observe({
    if("gene_genomic_res"%in%data_gen_enrich$omics){
      output$gen_dotplot <- renderPlotly({
        ans=gen_plot()
        ans=ans[[input$genomicTypeSelectEnrich]]
      })
    }else{
      output$gen_dotplot <- renderPlotly({
        gen_plot()
      })
    }
  })%>%bindEvent(input$genomicTypeSelectEnrich)
}

#######################################################################
########################################################################

.prepare_reactive_circos <- function(data, input, output) {
  reactive({
  if("class"%in%colnames(data)) data <- data[data$class==input$circosClass,]
  gr <- .circos_preprocess(data = data)
  tracks <- .create_tracks(data = data, gr = gr)
  width=800
  height=800
  if(input$circosLayout=="linear"){
    width=800
    height=100
  }
  composed_view <- .create_composed_view(tracks, height=height, width=width)
  if(input$circosType=="Gene"){
    ssel <- intersect(c("circos_genomic", "circos_met_gene", "circos_met_gene"),
                      names(composed_view))
  }
  if(input$circosType=="miRNA"){
    ssel <- intersect(c("circos_genomic_mirna"), names(composed_view))
  }
  composed_view <- composed_view[[ssel]]
  arranged_view <- arrange_views(title = 'Interactive Circos',
                                  subtitle = 'subtitle',
                                  views = composed_view,
                                  layout = input$circosLayout,
                                 xDomain = list(
                                    chromosome = input$circosChr
                                  ), assembly="hg38")

  return(arranged_view)
  })%>%bindEvent(input$circosClass,
                 input$circosLayout,
                 input$circosType,
                 input$circosChr)
}


