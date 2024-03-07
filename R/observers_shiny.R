#' @importFrom visNetwork renderVisNetwork
#' @importFrom visNetwork visHierarchicalLayout
#' @importFrom visNetwork  visIgraphLayout

.render_reactive_network <- function(reactive_network,
                                     input,
                                     output,
                                     deg = FALSE){
  observe({
    if(!deg){

        output$networkPlot <- renderVisNetwork({
        network_data <- reactive_network()
        network <- .build_network(nodes=network_data$nodes,
                                  edges=network_data$edges,
                                  legend_nodes=network_data$legend_nodes,
                                  legend_edges=network_data$legend_edges)
        network <- network%>%
          visHierarchicalLayout(enabled = input$layoutNetwork)%>%
          visIgraphLayout(physics = input$physics,
                          layout = "layout_with_fr",
                          randomSeed = 20)

        return(network)
      })
    } else {
      output$networkPlotDEG <- renderVisNetwork({
        network_data <- reactive_network()
        network <- .build_network(nodes=network_data$nodes,
                                  edges=network_data$edges,
                                  legend_nodes=network_data$legend_nodes,
                                  legend_edges=network_data$legend_edges)
        network <- network%>%
          visHierarchicalLayout(enabled = input$layoutNetworkDEG)%>%
          visIgraphLayout(physics = input$physicsDEG,
                          layout = "layout_with_fr",
                          randomSeed = 20)

        return(network)
      })
    }
  })%>%bindEvent(input$layoutNetwork,
                 input$physics)
}

################################################################
#################################################################

.select_network <- function(data_table,
                            network_data,
                            input,
                            output,
                            deg = FALSE){
  reactive({
    data_table <- data_table[data_table$omics%in%c("tf_res",
                                                   "tf_mirna_res",
                                                   "mirna_target_res"),]

    if(!deg){
      numInteractions <- input$numInteractions
      significativityCriteria <- input$significativityCriteriaNetwork
      pval <- input$pvalNetwork
      fdr <- input$fdrNetwork
      class <- input$classSelectNetwork
    }else{
      numInteractions <- input$numInteractionsDEG
      significativityCriteria <- input$significativityCriteriaNetworkDEG
      pval <- input$pvalNetworkDEG
      fdr <- input$fdrNetworkDEG
      class <- input$classSelectNetworkDEG
    }
      ans <- network_data
      ans$edges <- ans$edges[ans$edges$class==class,]
      nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
      ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges,]
      if(significativityCriteria == "pval"){
        ans$edges <- ans$edges[ans$edges$pval<= pval,]
      } else {
        ans$edges <- ans$edges[ans$edges$fdr<= fdr,]
      }
      if(deg){
        tto <- unique(data_table$response[data_table$deg])
        ans$edges <- ans$edges[ans$edges$to%in%tto,]
        }
      ans$edges <- ans$edges[seq_len(length.out = numInteractions),]
      nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
      ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges,]
      return(ans)
  })%>%bindEvent(input$layoutNetwork,
                 input$numInteractions,
                 input$significativityCriteriaNetwork,
                 input$pvalNetwork,
                 input$fdrNetwork,
                 input$physics,
                 input$classSelectNetwork,
                 input$layoutNetworkDEG,
                 input$numInteractionsDEG,
                 input$significativityCriteriaNetworkDEG,
                 input$pvalNetworkDEG,
                 input$fdrNetworkDEG,
                 input$physicsDEG,
                 input$classSelectNetworkDEG)
}

################################################################
#################################################################
.prepare_reactive_venn <- function(data_table,
                                   input,
                                   output,
                                   deg = FALSE){
  reactive_venn <- reactive({
    data_venn <- NULL
    if(deg==FALSE & "gene_genomic_res"%in%unique(data_table$omics)) {
      classSelect <- input$classSelectVenn
      pvalRange <- input$pvalRangeVenn
      fdrRange <- input$fdrRangeVenn
      significativityCriteria <- input$significativityCriteriaVenn
      data_table <- data_table[data_table$omics == "gene_genomic_res", ]
      if("class"%in%colnames(data_table)){
        data_table <- data_table[data_table$class == classSelect, ]
      }
      if(significativityCriteria == 'pval'){
        cnv_sign_genes <- subset(data_table,
                                 cnv_met == 'cnv' & pval >= pvalRange[1] &
                                   pval <= pvalRange[2],
                                 select = c('cov'))
        met_sign_genes <- subset(data_table,
                                 cnv_met == 'met' & pval >= pvalRange[1] &
                                   pval <= pvalRange[2],
                                 select = c('cov'))
      }else{
        cnv_sign_genes <- subset(data_table,
                                 cnv_met == 'cnv' & fdr >= fdrRange[1] &
                                   fdr <= fdrRange[2],
                                 select = c('cov'))
        met_sign_genes <- subset(data_table,
                                 cnv_met == 'met' & fdr >= fdrRange[1] &
                                   fdr <= fdrRange[2],
                                 select = c('cov'))
      }
      data_venn <- list(cnv_sign_genes = cnv_sign_genes,
                        met_sign_genes = met_sign_genes)
    }
    if(deg==FALSE & !"gene_genomic_res"%in%unique(data_table$omics)) {
      if("gene_cnv_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="gene_cnv_res",]
      }
      if("gene_met_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="gene_met_res",]
      }
    }
    if(deg==TRUE & "gene_genomic_res"%in%unique(data_table$omics)){
      classSelect <- input$classSelectVennDEG
      significativityCriteria <- input$significativityCriteriaVennDEG
      pvalRange <- input$pvalRangeVennDEG
      fdrRange <- input$fdrRangeVennDEG
      if('class'%in%colnames(data_table)){
        data_table <- data_table[data_table$class == classSelect, ]
        data_table <- data_table[data_table$deg,]
        if(significativityCriteria == 'pval'){
          cnv_sign_genes <- subset(data_table,
                                   cnv_met == 'cnv' & pval >= pvalRange[1] &
                                     pval <= pvalRange[2],
                                   select = c('cov'))
          met_sign_genes <- subset(data_table,
                                   cnv_met == 'met' & pval >= pvalRange[1] &
                                     pval <= pvalRange[2],
                                   select = c('cov'))
        }else{
          cnv_sign_genes <- subset(data_table,
                                   cnv_met == 'cnv' & fdr >= fdrRange[1] &
                                     fdr <= fdrRange[2],
                                   select = c('cov'))
          met_sign_genes <- subset(data_table,
                                   cnv_met == 'met' & fdr >= fdrRange[1] &
                                     fdr <= fdrRange[2],
                                   select = c('cov'))
        }
        data_venn <- list(cnv_sign_genes = cnv_sign_genes,
                          met_sign_genes = met_sign_genes)
      }else(return(NULL))
    }
    if(deg==TRUE & !"gene_genomic_res"%in%unique(data_table$omics)){
      if("gene_cnv_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="gene_cnv_res",]
      }
      if("gene_met_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="gene_met_res",]
      }
    }
    return(data_venn)
  })%>%bindEvent(input$classSelectVenn,
                 input$fdrRangeVenn,
                 input$pvalRangeVenn,
                 input$significativityCriteriaVenn,
                 input$classSelectVennDEG,
                 input$fdrRangeVennDEG,
                 input$pvalRangeVennDEG,
                 input$significativityCriteriaVennDEG)
}
#####################################################################
######################################################################

.prepare_reactive_volcano <- function(data_table,
                                      input,
                                      output,
                                      type = "genomic",
                                      deg = FALSE){
  reactive({
    if(type == "genomic"){
      integrationSelect <- input$genomicIntegrationSelectVolcano
      typeSelect <- input$genomicTypeSelectVolcano
      significativityCriteria <- input$genomicSignificativityCriteriaVolcano
      pvalRange <- input$genomicPvalRangeVolcano
      fdrRange <- input$genomicFdrRangeVolcano
    }
    if(type == "transcript") {
      integrationSelect <- input$transcriptIntegrationSelectVolcano
      significativityCriteria <- input$transcriptSignificativityCriteriaVolcano
      pvalRange <- input$transcriptPvalRangeVolcano
      fdrRange <- input$transcriptFdrRangeVolcano
    }
    if(type == "all" & deg == TRUE){
      integrationSelect <- input$integrationSelectVolcanoDEG
      typeSelect <- input$typeSelectVolcanoDEG
      significativityCriteria <- input$significativityCriteriaVolcanoDEG
      pvalRange <- input$pvalRangeVolcanoDEG
      fdrRange <- input$fdrRangeVolcanoDEG
      data_table <- data_table[data_table$deg,]
    }
    if(!"class"%in%colnames(data_table)) {data_table$class <- data_table$omics}
    data_table <- data_table[data_table$omics == integrationSelect,]
    if(integrationSelect == "gene_genomic_res"){
      data_table <- data_table[data_table$cnv_met == typeSelect,]
    }
    if(significativityCriteria == 'pval'){
      data_table["group"] <- "Not Significant"
      data_table[data_table$pval <= pvalRange, 'group'] <- "Significant"
      data_table$pval_fdr <- -log10(data_table$pval)
    }else{
      data_table["group"] <- "Not Significant"
      data_table[data_table$fdr <= fdrRange, 'group'] <- "Significant"
      data_table$pval_fdr <- -log10(data_table$fdr)
    }
    return(data_table)
  })%>%bindEvent(input$genomicIntegrationSelectVolcano,
                 input$genomicTypeSelectVolcano,
                 input$genomicSignificativityCriteriaVolcano,
                 input$genomicPvalRangeVolcano,
                 input$genomicFdrRangeVolcano,
                 input$transcriptIntegrationSelectVolcano,
                 input$transcriptSignificativityCriteriaVolcano,
                 input$transcriptPvalRangeVolcano,
                 input$transcriptFdrRangeVolcano,
                 input$integrationSelectVolcanoDEG,
                 input$significativityCriteriaVolcanoDEG,
                 input$pvalRangeVolcanoDEG,
                 input$fdrRangeVolcanoDEG,
                 input$typeSelectVolcanoDEG)
}
#######################################################################
#######################################################################
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap rowAnnotation

.prepare_reactive_heatmap <- function(data_table,
                                      multiomics_integration,
                                      input,
                                      output,
                                      session,
                                      deg = FALSE){
  observe({
    ht_opt$message = FALSE
    if(deg == FALSE){
      integrationSelect <- input$integrationSelectHeatmap
      numTopCNV <- input$numTopGenesHeatmapCNV
      numTopMET <- input$numTopGenesHeatmapMET
      numTopCNVonly <- input$numTopGenesHeatmapCNVonly
      numTopMETonly <- input$numTopGenesHeatmapMETonly
      numTopMiCNV <- input$numTopGenesHeatmapmirna_cnv
      classSelect <- input$classSelectHeatmap
      significativityCriteria <- input$significativityCriteriaHeatmap
      pvalRange <- input$pvalRangeHeatmap
      fdrRange <- input$FDRRangeHeatmap
      scale <- input$scaleHeatmap
    }
    if(deg == TRUE){
      integrationSelect <- input$integrationSelectHeatmapDEG
      numTopCNV <- input$numTopGenesHeatmapCNVDEG
      numTopMET <- input$numTopGenesHeatmapMETDEG
      numTopCNVonly <- input$numTopGenesHeatmapCNVonlyDEG
      numTopMETonly <- input$numTopGenesHeatmapMETonlyDEG
      numTopMiCNV <- input$numTopGenesHeatmapmirna_cnvDEG
      significativityCriteria <- input$significativityCriteriaHeatmapDEG
      pvalRange <- input$pvalRangeHeatmapDEG
      fdrRange <- input$FDRRangeHeatmapDEG
      classSelect <- input$classSelectHeatmapDEG
      scale <- input$scaleHeatmapDEG
    }
    if("class"%in%colnames(data_table) & deg == FALSE){
      df_heatmap <- multiomics_integration[[integrationSelect]][[
        classSelect]]$data$response_var
      data_table <- data_table[data_table$omics == integrationSelect,]
      data_table <- data_table[data_table$class == classSelect,]
    }
    if("class"%in%colnames(data_table) & deg == TRUE){
      df_heatmap <- multiomics_integration[[integrationSelect]][[
        classSelect]]$data$response_var
      data_table <- data_table[data_table$omics == integrationSelect,]
      data_table <- data_table[data_table$class == classSelect,]
      data_table <- data_table[data_table$deg,]
      df_heatmap <- df_heatmap[, unique(data_table$response)]
    }
    if(!"class"%in%colnames(data_table)){
      df_heatmap <- multiomics_integration[[
        integrationSelect]]$data$response_var
      data_table <- data_table[data_table$omics == integrationSelect,]
    }
    df_heatmap_t <- t(as.matrix(df_heatmap))
    if (integrationSelect == "gene_genomic_res"){
      ans <- .prepare_gen_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria=significativityCriteria,
                              pvalRange = pvalRange,
                              fdrRange = fdrRange,
                              numTopCNV = numTopCNV,
                              numTopMET = numTopMET,
                              scale = scale)
    }
    if(integrationSelect == "gene_cnv_res"){
      ans <- .prepare_cnv_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria=significativityCriteria,
                              pvalRange = pvalRange,
                              fdrRange = fdrRange,
                              numTopCNVonly = numTopCNVonly,
                              scale = scale)
    }
    if(integrationSelect == "gene_met_res"){
      ans <- .prepare_met_heatmap(data_table = data_table,
                               df_heatmap = df_heatmap,
                               df_heatmap_t = df_heatmap_t,
                               significativityCriteria=significativityCriteria,
                               pvalRange = pvalRange,
                               fdrRange = fdrRange,
                               numTopMETonly = numTopMETonly,
                               scale = scale)
    }
    if(integrationSelect == "mirna_cnv_res"){
      ans <- .prepare_mirna_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria=significativityCriteria,
                              pvalRange = pvalRange,
                              fdrRange = fdrRange,
                              numTopMiCNV = numTopMiCNV,
                              scale = scale)
    }
    if(is.null(ans)) return(NULL)
    if(deg == FALSE){
      ht <- makeInteractiveComplexHeatmap(input,output,session,ans,'heatmap')
    }else{
      ht <- makeInteractiveComplexHeatmap(input,output,session,ans,'heatmapDEG')
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
                 input$classSelectHeatmap,
                 input$integrationSelectHeatmapDEG,
                 input$numTopGenesHeatmapCNVDEG,
                 input$numTopGenesHeatmapMETDEG,
                 input$numTopGenesHeatmapCNVonlyDEG,
                 input$numTopGenesHeatmapMETonlyDEG,
                 input$numTopGenesHeatmapmirna_cnvDEG,
                 input$significativityCriteriaHeatmapDEG,
                 input$pvalRangeHeatmapDEG,
                 input$FDRRangeHeatmapDEG,
                 input$classSelectHeatmapDEG,
                 input$scaleHeatmap,
                 input$scaleHeatmapDEG)
}

######################################################################
######################################################################

.prepare_reactive_ridge <- function(data_table,
                                    input,
                                    output,
                                    type = "genomic",
                                    deg = FALSE) {
  reactive({
    df <- data_table
    if(type == "genomic"){
      integrationSelect <- input$genomicIntegrationSelectRidge
      classSelect <- input$genomicClassSelectRidge
      significativityCriteria <- input$genomicSignificativityCriteriaRidge
      pvalRange <- input$genomicPvalRangeRidge
      fdrRange <- input$genomicFdrRangeRidge
      typeSelect <- input$genomicTypeSelectRidge
    }
    if(type == "transcript"){
      integrationSelect <- input$transcriptIntegrationSelectRidge
      classSelect <- input$transcriptClassSelectRidge
      significativityCriteria <- input$transcriptSignificativityCriteriaRidge
      pvalRange <- input$transcriptPvalRangeRidge
      fdrRange <- input$transcriptFdrRangeRidge
    }
    if(type == "all" & deg == TRUE & "class"%in%colnames(data_table)){
      integrationSelect <- input$integrationSelectRidgeDEG
      significativityCriteria <- input$significativityCriteriaRidgeDEG
      pvalRange <- input$pvalRangeRidgeDEG
      fdrRange <- input$fdrRangeRidgeDEG
      typeSelect <- input$typeSelectRidgeDEG
      classSelect <- input$classSelectRidgeDEG
      df <- df[df$deg,]
      df <- df[df$class == classSelect,]
    }
    if("class"%in%colnames(data_table)){df <- df[df$class == classSelect,]}
    if(integrationSelect == "gene_genomic_res"){
      df <- df[df$cnv_met == typeSelect,]
      df <- df[!is.na(df$cnv_met), ]
    }
    df <- df[df$omics == integrationSelect,]
    if (significativityCriteria == 'pval'){
      df$significance <- ifelse(df$pval >= pvalRange[1] & df$pval <= pvalRange[2],
                                "Significant", "Not Significant")
    }else{
      df$significance <- ifelse(df$fdr >= fdrRange[1] & df$fdr <= fdrRange[2],
                                "Significant", "Not Significant")
    }
    if (nrow(df) == 0){
      return(NULL)}
    lower_quantile <- quantile(df$coef, 0.001)
    upper_quantile <- quantile(df$coef, 0.999)
    ans <- list(df = df, quantiles = c(lower_quantile, upper_quantile))
    return(ans)
  }) %>% bindEvent(input$genomicIntegrationSelectRidge,
                   input$genomicClassSelectRidge,
                   input$genomicSignificativityCriteriaRidge,
                   input$genomicPvalRangeRidge,
                   input$genomicFdrRangeRidge,
                   input$genomicTypeSelectRidge,
                   input$transcriptIntegrationSelectRidge,
                   input$trascriptClassSelectRidge,
                   input$transcriptSignificativityCriteriaRidge,
                   input$transcriptPvalRangeRidge,
                   input$transcriptFdrRangeRidge,
                   input$integrationSelectRidgeDEG,
                   input$significativityCriteriaRidgeDEG,
                   input$pvalRangeRidgeDEG,
                   input$fdrRangeRidgeDEG,
                   input$typeSelectRidgeDEG,
                   input$classSelectRidgeDEG)
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
    df <- mutate_if(df, is.numeric, ~ round(., 3))
    if(type == "genomic"){
      integrationSelect <- input$genomicIntegrationSelectRidge
      classSelect <- input$genomicClassSelectRidge
      significativityCriteria <- input$genomicSignificativityCriteriaRidge
      pvalRange <- input$genomicPvalRangeRidge
      fdrRange <- input$genomicFdrRangeRidge
      typeSelect <- input$genomicTypeSelectRidge
    }
    if(type == "transcript"){
      integrationSelect <- input$transcriptIntegrationSelectRidge
      classSelect <- input$transcriptClassSelectRidge
      significativityCriteria <- input$transcriptSignificativityCriteriaRidge
      pvalRange <- input$transcriptPvalRangeRidge
      fdrRange <- input$transcriptFdrRangeRidge
    }
    if(type == "all" & deg == TRUE){
      integrationSelect <- input$integrationSelectRidgeDEG
      significativityCriteria <- input$significativityCriteriaRidgeDEG
      pvalRange <- input$pvalRangeRidgeDEG
      fdrRange <- input$fdrRangeRidgeDEG
      typeSelect <- input$typeSelectRidgeDEG
      df <- df[df$deg,]
    }
    if("class"%in%colnames(df) & deg == FALSE){
      df <- df[df$class == classSelect,]
    }
    if(integrationSelect == "gene_genomic_res"){
      df <- df[df$cnv_met == typeSelect,]
      df <- df[!is.na(df$cnv_met), ]
    }
    df <- df[df$omics == integrationSelect,]
    if(significativityCriteria == 'pval'){
      df$significance <- ifelse(df$pval >= pvalRange[1] & df$pval <= pvalRange[2],
                                "Significant", "Not Significant")
    }else{
      df$significance <- ifelse(df$fdr >= fdrRange[1] & df$fdr <= fdrRange[2],
                                "Significant", "Not Significant")
    }
    if (nrow(df) == 0){
      return(NULL)}else{
    return(df)}
  }) %>% bindEvent(input$genomicIntegrationSelectRidge,
                   input$genomicClassSelectRidge,
                   input$genomicSignificativityCriteriaRidge,
                   input$genomicPvalRangeRidge,
                   input$genomicFdrRangeRidge,
                   input$genomicTypeSelectRidge,
                   input$transcriptIntegrationSelectRidge,
                   input$transcriptClassSelectRidge,
                   input$transcriptSignificativityCriteriaRidge,
                   input$transcriptPvalRangeRidge,
                   input$transcriptFdrRangeRidge,
                   input$integrationSelectRidgeDEG,
                   input$significativityCriteriaRidgeDEG,
                   input$pvalRangeRidgeDEG,
                   input$fdrRangeRidgeDEG,
                   input$typeSelectRidgeDEG)
}

#######################################################################
########################################################################
#' @importFrom gtools mixedsort

.prepare_reactive_histo <- function(data_table,
                                    input,
                                    output,
                                    type = "genomic",
                                    deg = FALSE){
  reactive({
    if(type == "genomic"){
      integrationSelect <- input$genomicIntegrationSelectHisto
      typeSelect <- input$genomicTypeSelect
      classSelect <- input$genomicClassSelectHisto
      chrSelect <- input$genomicChrSelectHisto
      significativityCriteria <- input$genomicSignificativityCriteriaHisto
      pvalRange <- input$genomicPvalRangeHisto
      fdrRange <- input$genomicFdrRangeHisto
    }
    if(type == "transcript"){
      integrationSelect <- input$transcriptIntegrationSelectHisto
      classSelect <- input$transcriptClassSelectHisto
      chrSelect <- input$transcriptChrSelectHisto
      significativityCriteria <- input$transcriptSignificativityCriteriaHisto
      pvalRange <- input$transcriptPvalRangeHisto
      fdrRange <- input$transcriptFdrRangeHisto
    }
    if(type == "all" & deg == TRUE){
      integrationSelect <- input$integrationSelectHistoDEG
      classSelect <- input$classSelectHistoDEG
      chrSelect <- input$chrSelectHistoDEG
      significativityCriteria <- input$significativityCriteriaHistoDEG
      pvalRange <- input$pvalRangeHistoDEG
      fdrRange <- input$fdrRangeHistoDEG
      typeSelect <- input$typeSelectDEG
    }
    data_table <- data_table[!is.na(data_table$chr_cov),]
    data_table <- data_table[!is.na(data_table$chr_response),]
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    data_table <- data_table[data_table$omics == integrationSelect,]
    if(integrationSelect == "gene_genomic_res"){
      data_table <- data_table[data_table$cnv_met == typeSelect,]
    }
    if(deg==TRUE){
      data_table <- data_table[data_table$deg,]
    }
    data_table <- data_table[!is.na(data_table$chr_cov),]
    if(chrSelect != "All"){
      data_table <- data_table[data_table$chr_cov == chrSelect,]
    }
    if("class"%in%colnames(data_table)){
      data_table <- data_table[data_table$class == classSelect,]
    }
    if(significativityCriteria == 'pval'){
      data_table$significance <- ifelse(data_table$pval >= pvalRange[1] &
                                          data_table$pval <= pvalRange[2],
                                        "Significant",
                                        "Not Significant")
    }else{
      data_table$significance <- ifelse(data_table$fdr >= fdrRange[1] &
                                          data_table$fdr <= fdrRange[2],
                                        "Significant",
                                        "Not Significant")
    }
    if (nrow(data_table) == 0){
      return(NULL)}else{
    return(data_table)}
  })%>%bindEvent(input$genomicIntegrationSelectHisto,
                 input$genomicTypeSelect,
                 input$genomicClassSelectHisto,
                 input$genomicChrSelectHisto,
                 input$genomicSignificativityCriteriaHisto,
                 input$genomicPvalRangeHisto,
                 input$genomicFdrRangeHisto,
                 input$transcriptIntegrationSelectHisto,
                 input$transcriptClassSelectHisto,
                 input$transcriptChrSelectHisto,
                 input$transcriptSignificativityCriteriaHisto,
                 input$transcriptPvalRangeHisto,
                 input$transcriptFdrRangeHisto,
                 input$integrationSelectHistoDEG,
                 input$classSelectHistoDEG,
                 input$chrSelectHistoDEG,
                 input$significativityCriteriaHistoDEG,
                 input$pvalRangeHistoDEG,
                 input$fdrRangeHistoDEG,
                 input$typeSelectDEG)
}

#######################################################################
########################################################################
#' @importFrom gtools mixedsort

.prepare_reactive_histo_table <- function(data_table,
                                          input,
                                          output,
                                          type = "genomic",
                                          deg = FALSE){
  reactive({
    data_table <- mutate_if(data_table, is.numeric, ~ round(., 3))
    if(type == "genomic"){
      integrationSelect <- input$genomicIntegrationSelectHisto
      typeSelect <- input$genomicTypeSelect
      classSelect <- input$genomicClassSelectHisto
      chrSelect <- input$genomicChrSelectHisto
      significativityCriteria <- input$genomicSignificativityCriteriaHisto
      pvalRange <- input$genomicPvalRangeHisto
      fdrRange <- input$genomicFdrRangeHisto
    }
    if(type == "transcript"){
      integrationSelect <- input$transcriptIntegrationSelectHisto
      classSelect <- input$transcriptClassSelectHisto
      chrSelect <- input$transcriptChrSelectHisto
      significativityCriteria <- input$transcriptSignificativityCriteriaHisto
      pvalRange <- input$transcriptPvalRangeHisto
      fdrRange <- input$transcriptFdrRangeHisto
    }
    if(type == "all" & deg == TRUE){
      integrationSelect <- input$integrationSelectHistoDEG
      classSelect <- input$classSelectHistoDEG
      chrSelect <- input$chrSelectHistoDEG
      significativityCriteria <- input$significativityCriteriaHistoDEG
      pvalRange <- input$pvalRangeHistoDEG
      fdrRange <- input$fdrRangeHistoDEG
      typeSelect <- input$typeSelectDEG
    }
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    data_table <- data_table[data_table$omics == integrationSelect,]
    if(integrationSelect == "gene_genomic_res"){
      data_table <- data_table[data_table$cnv_met == typeSelect,]
    }
    if(deg==TRUE){
      data_table <- data_table[data_table$deg,]
    }
    data_table <- data_table[!is.na(data_table$chr_cov),]
    if(chrSelect != "All"){
      data_table <- data_table[data_table$chr_cov == chrSelect,]
    }
    if("class"%in%colnames(data_table)){
      data_table <- data_table[data_table$class == classSelect,]
    }
    if(significativityCriteria == 'pval'){
      data_table$significance <- ifelse(data_table$pval >= pvalRange[1] &
                                          data_table$pval <= pvalRange[2],
                                        "Significant",
                                        "Not Significant")
    }else{
      data_table$significance <- ifelse(data_table$fdr >= fdrRange[1] &
                                          data_table$fdr <= fdrRange[2],
                                        "Significant",
                                        "Not Significant")
    }
    if (nrow(data_table) == 0){
      return(NULL)}else{
    return(data_table)}
  })%>%bindEvent(input$genomicIntegrationSelectHisto,
                 input$genomicTypeSelect,
                 input$genomicClassSelectHisto,
                 input$genomicChrSelectHisto,
                 input$genomicSignificativityCriteriaHisto,
                 input$genomicPvalRangeHisto,
                 input$genomicFdrRangeHisto,
                 input$transcriptIntegrationSelectHisto,
                 input$transcriptClassSelectHisto,
                 input$transcriptChrSelectHisto,
                 input$transcriptSignificativityCriteriaHisto,
                 input$transcriptPvalRangeHisto,
                 input$transcriptFdrRangeHisto,
                 input$integrationSelectHistoDEG,
                 input$classSelectHistoDEG,
                 input$chrSelectHistoDEG,
                 input$significativityCriteriaHistoDEG,
                 input$pvalRangeHistoDEG,
                 input$fdrRangeHistoDEG,
                 input$typeSelectDEG)
}
#######################################################################
########################################################################
#' @importFrom gtools mixedsort

.prepare_reactive_histo_tf <- function(data_table,
                                       input,
                                       output,
                                       deg = FALSE){
  reactive({
    genes_count_df <- NULL
    if(deg==FALSE){
      integrationSelect <- input$transcriptIntegrationSelectHisto
      classSelect <- input$transcriptClassSelectHisto
      chrSelect <- input$transcriptChrSelectHisto
      significativityCriteria <- input$transcriptSignificativityCriteriaHisto
      pvalRange <- input$transcriptPvalRangeHisto
      fdrRange <- input$transcriptFdrRangeHisto
    }
    if(deg==TRUE){
      integrationSelect <- input$integrationSelectHistoDEG
      classSelect <- input$classSelectHistoDEG
      chrSelect <- input$chrSelectHistoDEG
      significativityCriteria <- input$significativityCriteriaHistoDEG
      pvalRange <- input$pvalRangeHistoDEG
      fdrRange <- input$fdrRangeHistoDEG
    }
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    if(integrationSelect == "tf_res"){
    data_table <- data_table[data_table$omics == 'tf_res', ]
    data_table <- data_table[,
                             colnames(data_table)%in%c("response", "cov",
                                                       "pval","fdr","chr_cov",
                                                       "deg","class")]
    if('class' %in% colnames(data_table)){
      data_table <- data_table[
        data_table$class == classSelect,]
    }
    if(deg==TRUE){
      data_table <- data_table[data_table$deg,]
    }
    if(significativityCriteria == "pval"){
      data_table <- data_table[data_table$pval <= pvalRange,]
    }else{
      data_table <- data_table[data_table$fdr <= fdrRange,]
    }
    genes_count <- table(data_table$cov, data_table$chr_cov)
    genes_count_df <- as.data.frame.table(genes_count)
    genes_count_df <- subset(genes_count_df, Freq != 0)
    colnames(genes_count_df) <- c("TF", "Chromosome", "Count")
    genes_count_df <- genes_count_df[order(-genes_count_df$Count),]
    }
    return(genes_count_df)
  })%>%bindEvent(input$transcriptIntegrationSelectHisto,
                 input$transcriptClassSelectHisto,
                 input$transcriptChrSelectHisto,
                 input$transcriptSignificativityCriteriaHisto,
                 input$transcriptPvalRangeHisto,
                 input$transcriptFdrRangeHisto,
                 input$integrationSelectHistoDEG,
                 input$classSelectHistoDEG,
                 input$chrSelectHistoDEG,
                 input$significativityCriteriaHistoDEG,
                 input$pvalRangeHistoDEG,
                 input$fdrRangeHistoDEG
                 )
}

#######################################################################
########################################################################
#' @importFrom gtools mixedsort

.prepare_reactive_table <- function(data_table,
                                    input,
                                    output){
  reactive({
    data_table <- mutate_if(data_table, is.numeric, ~ round(., 3))
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)

    data_table <- data_table[!is.na(data_table$chr_cov),]
    data_table <- data_table[data_table$omics == input$integrationSelectTable,]
    if('class' %in% colnames(data_table)){
      data_table <- data_table[data_table$class == input$classSelectTable,]}

    if(input$significativityCriteriaTable == 'pval'){
      data_table <- data_table[data_table$pval >= input$pvalRangeTable[1] &
                                 data_table$pval <= input$pvalRangeTable[2],]
    }else{
      data_table <- data_table[data_table$fdr >= input$FDRRangeTable[1] &
                                 data_table$fdr <= input$FDRRangeTable[2],]
    }
    if(input$degSelectTable == 'Only DEGs'){
      data_table <- data_table[data_table$deg,]
      }
    if(input$chrSelectTable != "All"){
      data_table <- data_table[data_table$chr_cov == input$chrSelectTable,]
    }
    return(data_table)
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
#' @importFrom shiny renderText
#' @importFrom DT dataTableOutput

.background_srv <- function(input,
                            output,
                            session,
                            data_gen_enrich,
                            data_tf_enrich,
                            name){

  gen_enr <- .run_bg(FFUN = run_genomic_enrich,
                     input = input,
                     output = output,
                     args = list(model_results = NULL,
                                 qvalueCutoff = 1,
                                 pvalueCutoff = 0.05,
                                 extracted_data = data_gen_enrich,
                                 pAdjustMethod="none"))
  check <- gINTomics:::.check_reactive_bg_enrich(bg_enrich = gen_enr,
                                                 input = input,
                                                 output = output,
                                                 session = session)


  output$gen_enrichment <- renderText({
    check()
  })

  gen_plot <- gINTomics:::.reactive_gen_enrich(bg_enrich = gen_enr,
                                               input = input,
                                               output = output,
                                               session = session)
  output$gen_dotplot <- renderPlotly({
    gen_plot()[["plot"]]
  })
  output$gen_enrich_table <- DT::renderDataTable({
    ans <- gen_plot()[["table"]]
    if(!gen_enr$is_alive()){
    ans <- mutate_if(ans, is.numeric, ~ round(., 3))
    }
  })

  tf_enr <- gINTomics:::.run_bg(FFUN = run_tf_enrich,
                                input = input,
                                output = output,
                                args = list(model_results = NULL,
                                            qvalueCutoff = 1,
                                            pvalueCutoff = 0.05,
                                            extracted_data = data_tf_enrich,
                                            pAdjustMethod="none"))
  check_tf <- gINTomics:::.check_reactive_bg_enrich(bg_enrich = tf_enr,
                                                    input = input,
                                                    output = output,
                                                    session = session)

  output$tf_enrichment <- renderText({
    check_tf()
  })
  tf_plot <- gINTomics:::.reactive_tf_enrich(bg_enrich = tf_enr,
                                             input = input,
                                             output = output,
                                             session = session)

  ns <- NS(name)
  output$tf_dotplot <- renderUI({
    plots <- tf_plot()
    plot_list <- lapply(seq_along(plots), function(i) {
      list(plotlyOutput(ns(paste0(names(plots)[i], "_plot"))),
           HTML(paste0(rep("<br>", 20), collapse = "")),
           tags$div(
             style = 'overflow-x: auto;',
           dataTableOutput(ns(paste0(names(plots)[i], "_table")))),
           HTML(paste0(rep("<br>", 20), collapse = "")))
    })
    plot_list <- as.list(unlist(plot_list, recursive = F))
    do.call(tagList, plot_list)
    return(plot_list)
  })
  observe({
    plots <- tf_plot()
    for (i in 1:length(plots)) {
      output[[paste0(names(plots)[i], "_plot")]] <- renderPlotly({
        plots[[i]][["plot"]]
      })
      output[[paste0(names(plots)[i], "_table")]] <- renderDataTable({
       ans <- plots[[i]][["table"]]
        if(!tf_enr$is_alive()){
          ans <- mutate_if(ans, is.numeric, ~ round(., 3))
        }
      })

    }
  })


    onStop(function(){
      cat(sprintf("Session was closed, killing all background processes"))
        tf_enr$kill_tree()
        gen_enr$kill_tree()
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
                     supervise = TRUE)
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
#' @importFrom clusterProfiler dotplot

.reactive_gen_enrich <- function(bg_enrich,
                                 input,
                                 output,
                                 session){
  reactive({
    if (bg_enrich$is_alive()){
      invalidateLater(millis = 1000, session = session)
    }
    class <- input$genomicClassSelectEnrich
    db <- input$genomicDBSelectEnrich
    type <- input$genomicTypeSelectEnrich
    if (!bg_enrich$is_alive()) {
      data <- bg_enrich$get_result()
      if(sum(c("cnv", "met")%in%names(data))==2){
        ans <- data[[type]][[class]][[db]]
        ans2 <- dot_plotly(ans)
        ans <- list(plot=ans2, table=ans@result)
      }else{
        ans <- data[[class]][[db]]
        ans2 <- dot_plotly(ans)
        ans <- list(plot=ans2, table=ans@result)
      }
      return(ans)
    }
  })
}
#######################################################################
########################################################################
#' @importFrom plotly renderPlotly

.reactive_tf_enrich <- function(bg_enrich,
                                input,
                                output,
                                session){
  reactive({
    if (bg_enrich$is_alive()){
      invalidateLater(millis = 1000, session = session)
    }
    class <- input$transcriptionalClassSelectEnrich
    db <- input$transcriptionalDBSelectEnrich
    if (!bg_enrich$is_alive()) {
      data <- bg_enrich$get_result()
      tf <- names(data[[class]])
      ans <- lapply(tf, function(x){
        ans <- data[[class]][[x]][[db]]
        if(!is.null(ans)){
          ans2 <- dot_plotly(ans)
          ans <- list(plot=ans2, table=ans@result)
          return(ans)
        }else{return(NULL)}
      })
      names(ans) <- tf
      return(ans)
    }
  })
}
#######################################################################
########################################################################
#' @importFrom shiny.gosling arrange_views

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
                                   ),
                                   assembly="hg38")

    return(arranged_view)
  })%>%bindEvent(input$circosClass,
                 input$circosLayout,
                 input$circosType,
                 input$circosChr)
}
