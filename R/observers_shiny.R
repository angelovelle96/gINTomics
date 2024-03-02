#' @importFrom visNetwork renderVisNetwork
#' @importFrom visNetwork visHierarchicalLayout
#' @importFrom visNetwork  visIgraphLayout

.render_reactive_network <- function(reactive_network,
                                     input,
                                     output,
                                     deg = FALSE){
  observe({
    if(deg==FALSE){
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
    } else {
      output$networkPlotDEG <- renderVisNetwork({
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
    }
  })%>%bindEvent(input$layoutNetwork,
                 input$physics)
}

################################################################
#################################################################

.select_deg_network <- function(data_table,
                                network_data,
                                input,
                                output,
                                deg = FALSE){
  reactive({
    if(deg==FALSE){
      numNodes <- input$numNodes
      significativityCriteria <- input$significativityCriteriaNetwork
      pval <- input$pvalNetwork
      fdr <- input$fdrNetwork

      ans <- network_data
      #data_table <- data_table[data_table$cov != "(Intercept)",]
      if(significativityCriteria == "pval"){
        data_table <- data_table[data_table$pval <= pval,]
      } else {
        data_table <- data_table[data_table$fdr <= fdr,]
      }
      ssel <- unique(data_table$response)
      ssel <- head(ssel, numNodes)
      ans$nodes <- ans$nodes[ans$nodes$id%in%ssel,]
      ans$edges <- ans$edges[ans$edges$from%in%ssel,]
      ans$edges <- ans$edges[ans$edges$to%in%ssel,]
      nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
      ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges,]
    }

    if(deg==TRUE & "class"%in%colnames(data_table)){
      numNodes <- input$numNodesDEG
      significativityCriteria <- input$significativityCriteriaNetworkDEG
      pval <- input$pvalNetworkDEG
      fdr <- input$fdrNetworkDEG

      ans <- network_data
      #data_table <- data_table[data_table$cov != "(Intercept)",]
      if(significativityCriteria == "pval"){
        data_table <- data_table[data_table$pval <= pval,]
      } else {
        data_table <- data_table[data_table$fdr <= fdr,]
      }
      ssel <- unique(data_table$response[data_table$deg])
      ssel <- head(ssel, numNodes)
      ans$nodes <- ans$nodes[ans$nodes$id%in%ssel,]
      ans$edges <- ans$edges[ans$edges$from%in%ssel,]
      ans$edges <- ans$edges[ans$edges$to%in%ssel,]
      # nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
      # ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges,]
    }else{return(NULL)}

    return(ans)
  })%>%bindEvent(input$layoutNetwork,
                 input$numNodes,
                 input$significativityCriteriaNetwork,
                 input$pvalNetwork,
                 input$fdrNetwork,
                 input$physics,
                 input$layoutNetworkDEG,
                 input$numNodesDEG,
                 input$significativityCriteriaNetworkDEG,
                 input$pvalNetworkDEG,
                 input$fdrNetworkDEG,
                 input$physicsDEG)
}

################################################################
#################################################################
.prepare_reactive_venn <- function(data_table,
                                   input,
                                   output,
                                   deg = FALSE){
  reactive_venn <- reactive({
    if(deg==FALSE & "gene_genomic_res"%in%unique(data_table$omics)) {
      classSelect <- input$classSelectVenn
      pvalRange <- input$pvalRangeVenn
      fdrRange <- input$fdrRangeVenn
      significativityCriteria <- input$significativityCriteriaVenn
      data_table <- data_table[data_table$omics == 'gene_genomic_res', ]
      if('class'%in%colnames(data_table)){
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
      if("cnv_genomic_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="cnv_genomic_res",]
      }
      if("met_genomic_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="met_genomic_res",]
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
      if("cnv_genomic_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="cnv_genomic_res",]
      }
      if("met_genomic_res"%in%unique(data_table$omics)){
        data_venn <- data_table[data_table$omics=="met_genomic_res",]
      }
    }
    # if(!"gene_genomic_res"%in%unique(data_table$omics)){return(NULL)}
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
      top_peaks <- data_table[with(data_table,
                                   order(pval, coef)),]
      data_table <- rbind(top_peaks, data_table[with(data_table,
                                                     order(-coef, pval)), ][1:10,])
      data_table$pval_fdr <- -log10(data_table$pval)
    }else{
      data_table["group"] <- "Not Significant"
      data_table[data_table$fdr <= fdrRange, 'group'] <- "Significant"
      top_peaks <- data_table[with(data_table,
                                   order(fdr, coef)),]
      data_table <- rbind(top_peaks, data_table[with(data_table,
                                                     order(-coef, fdr)),][1:10,])
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
    }
    if("class"%in%colnames(data_table) & deg == FALSE){
      df_heatmap <- multiomics_integration[[integrationSelect]][[classSelect]]$data$response_var
      data_table <- data_table[data_table$omics == integrationSelect,]
      data_table <- data_table[data_table$class == classSelect,]
    }
    if(!"class"%in%colnames(data_table) & deg == FALSE){
      df_heatmap <- multiomics_integration[[integrationSelect]]$data$response_var
      data_table <- data_table[data_table$omics == integrationSelect,]
    }else{return(NULL)}
    if("class"%in%colnames(data_table) & deg == TRUE){
      df_heatmap <- multiomics_integration[[integrationSelect]][[classSelect]]$data$response_var
      data_table <- data_table[data_table$omics == integrationSelect,]
      data_table <- data_table[data_table$class == classSelect,]
      data_table <- data_table[data_table$deg,]
      df_heatmap <- df_heatmap[, unique(data_table$response)]
    }
    df_heatmap_t <- t(as.matrix(df_heatmap))
    if (integrationSelect == "gene_genomic_res"){
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
      if(significativityCriteria == 'pval'){
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= pvalRange |
                                       df_heatmap_t$pval_met <= pvalRange,]
      }else{
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= fdrRange |
                                       df_heatmap_t$fdr_met <= fdrRange,]
      }
      if (nrow(df_heatmap_t) == 0){
        return(NULL)
      }
      top_met <- df_heatmap_t %>% arrange(desc(abs(met))) %>% head(numTopCNV)
      top_cnv <- df_heatmap_t %>% arrange(desc(abs(cnv))) %>% head(numTopMET)
      expr_top <- rbind(top_cnv, top_met[!rownames(top_met)%in%rownames(top_cnv),])
      expr_top_subset <- expr_top[, -c((ncol(expr_top) - 5):ncol(expr_top))]
      set.seed(123)
      row_ha <- rowAnnotation(coef_cnv = expr_top$cnv, coef_met = expr_top$met)
      ht <- Heatmap(expr_top_subset, right_annotation = row_ha)
      ht <- draw(ht)
      if(deg == FALSE){
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')
      }else{
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmapDEG')
      }
    }
    if(integrationSelect == "cnv_genomic_res"){
      tmp <-  data.frame(cnv=data_table$coef[data_table$omics == 'cnv'],
                         pval_cnv=data_table$pval[data_table$omics == 'cnv'],
                         fdr_cnv=data_table$fdr[data_table$omics == 'cnv'],
                         row.names = data_table$response[data_table$omics == 'cnv'])
      tmp2 <- unique(rownames(tmp))
      data_table <- cbind(cnv=tmp[tmp2,])
      colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
      rownames(data_table) <- tmp2
      df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
      if (significativityCriteria == 'pval'){
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= pvalRange,]
      }else{
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= fdrRange,]
      }
      if (nrow(df_heatmap_t) == 0){
        return(NULL)}
      top_cnv <- df_heatmap_t %>%
        arrange(desc(abs(cnv))) %>%
        head(numTopCNVonly)
      expr_top <- top_cnv
      set.seed(123)
      row_ha <- rowAnnotation(coef_cnv=expr_top$cnv)
      ht <- Heatmap(as.matrix(expr_top, right_annotation=row_ha))
      ht = draw(ht)
      if(deg == FALSE){
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')
      }else{
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmapDEG')
      }
    }
    if(integrationSelect == "met_genomic_res"){
      tmp <-  data.frame(met=data_table$coef[data_table$omics == 'met'],
                         pval_met=data_table$pval[data_table$omics == 'met'],
                         fdr_met=data_table$fdr[data_table$omics == 'met'],
                         row.names = data_table$response[data_table$omics == 'met'])
      tmp2 <- unique(rownames(tmp))
      data_table <- cbind(met=tmp[tmp2,])
      colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
      rownames(data_table) <- tmp2
      df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
      if(significativityCriteria == 'pval'){
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_met <= pvalRange,]
      }else{
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_met <= fdrRange,]
      }
      if (nrow(df_heatmap_t) == 0){
        return(NULL)}
      top_met <- df_heatmap_t %>%
        arrange(desc(abs(met))) %>%
        head(numTopMETonly)
      expr_top <- top_met
      set.seed(123)
      row_ha <- rowAnnotation(coef_met=expr_top$met)
      ht <- Heatmap(as.matrix(expr_top, right_annotation=row_ha))
      ht = draw(ht)
      if(deg == FALSE){
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')
      }else{
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmapDEG')
      }
    }
    if(integrationSelect == "mirna_cnv_res"){
      tmp <-  data.frame(mirna_cnv=data_table$coef[data_table$omics == 'mirna_cnv_res'],
                         pval_mirna_cnv=data_table$pval[data_table$omics == 'mirna_cnv_res'],
                         fdr_mirna_cnv=data_table$fdr[data_table$omics == 'mirna_cnv_res'],
                         row.names = data_table$response[data_table$omics == 'mirna_cnv_res'])
      tmp2 <- unique(rownames(tmp))
      data_table <- cbind(mirna_cnv=tmp[tmp2,])
      colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
      rownames(data_table) <- tmp2
      df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
      if(significativityCriteria == 'pval'){
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_mirna_cnv <= pvalRange,]
      }else{
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_mirna_cnv <= fdrRange,]
      }
      if (nrow(df_heatmap_t) == 0){
        return(NULL)}
      df_heatmap_t <- as.data.frame(df_heatmap_t)
      top_mirna_cnv <- df_heatmap_t %>%
        arrange(desc(abs(mirna_cnv)))  %>%
        head(numTopMiCNV)
      expr_top <- top_mirna_cnv
      expr_top_subset <- expr_top[, -c((ncol(expr_top) - 3):ncol(expr_top))]
      set.seed(123)
      row_ha <- rowAnnotation(coef_mirna=expr_top$mirna_cnv)
      ht <- Heatmap(expr_top_subset, right_annotation=row_ha)
      ht = draw(ht)
      if(deg == FALSE){
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')
      }else{
        ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmapDEG')
      }
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
                 input$classSelectHeatmapDEG)
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
    chr_order <- gtools::mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    data_table <- data_table[data_table$omics == integrationSelect,]
    if(integrationSelect == "gene_genomic_res"){
      data_table <- data_table[data_table$omics == integrationSelect,]
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
      data_table <- data_table[data_table$omics == integrationSelect,]
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
                                       output){
  reactive({
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    data_table <- data_table[data_table$omics == 'tf_res', ]
    data_table <- data_table[,
                             colnames(data_table)%in%c("response", "cov",
                                                       "pval","fdr","chr_cov",
                                                       "deg","class")]
    df_filtered_histo_tf <- data_table
    if('class' %in% colnames(df_filtered_histo_tf)){
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
#' @importFrom shiny renderText

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
    gen_plot()[["table"]]
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
           DT::dataTableOutput(ns(paste0(names(plots)[i], "_table"))),
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
        plots[[i]][["table"]]
      })

    }
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
