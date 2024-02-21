
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
    if(input$deg==TRUE){
      ans <- network_data
      ssel <- unique(data_table$response[data_table$deg==TRUE])
      ans$nodes <- ans$nodes[ans$nodes$id%in%ssel,]
      ans$edges <- ans$edges[ans$edges$from%in%ssel,]
      ans$edges <- ans$edges[ans$edges$to%in%ssel,]
    }else{ans <- network_data}



    return(ans)
  })%>%bindEvent(input$layoutNetwork,
                 input$deg,
                 input$physics)
}


################################################################
#################################################################
.prepare_reactive_venn <- function(data_table,
                                   input,
                                   output){
  reactive_venn <- reactive({
    if(input$integrationSelectVenn=='gene_genomic_res') {
      data_table <- data_table[data_table$omics == 'gene_genomic_res', ]
    }else{
      return(NULL)}
    if('class' %in% names(data_table)){data_table <- data_table[data_table$class == input$classSelectVenn, ]}
    # Leggi il range di pval e in base a questo prendi i cnv e met significativi
    if(input$degSelectVenn == 'Only DEGs' & 'class' %in% names(data_table)){data_table <- data_table[data_table$deg == 'TRUE', ] }

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
                                         input$integrationSelectVenn,
                                         input$significativityCriteriaVenn,
                                         input$degSelectVenn)
}
#####################################################################
######################################################################

.prepare_reactive_volcano <- function(data_table,
                                      input,
                                      output){
  reactive({
    if(!"class"%in%colnames(data_table)) data_table$class <- data_table$omics
    if (input$integrationSelectVolcano == 'gene_genomic_res'){
      data_volcano <- data_table[data_table$omics == 'gene_genomic_res',]
      data_volcano <- data_volcano[
        data_volcano$cnv_met == input$genomicTypeSelectVolcano,]
    } else {
      data_volcano <- data_table[
        data_table$omics == input$integrationSelectVolcano,]
    }
    if(input$degSelectVolcano == 'Only DEGs' & 'class' %in% names(data_table)){

      data_volcano <- data_volcano[data_volcano$deg == TRUE, ]

    }

    if(input$significativityCriteriaVolcano == 'pval'){

      data_volcano["group"] <- "NotSignificant"

      data_volcano[data_volcano$pval <= input$pvalRangeVolcano, 'group'] <- "Significant"

      top_peaks <- data_volcano[with(data_volcano,
                                     order(pval, coef)),][1:input$numTopGenesVolcano,]
      data_volcano <- rbind(top_peaks, data_volcano[with(data_volcano,
                                                         order(-coef, pval)), ][1:10,])


      data_volcano$pval_fdr <- -log10(data_volcano$pval)
    }else{

      data_volcano[data_volcano$fdr <= input$FDRRangeVolcano, 'group'] <- "Significant"

      top_peaks <- data_volcano[with(data_volcano,
                                     order(fdr, coef)),][1:input$numTopGenesVolcano,]
      data_volcano <- rbind(top_peaks, data_volcano[with(data_volcano,
                                                         order(-coef, fdr)),][1:10,])


      data_volcano$pval_fdr <- -log10(data_volcano$fdr)
    }

    return(data_volcano)
  })%>%bindEvent(input$genomicTypeSelectVolcano,
                 input$integrationSelectVolcano,
                 input$FDRRangeVolcano,
                 input$numTopGenesVolcano,
                 input$pvalRangeVolcano,
                 input$degSelectVolcano,
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

  if(!is.null(input$classSelectHeatmap)){ df_heatmap <- multiomics_integration[[input$integrationSelectHeatmap]][[input$selectClassHeatmap]]$data$response_var
                                          data_table <- data_table[data_table$class == input$selectClassHeatmap,]}
  if(input$degSelectHeatmap == "Only DEGs"){data_table <- data_table[data_table$deg == TRUE,]}
  if(is.null(input$classSelectHeatmap)){ df_heatmap <- multiomics_integration[[input$integrationSelectHeatmap]]$data$response_var}
  data_table <- data_table[data_table$cov != '(Intercept)',]
  df_heatmap_t <- t(df_heatmap)

  if(input$integrationSelectHeatmap == "gene_genomic_res"){
    data_table <- data_table[data_table$omics == 'gene_genomic_res',]
    tmp <-  data.frame(cnv=data_table$coef[data_table$cnv_met == 'cnv'],
                       pval_cnv=data_table$pval[data_table$cnv_met == 'cnv'],
                       fdr_cnv=data_table$fdr[data_table$cnv_met == 'cnv'],
                       row.names = data_table$response[data_table$cnv_met == 'cnv'])

    tmp2 <- data.frame(met=data_table$coef[data_table$cnv_met == 'met'],
                       pval_met=data_table$pval[data_table$cnv_met == 'met'],
                       fdr_met=data_table$fdr[data_table$cnv_met == 'met'],
                       row.names = data_table$response[data_table$cnv_met == 'met'])

    tmp3 <- unique(c(rownames(tmp), rownames(tmp2)))

    data_table <- cbind(cnv=tmp[tmp3,], met=tmp2[tmp3,])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp3
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
    if (input$significativityCriteriaHeatmap == 'pval'){
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= input$pvalRangeHeatmap &
                                     df_heatmap_t$pval_met <= input$pvalRangeHeatmap, ]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= input$FDRRangeHeatmap &
                                     df_heatmap_t$fdr_met <= input$FDRRangeHeatmap, ]
    }
    top_met <- df_heatmap_t %>%
      arrange(desc(abs(met))) %>%
      head(input$numTopGenesHeatmap/2)
    top_cnv <- df_heatmap_t %>%
      arrange(desc(abs(cnv))) %>%
      head(input$numTopGenesHeatmap/2)
    expr_top <- rbind(top_cnv, top_met)

    row_ha <- rowAnnotation(coef_cnv=expr_top$cnv, coef_met=expr_top$met)
    ht <- ComplexHeatmap:::Heatmap(expr_top, right_annotation=row_ha)
    ht = draw(ht)
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
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= input$pvalRangeHeatmap &
                                     df_heatmap_t$pval_cnv <= input$pvalRangeHeatmap, ]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= input$FDRRangeHeatmap &
                                     df_heatmap_t$fdr_cnv <= input$FDRRangeHeatmap, ]
    }

    top_cnv <- df_heatmap_t %>%
      arrange(desc(abs(cnv))) %>%
      head(input$numTopGenesHeatmap)
    expr_top <- top_cnv

    row_ha <- rowAnnotation(foo1=runif(10), foo2=runif(10))
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
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_met <= input$pvalRangeHeatmap &
                                     df_heatmap_t$pval_met <= input$pvalRangeHeatmap, ]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_met <= input$FDRRangeHeatmap &
                                     df_heatmap_t$fdr_met <= input$FDRRangeHeatmap, ]
    }

    top_met <- df_heatmap_t %>%
      arrange(desc(abs(met))) %>%
      head(input$numTopGenesHeatmap)
    expr_top <- top_met

    row_ha <- rowAnnotation(foo1=runif(10), foo2=runif(10))
    ht <- Heatmap(as.matrix(expr_top, right_annotation=row_ha))
    ht = draw(ht)
    ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')

  }
  if(input$integrationSelectHeatmap == "mirna_cnv_res"){

    data_table <- data_table[data_table$cov != "(Intercept)",]
    data_table <- data_table[data_table$omics == 'mirna_cnv_res',]
    tmp <-  data.frame(mirna_cnv=data_table$coef[data_table$omics == 'mirna_cnv_res'],
                       pval_mirna_cnv=data_table$pval[data_table$omics == 'mirna_cnv_res'],
                       fdr_mirna_cnv=data_table$fdr[data_table$omics == 'mirna_cnv_res'],
                       row.names = data_table$response[data_table$omics == 'mirna_cnv_res'])
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(met=tmp[tmp2,])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t),])
    if (input$significativityCriteriaHeatmap == 'pval'){
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_mirna_cnv <= input$pvalRangeHeatmap &
                                     df_heatmap_t$pval_mirna_cnv <= input$pvalRangeHeatmap, ]
    } else {
      df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_mirna_cnv <= input$FDRRangeHeatmap &
                                     df_heatmap_t$fdr_mirna_cnv <= input$FDRRangeHeatmap, ]
    }

    df_heatmap_t <- as.data.frame(df_heatmap_t)
    top_mirna_cnv <- df_heatmap_t %>%
    arrange(desc(abs(mirna_cnv)))  %>%
    head(input$numTopGenesHeatmap)
    expr_top <- top_mirna_cnv

    row_ha <- rowAnnotation(coef_mirna=expr_top$mirna_cnv)
    ht <- ComplexHeatmap:::Heatmap(expr_top, right_annotation=row_ha)
    ht = draw(ht)
    ht2 <- makeInteractiveComplexHeatmap(input, output, session, ht, 'heatmap')

  }
  })%>%bindEvent(input$integrationSelectHeatmap,
                 input$numTopGenesHeatmap,
                 input$FDRRangeHeatmap,
                 input$pvalRangeHeatmap,
                 input$significativityCriteriaHeatmap,
                 input$selectClassHeatmap)
}

######################################################################
######################################################################

.prepare_reactive_ridge <- function(data_table,
                                      input,
                                      output){

  reactive({
    df <- data_table
    if ('class' %in% names(data_table)) {
      df <- df[df$class == input$classSelectRidge, ]
    }
    if('gene_genomic_res' %in% colnames(data_table)){
      df <- df[df$cnv_met == input$genomicTypeSelectRidge, ]
    }
    if(input$degSelectRidge == 'Only DEGs'){df <- df[df$deg == 'TRUE', ]}
    df <- df[df$omics == input$integrationSelectRidge, ]
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
                 input$integrationSelectRidge,
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
    chr_order <- gtools::mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    df_filtered_histo <- data_table[data_table$omics == input$integrationSelectHisto, ]
    if(input$chrSelectHisto != "All"){
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$chr_cov == input$chrSelectHisto, ]
    }
    if('class' %in% names(df_filtered_histo)){df_filtered_histo <- df_filtered_histo[df_filtered_histo$class == input$classSelectHisto, ]
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
      df_filtered_histo <- df_filtered_histo[df_filtered_histo$deg == 'TRUE', ]
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
      df_filtered_histo_tf <- df_filtered_histo_tf[df_filtered_histo_tf$deg == TRUE, ]
    }
    df_filtered_histo_tf <- df_filtered_histo_tf[df_filtered_histo_tf$pval <= 0.05, ]
    genes_count <- table(df_filtered_histo_tf$cov, df_filtered_histo_tf$chr_cov)
    genes_count_df <- as.data.frame.table(genes_count)
    genes_count_df <- subset(genes_count_df, Freq != 0)
    colnames(genes_count_df) <- c("TF", "Chromosome", "Count") # Rinomina le colonne
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
      filtered_df <- data_table[data_table$omics == input$integrationSelectTable,]
      if('class' %in% names(filtered_df)){
        filtered_df <- filtered_df[filtered_df$class == input$classSelectTable,]}
      if(input$significativityCriteriaTable == 'pval'){
        filtered_df <- filtered_df[filtered_df$pval >= input$pvalRangeTable[1] &
                                     filtered_df$pval <= input$pvalRangeTable[2],]
      }else{
        filtered_df <- filtered_df[filtered_df$fdr >= input$FDRRangeTable[1] &
                                     filtered_df$fdr <= input$FDRRangeTable[2],]
      }
      if(input$degSelectTable == 'Only DEGs'){
        filtered_df <- filtered_df[filtered_df$deg == 'TRUE',]}
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
                            data_gen_enrich){

  gen_enr <- gINTomics:::.reactive_bg(FFUN = run_genomic_enrich,
                                  data_gen_enrich = data_gen_enrich,
                                  input = input,
                                  output = output,
                                  args = list(model_results = NULL,
                                              qvalueCutoff = 1,
                                              pvalueCutoff = 1,
                                              extracted_data = data_gen_enrich))
  check <- gINTomics:::.check_reactive_bg_enrich(reactive_enrich = gen_enr,
                                                 input = input,
                                                 output = output,
                                                 session = session)

  output$gen_enrichment <- renderText({
    check()
  })

  gen_plot <- gINTomics:::.reactive_gen_dotplot(reactive_enrich = gen_enr,
                                                 input = input,
                                                 output = output,
                                                 session = session)
  gINTomics:::.render_gen_dotplot(gen_plot=gen_plot,
                                  data_gen_enrich = data_gen_enrich,
                                  input = input,
                                  output = output,
                                  session = session)

}


#######################################################################
########################################################################

.reactive_bg <- function(FFUN,
                          args,
                          data_gen_enrich,
                          input,
                          output){
  reactive({
    ans <- callr::r_bg(func = FFUN,
                           args = args,
                           supervise = T)
    return(ans)
  })
}

#######################################################################
########################################################################

.check_reactive_bg_enrich <- function(reactive_enrich,
                                        input,
                                        output,
                                        session){
  reactive({
    invalidateLater(millis = 1000, session = session)
    if (reactive_enrich()$is_alive()) {
      x <- "Enrichment running in background, this may take several minutes"
    } else {
      x <- "Enrichment completed"
    }
    return(x)
  })
}

#######################################################################
########################################################################

.reactive_gen_dotplot <- function(reactive_enrich,
                                      input,
                                      output,
                                      session){
  reactive({
    if (reactive_enrich()$is_alive()){
      invalidateLater(millis = 1000, session = session)
      }
    if (!reactive_enrich()$is_alive()) {
      data <- reactive_enrich()$get_result()
      if(sum(c("cnv", "met")%in%names(data))){
      cnv <- data[["cnv"]][[1]][[1]]
      cnv <- clusterProfiler::dotplot(cnv)
      cnv <- ggplotly(cnv)
      met <- data[["met"]][[1]][[1]]
      met <- clusterProfiler::dotplot(met)
      met <- ggplotly(met)
      ans <- list(cnv=cnv, met=met)
      }else{
        ans <- data[[1]][[1]]
        ans <- clusterProfiler::dotplot(ans)
        cnv <- ggplotly(ans)
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

    arranged_views <- list()    
  tracks <- .create_tracks(data_table, data)

if(c("cnv_track","met_track") %in% tracks){

  composed_view_genomic <- .create_composed_view(cnv_track,
                                                met_track,
                                                cyto_track,
                                                expr_track)
  arranged_view_circos_genomic <- arrange_views(
    title = 'Interactive Circos',
    subtitle = 'subtitle',
    views = composed_view_genomic,
    # xDomain = list(chromosome = 'chr1')
  )
  arranged_views <- append(arranged_views, arranged_view_circos_genomic)

  }
if(!(cnv_track %in% tracks)){composed_view_met <- .create_composed_view(met_track,
                                                                         cyto_track,
                                                                         expr_track)
arranged_view_circos_cnv_gene <- arrange_views(
  title = 'Interactive Circos',
  subtitle = 'subtitle',
  views = composed_view_cnv_gene,
  # xDomain = list(chromosome = 'chr1')
)
arranged_views <- append(arranged_views, arranged_view_circos_cnv_gene)
}

if(!(met_track %in% tracks)){composed_view_cnv <- .create_composed_view(cnv_track,
                                                                        cyto_track,
                                                                        expr_track)
arranged_view_circos_met_gene <- arrange_views(
  title = 'Interactive Circos',
  subtitle = 'subtitle',
  views = composed_view_met_gene,
  # xDomain = list(chromosome = 'chr1')
)
arranged_views <- append(arranged_views, arranged_view_circos_met_gene)
}


    return(arranged_view_circos)
  })
}

.create_composed_view <- function(tracks, layout) {

  composed_views <- list()

  if(c(track_expr, track_cnv, track_met)%in%tracks){

    composed_view_circos_genomic <- compose_view(
      multi = TRUE,
      layout = layout,
      tracks = add_multi_tracks(track_expr,
                                track_cnv,
                                track_met,
                                track_cyto),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail"
    )
    composed_views <- append(composed_views, composed_view_circos_genomic)
  }

if(!(track_cnv%in%tracks)){

  composed_view_circos_met_gene <- compose_view(
    multi = TRUE,
    layout = layout,
    tracks = add_multi_tracks(track_expr,
                              track_met,
                              track_cyto),
    alignment = 'stack',
    spacing = 0.01,
    linkingId = "detail"
  )
  composed_views <- append(composed_views, composed_view_circos_met_gene)
}
  if(!(track_met%in%tracks)){

    composed_view_circos_cnv_gene <- compose_view(
      multi = TRUE,
      layout = layout,
      tracks = add_multi_tracks(track_expr,
                                track_cnv,
                                track_cyto),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail"
    )
    composed_views <- append(composed_views, composed_view_circos_cnv_gene)
  }

  return(composed_views)
}
