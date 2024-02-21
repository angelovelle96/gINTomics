.prepare_network <- function(data_table){

  all <- subset(data_table, omics %in% c('tf_res',
                                         'tf_mirna_res',
                                         'mirna_target_res'),
                select = c('response',
                           'cov',
                           'coef',
                           'pval',
                           'omics'))

  all <- all[order(all$pval), ] # modificare
  all <- head(all, 600)
  # ottengo tutti i geni(mirna/tf/targets) in modo da avere un nodo per ogni elemento unico
  nodes <- data.frame(gene = c(all$cov, all$response))  #controllare doppioni
  #ottengo gli edges ovvero le coppie di interazioni (from-to)
  edges <- subset(all, select = c('cov', 'response', 'pval', 'coef'))

  nodes <- unique(nodes)
  names(nodes) <- 'id'
  names(edges) <- c('from', 'to', 'pval', 'coef')

  # ricordare che un tf puó anche essere un target
  tf_list <- all[all$omics=='tf_res', 'cov']
  tf_list2 <- all[all$omics=='tf_mirna_res', 'cov']
  mirna_list <- all[all$omics=='tf_mirna_res', 'response']
  mirna_list2 <- all[all$omics=='mirna_target_res', 'cov']
  target_list <- all[all$omics=='tf_res', 'response']
  target_list2 <- all[all$omics=='tf_mirna_res', 'cov']

  #unisco le liste
  tf_list <- unique(c(tf_list, tf_list2))
  mirna_list <- unique(c(mirna_list, mirna_list2))
  target_list <- unique(c(target_list, target_list2))

  #assegnazione dinamica in base ai dati che ho nel dataframe nodes

  nodes$label <- nodes$id
  nodes$shape <- ifelse(nodes$id %in% tf_list,
                        'diamond',
                        ifelse(nodes$id %in% mirna_list,
                               'triangle',
                               ifelse(nodes$id %in% target_list,
                                      'circle',
                                      'diamond')))
  nodes$title <- nodes$id
  nodes$shadow = TRUE
  nodes$color <- ifelse(nodes$id %in% tf_list,
                        '#FFBA01',
                        ifelse(nodes$id %in% mirna_list,
                               '#9969C7',
                               ifelse(nodes$id %in% target_list,
                                      '#CCE7C9',
                                      '#FFBA01')))
  nodes$width <- 10

  # define edges (TF/mirna and linked genes)
  scaling_factor <- 10
  edges$width <- ifelse(abs(edges$coef) > 0,
                        abs(edges$coef) * scaling_factor,
                        4) #in base al pval o effect size
  edges$color <- ifelse(edges$coef > 0,
                        '#4169E1',
                        ifelse(edges$coef < 0,
                               '#ED5564',
                               'black'))
  edges$length <- 500
  edges$title <- paste('pval:',
                       edges$pval,
                       'coef:', edges$coef)

  legend_nodes <- data.frame(label = c('TF',
                                       'Target',
                                       'miRNA'),
                             shape = c('diamond',
                                       'circle',
                                       'triangle'),
                             color = c('#FFBA01',
                                       '#CCE7C9',
                                       '#9969C7'))

  legend_edges <- data.frame(color = c('#4169E1',
                                       '#ED5564',
                                       'black'),
                             label = c('UPregulate',
                                       'DOWNregulate',
                                       'no-effect'),
                             arrows = c('to',
                                        'to',
                                        'to'))
  return(list(nodes=nodes,
              edges=edges,
              legend_edges=legend_edges,
              legend_nodes=legend_nodes))
}

#########################################################################
#########################################################################
.build_network <- function(nodes,
                           edges,
                           legend_nodes,
                           legend_edges){
      visNetwork(nodes,
                 edges,
                 width = "100%") %>%
        visGroups(groupname = 'TF',
                  color = '#FFBA01') %>%
        visGroups(groupname = 'Target',
                  color = '#CCE7C9') %>%
        visGroups(groupname = 'miRNA',
                  color = '#9969C7') %>%
        visLegend(addEdges = legend_edges,
                  addNodes = legend_nodes,
                  useGroups = FALSE,
                  width = 0.2,
                  position = 'right',
                  main = 'Legend') %>%
        visEdges(shadow = TRUE,
                 smooth = FALSE,
                 arrows = list(to = list(enabled = TRUE)),
                 color = list(color = "black", highlight = "red")) %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           degree = 2,
                                           hover = TRUE),
        nodesIdSelection = TRUE,
        manipulation = TRUE) %>%
        visLayout(randomSeed = 20, improvedLayout = TRUE)
}

#########################################################################
#########################################################################

.build_venn <- function(venn_data){
    cnv_sign_genes <- venn_data$cnv_sign_genes
    met_sign_genes <- venn_data$met_sign_genes
    list_venn <- list(CNV = cnv_sign_genes$cov,
                      MET = met_sign_genes$cov)
    venn_diagram <- ggVennDiagram(list_venn, show_intersect = TRUE)
    return(venn_diagram)
}
#########################################################################
#########################################################################

.render_venn <- function(reactive_venn){
  renderPlotly({
    venn_data <- reactive_venn()
    venn_diagram <- .build_venn(venn_data)
    plotly_venn <- ggplotly(venn_diagram)
    return(plotly_venn)
  })
}

##########################################################################
##########################################################################

.render_venn_table <- function(reactive_venn){
  renderDataTable({
    venn_data <- reactive_venn()
    common_genes <- base::intersect(unlist(venn_data$cnv_sign_genes), unlist(venn_data$met_sign_genes))
    data.frame(Genes = common_genes)
  })
}
############################################################################
############################################################################

.build_volcano <- function(volcano_data){
    plot_ly(volcano_data,
            x = ~coef,
            y = ~pval_fdr,
            mode = 'markers',
            color =  ~group,
            symbol = ~class,
            text =  ~paste("Group:", group, "<br>",
                           "Class:", class,"<br>",
                           "Name:", cov, "<br>",
                           "Pval/FDR:", pval, "<br>",
                           "coef", coef),
            textposition = 'top right') %>%
      layout(title = "Volcano Plot") %>%
      layout(width = 1000,
             height = 700)
}
#######################################################################
#######################################################################

.render_volcano <- function(reactive_volcano){
  renderPlotly({
    volcano_data <- reactive_volcano()
    volcano_plot <- .build_volcano(volcano_data)
    return(volcano_plot)
  })
}


############################################################################
############################################################################

.build_ridge <- function(ridge_data,
                         quantiles){

  ggplot(ridge_data,
         aes(x = coef,
             y = significance,
             fill = factor(significance))) +
    geom_density_ridges(jittered_points = TRUE,
                        quantile_lines = TRUE,
                        vline_size = 1,
                        vline_color = "red",
                        point_size = 0.4,
                        point_alpha = 1,
                        position = position_raincloud(adjust_vlines = TRUE)) +
    labs(title = "Ridgeline Plot",
         x = "Value",
         y = "Significativity") +
    theme_minimal() +
    theme_ridges() +
    scale_x_continuous(limits = quantiles)
}

############################################################################
############################################################################

.render_ridge <- function(reactive_ridge){
  renderPlot({
    tmp <- reactive_ridge()
    ridge_data <- tmp$df
    quantiles <- tmp$quantiles
    ridge_plot <- .build_ridge(ridge_data = ridge_data,
                               quantiles = quantiles)
    return(ridge_plot)
  })
}


############################################################################
############################################################################

.build_histo <- function(histo_data){
  ggplot(histo_data,
         aes(x = factor(chr_cov), fill = significance)) +
    geom_bar() +
    labs(title = "Number of Genes with Significant Coefficients by Chromosome",
         x = "Chromosome",
         y = "Count") +
    theme_minimal()
}

############################################################################
############################################################################

.render_histo <- function(reactive_histo){
  renderPlotly({
    histo_data <- reactive_histo()
    histo_plot <- ggplotly(.build_histo(histo_data = histo_data))
    return(histo_plot)
    })
}

############################################################################
############################################################################

.build_histo_TF <- function(histo_data){
  plot_ly(histo_data,
          x = ~Chromosome, y = ~Count, type = 'bar') %>%
    layout(title = "Distribution of target genes by chromosomes",
           xaxis = list(title = "Chromosome"),
           yaxis = list(title = "Number of target genes"))
}

############################################################################
############################################################################

.build_histo_TFbyChr <- function(histo_data){
  plot_ly(histo_data,
          x = ~TF, y = ~Count, type = 'bar', color = ~Chromosome) %>%
    layout(title = "Number of genes targeted by TFs",
           xaxis = list(title = "TF"),
           yaxis = list(title = "Number of target genes"),
           barmode = 'group')
}

############################################################################
############################################################################

.render_histo_TF <- function(reactive_histo,
                             by_chr=F){
  renderPlotly({
    histo_data <- reactive_histo()
    if(by_chr){histo_plot <- .build_histo_TFbyChr(histo_data = histo_data)
    }else{histo_plot <- .build_histo_TF(histo_data = histo_data)}
    return(histo_plot)
  })
}

############################################################################
############################################################################

.build_table <- function(table_data){
  DT::datatable(table_data,
                options = list(orderClasses = TRUE))
}

############################################################################
############################################################################

.render_table <- function(reactive_table){
  DT::renderDataTable({
    table_data <- reactive_table()
    ttable <- .build_table(table_data)
    return(ttable)
  })


}

############################################################################
############################################################################

.render_circos <- function(circos_reactive){

  renderUI({
    req(arranged_view_circos())
  })
}

############################################################################
############################################################################

run_shiny <- function(multiomics_integration){
  data <- extract_model_res(multiomics_integration)
  data <- .shiny_preprocess(data)
  data_table <- data$data_table
  ui <- .create_ui(data_table)
  server <- function(input, output, session) {

    # ---------------------- NETWORK SERVER -----------------------------
    nnet <- gINTomics:::.prepare_network(data_table)
    reactive_network <- gINTomics:::.select_deg_network(data_table = data_table,
                                                        input = input,
                                                        output = output,
                                                        network_data = nnet)
    gINTomics:::.render_reactive_network(reactive_network = reactive_network,
                                         input = input,
                                         output = output)


    ### ------------------------ VENN SERVER ----------------------
    reactive_venn <- gINTomics:::.prepare_reactive_venn(data_table = data_table, input = input, output = output)
    output$venn_plot <- gINTomics:::.render_venn(reactive_venn)
    output$common_genes_table <- gINTomics:::.render_venn_table(reactive_venn)

    ## -------------------------- VOLCANO SERVER ------------------------
    reactive_volcano <- gINTomics:::.prepare_reactive_volcano(data_table, input = input, output = output)
    output$volcanoPlot <- gINTomics:::.render_volcano(reactive_volcano)

    ## -------------------------- HEATMAP SERVER ------------------------
    gINTomics:::.prepare_reactive_heatmap(data_table=data_table, multiomics_integration = multiomics_integration, input=input, output=output, session = session)

    ## ---------------------- RIDGE SERVER ------------------------
    reactive_ridge <- gINTomics:::.prepare_reactive_ridge(data_table, input = input, output = output)
    output$ridgelinePlot <- gINTomics:::.render_ridge(reactive_ridge)
    ## ----------------------- HISTO SERVER --------------------------
    reactive_histo <- gINTomics:::.prepare_reactive_histo(data_table, input = input, output = output)
    output$histogramPlot <- gINTomics:::.render_histo(reactive_histo)

    ## ----------------------- HISTO SERVER TF --------------------------
    reactive_histo_tf <- gINTomics:::.prepare_reactive_histo_tf(data_table, input = input, output = output)
    output$histogramPlotTFs <- gINTomics:::.render_histo_TF(reactive_histo_tf, by_chr = F)
    output$histogramPlotTFsByChromosome <- gINTomics:::.render_histo_TF(reactive_histo_tf, by_chr = T)
    #### ------------------- TABLE SERVER ----------------------------
    reactive_table <- gINTomics:::.prepare_reactive_table(data_table, input = input, output = output)
    output$res_table <- gINTomics:::.render_table(reactive_table)
    #### ------------------- ENRICHMENT SERVER ----------------------------
    data_gen_enrich <- data_table[data_table$omics=="gene_genomic_res",]
    callModule(gINTomics:::.background_srv, id = "prova", data_gen_enrich=data_gen_enrich)
  }


  shiny::shinyApp(ui, server)
}


