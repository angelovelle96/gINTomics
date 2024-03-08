.prepare_network <- function(data_table){
  tmp <- intersect(c('cov', 'response', 'pval', 'fdr', 'coef', 'class', 'omics'), colnames(data_table))
  all <- subset(data_table, omics %in% c('tf_res',
                                         'tf_mirna_res',
                                         'mirna_target_res'),
                select = tmp)
  nodes <- data.frame(gene = c(all$cov, all$response))
  edges <- subset(all, select = tmp)
  nodes <- unique(nodes)
  names(nodes) <- 'id'
  names(edges) <- gsub("cov", "from", names(edges))
  names(edges) <- gsub("response", "to", names(edges))
  tf_list <- all[all$omics=='tf_res', 'cov']
  tf_list2 <- all[all$omics=='tf_mirna_res', 'cov']
  mirna_list <- all[all$omics=='tf_mirna_res', 'response']
  mirna_list2 <- all[all$omics=='mirna_target_res', 'cov']
  target_list <- all[all$omics=='tf_res', 'response']
  target_list2 <- all[all$omics=='tf_mirna_res', 'cov']
  tf_list <- unique(c(tf_list, tf_list2))
  mirna_list <- unique(c(mirna_list, mirna_list2))
  target_list <- unique(c(target_list, target_list2))
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
  edges$width <- abs(edges$coef) * 10
  edges$width[edges$width>40] <- 40
  edges$color <- ifelse(edges$coef > 0,
                        '#4169E1',
                        ifelse(edges$coef < 0,
                               '#ED5564',
                               'black'))
  edges$length <- 500
  edges$pval <- round(edges$pval, 3)
  edges$fdr <- round(edges$fdr, 3)
  edges$coef <- round(edges$coef, 3)
  edges$title <- paste('pval:',
                       edges$pval,
                       ',coef:', edges$coef,
                       ',FDR:', edges$fdr)

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
  edges <- edges[order(abs(edges$coef), decreasing = TRUE),]
  return(list(nodes=nodes,
              edges=edges,
              legend_edges=legend_edges,
              legend_nodes=legend_nodes))
}

#########################################################################
#########################################################################
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visGroups
#' @importFrom visNetwork visLegend
#' @importFrom visNetwork visEdges
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visLayout

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
    visLayout(randomSeed = 20, improvedLayout = TRUE)%>%
    visExport(type = "pdf", name = "export-network",
              float = "left", label = "Save network")
}

#########################################################################
#########################################################################
#' @importFrom ggvenn ggvenn
#' @importFrom RColorBrewer brewer.pal

.build_venn <- function(venn_data){
  cnv_sign_genes <- venn_data$cnv_sign_genes
  met_sign_genes <- venn_data$met_sign_genes
  list_venn <- list(CNV = cnv_sign_genes$cov,
                    MET = met_sign_genes$cov)
  ccol <- brewer.pal(name = "Set3", n = max(3,length(list_venn)))
  venn_diagram <- ggvenn(list_venn, fill_color = ccol)
  return(venn_diagram)
}
#########################################################################
#########################################################################
#' @importFrom plotly ggplotly

.render_venn <- function(reactive_venn){
  renderPlotly({
    venn_data <- reactive_venn()
    if(!is.null(venn_data)){
    venn_diagram <- .build_venn(venn_data)
    venn_diagram <- ggplotly(venn_diagram)
    return(venn_diagram)}else{
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "No data available")
    }
  })
}

##########################################################################
##########################################################################

.render_venn_table <- function(reactive_venn) {
  renderDataTable({
    venn_data <- reactive_venn()

    if (!is.null(venn_data$cnv_sign_genes) && !is.null(venn_data$met_sign_genes)){
      cnv_genes <- unlist(venn_data$cnv_sign_genes)
      met_genes <- unlist(venn_data$met_sign_genes)

      if(length(cnv_genes) > 0 & length(met_genes) > 0){
        common_genes <- base::intersect(cnv_genes, met_genes)
        data.frame(Genes = common_genes)
      }else if(length(cnv_genes) == 1 & length(met_genes) == 1){
        data.frame(CNV_Genes = cnv_genes, MET_Genes = met_genes, stringsAsFactors = FALSE)
      }else{
        data.frame(Genes = character(), stringsAsFactors = FALSE)
      }
    }else if(!is.null(venn_data$cnv_sign_genes)){
      data.frame(CNV_Genes = unlist(venn_data$cnv_sign_genes), stringsAsFactors = FALSE)
    }else if(!is.null(venn_data$met_sign_genes)){
      data.frame(MET_Genes = unlist(venn_data$met_sign_genes), stringsAsFactors = FALSE)
    }else{
      data.frame(Genes = character(), stringsAsFactors = FALSE)
    }
  })
}

############################################################################
############################################################################
#' @importFrom plotly plot_ly

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
                         "Pval/FDR:", pval_fdr, "<br>",
                         "coef", coef),
          textposition = 'top right') %>%
    layout(title = "Volcano Plot") %>%
    layout(width = 1000,
           height = 700) #%>%
  #add_text(data = volcano_data, text=~cov)

}
#######################################################################
#######################################################################

.render_volcano <- function(reactive_volcano, annotations){
  renderPlotly({
    volcano_data <- reactive_volcano()
    if(!is.null(volcano_data)){
    volcano_plot <- .build_volcano(volcano_data)
    return(volcano_plot)}else{
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "No data available")
    }
  })
}

############################################################################
############################################################################
#' @importFrom ggplot2 ggplot

.build_ridge <- function(ridge_data,
                         quantiles){

  ggplot(ridge_data,
         aes(x = coef,
             y = significance,
             fill = factor(significance))) +
    geom_density_ridges(jittered_points = TRUE,
                        quantile_lines = TRUE,
                        vline_width = 1,
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
#' @importFrom shiny renderPlot

.render_ridge <- function(reactive_ridge) {
  renderPlot({
    tmp <- reactive_ridge()
    if (!is.null(tmp)) {
      ridge_data <- tmp$df
      quantiles <- tmp$quantiles
      ridge_plot <- .build_ridge(ridge_data = ridge_data,
                                 quantiles = quantiles)
      return(ridge_plot)
    } else {
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "No data available")
    }
  })
}



############################################################################
############################################################################
#' @importFrom ggplot2 ggplot

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
#' @importFrom plotly ggplotly

.render_histo <- function(reactive_histo){
  renderPlotly({
    histo_data <- reactive_histo()
    if(!is.null(histo_data)){
    histo_plot <- ggplotly(.build_histo(histo_data = histo_data))
    return(histo_plot)}else{
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "No data available")
    }
  })
}

############################################################################
############################################################################
#' @importFrom DT renderDataTable

.render_histo_table <- function(reactive_histo_table){
  renderDataTable({
    table_data <- reactive_histo_table()
    if(!is.null(table_data)){
    ttable <- .build_table(table_data)
    return(ttable)}else{
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "No data available")
    }
  })
}

############################################################################
############################################################################
#' @importFrom plotly plot_ly

.build_histo_TFbyChr <- function(histo_data){
  if(nrow(histo_data)==0) return(NULL)
  plot_ly(histo_data,
          x = ~reorder(TF, -Count), y = ~Count, type = 'bar', color = ~Chromosome) %>%
    layout(title = "Number of genes targeted by TFs/miRNAs",
           xaxis = list(title = "TF/miRNAs"),
           yaxis = list(title = "Number of targets"),
           barmode = 'group')
}

############################################################################
############################################################################

.render_histo_TF <- function(reactive_histo){
  renderPlotly({
    histo_data <- reactive_histo()
    if(!is.null(histo_data)){
    histo_plot <- .build_histo_TFbyChr(histo_data=histo_data)
    return(histo_plot)}else{
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), main = "No data available")
    }
  })
}

#' @importFrom DT renderDataTable

.render_histo_tf_table <- function(reactive_histo_tf_table){
  renderDataTable({
    table_data <- reactive_histo_tf_table()
    if(!is.null(table_data)){
    ttable <- .build_table(table_data)
    return(ttable)}else{
      return(NULL)
    }
  })
}
############################################################################
############################################################################
#' @importFrom DT renderDataTable

.render_ridge_table <- function(reactive_ridge_table){
  renderDataTable({
    table_data <- reactive_ridge_table()
    if(!is.null(table_data)){
    ttable <- .build_table(table_data)
    return(ttable)}else{
      return(NULL)
    }
  })
}
############################################################################
############################################################################
#' @importFrom DT datatable

.build_table <- function(table_data){
  datatable(table_data,
            options = list(orderClasses = TRUE))
}
############################################################################
############################################################################
#' @importFrom DT renderDataTable

.render_table <- function(reactive_table){
  renderDataTable({
    table_data <- reactive_table()
    if(!is.null(table_data)){
    ttable <- .build_table(table_data)
    return(ttable)}else{
      return(NULL)
    }
  })
}
############################################################################
############################################################################

download_csv <- function(type = "genomic",
                         deg = FALSE,
                         plotType = "histo",
                         input,
                         output,
                         data_table){
  if(plotType == "histo"){
    handler <- downloadHandler(
      filename = function(){
        "data_table_histo.csv"
      },
      content = function(file){
        reactive_histo_table <- gINTomics:::.prepare_reactive_histo_table(type=type,
                                                                          deg=deg,
                                                                          data_table=data_table,
                                                                          input=input,
                                                                          output=output)
        write.csv(reactive_histo_table(), file, row.names = FALSE)
      }
    )
  }
  if(plotType == "ridge"){
    handler <- downloadHandler(
      filename = function(){
        "data_ridge.csv"
      },
      content = function(file){
        reactive_ridge_table <- gINTomics:::.prepare_reactive_ridge_table(type=type,
                                                                          deg=deg,
                                                                          data_table=data_table,
                                                                          input=input,
                                                                          output=output)
        write.csv(reactive_ridge_table(), file, row.names = FALSE)
      }
    )
  }
  if(plotType == "table"){
    handler <- downloadHandler(
      filename = function(){
        "data_table.csv"
      },
      content = function(file){
        reactive_table <- gINTomics:::.prepare_reactive_table(data_table=data_table,
                                                              input=input,
                                                              output=output)
        write.csv(reactive_table(), file, row.names = FALSE)
      }
    )
  }
  if(plotType == "dotPlot"){
    handler <- downloadHandler(
      filename = function(){
        "enrich_table.csv"
      },
      content = function(file){
        reactive_enr <- gINTomics:::.prepare_reactive_table(deg=deg,
                                                            data_table=data_table,
                                                            input=input,
                                                            output=output)
        write.csv(reactive_enr(), file, row.names = FALSE)
      }
    )
  }
  if(plotType == "venn"){
    handler <- downloadHandler(
      filename = function(){
        "venn_table.csv"
      },
      content = function(file){
        reactive_venn <- gINTomics:::.prepare_reactive_venn(deg=deg,
                                                            data_table=data_table,
                                                            input=input,
                                                            output=output)
        venn_data <- reactive_venn()
        cnv_genes <- unlist(venn_data$cnv_sign_genes)
        met_genes <- unlist(venn_data$met_sign_genes)
        common_genes <- base::intersect(cnv_genes, met_genes)
        venn_data <- data.frame(Genes = common_genes)
        write.csv(venn_data, file, row.names = FALSE)
      }
    )
  }
  return(handler)
}
############################################################################
############################################################################
#' @importFrom shiny renderUI

.render_circos <- function(circos_reactive){

  renderUI({
    req(arranged_view_circos())
  })
}
##########################################################################
##########################################################################
.circos_preprocess <- function(data){

  library(GenomicRanges)
  dataframes <- lapply(unique(data$omics), function(x) {
    single_omic_df <- data[data$omics==x,]
    if(max(single_omic_df$response_value, na.rm = T)>30){
      single_omic_df$response_value <- log2(single_omic_df$response_value+1)
    }
    return(single_omic_df)
  })
  names(dataframes) <- paste0("df_", unique(data$omics))
  if("df_gene_genomic_res"%in%names(dataframes)){
    dataframes$df_cnv <- filter(dataframes$df_gene_genomic_res,
                                cnv_met == 'cnv')
    dataframes$df_met <- filter(dataframes$df_gene_genomic_res,
                                cnv_met == 'met')
    tmp <- lapply(which(names(dataframes)!="df_gene_genomic_res"), function(x){
      ans <- dataframes[[x]]
      ans <- ans[ans$cov!="(Intercept)",]
      ans$cov_value_original <- ans$cov_value
      ans$cov_value <- abs(ans$cov_value)
      ans$coef_original <- ans$coef
      ans$coef <- abs(ans$coef)
      return(ans)
    })
    names(tmp) <- names(dataframes)[
      which(names(dataframes)!="df_gene_genomic_res")]
    dataframes <- tmp
  }
  gr <- lapply(dataframes, function(x){
    ans <- makeGRangesFromDataFrame(x,
                                    seqnames.field = 'chr',
                                    start.field = 'start',
                                    end.field = 'end',
                                    keep.extra.columns = TRUE,
                                    na.rm = TRUE)
    return(ans)
  })
  return(gr)
}


############################################################################
#############################################################################
#' @importFrom RColorBrewer brewer.pal
.circos_track_cols <- function(vvalues,
                               n=50,
                               colors=NULL){

if(is.null(colors)){
  colors <- rev(brewer.pal(name = "RdBu", n = 11))
}
vvalues <- as.numeric(vvalues)
if(min(vvalues)<0){
  breaks <- c(seq(min(vvalues, na.rm = T),0, len=n/2),
               seq(0.00001, max(vvalues, na.rm = T), len=n/2))
}else{
  breaks <- c(seq(min(vvalues, na.rm = T),mean(vvalues, na.rm = T), len=n/2),
               seq(mean(vvalues, na.rm = T)+0.00001,
                   max(vvalues, na.rm = T), len=n/2))
}
ii <- cut(unique(vvalues),
          breaks = breaks,
          include.lowest = TRUE)
ccol <- colorRampPalette(colors)(n-1)[ii]
return(ccol)

}


############################################################################
#############################################################################
.create_single_track <- function(data,
                                 dataValue=NULL,
                                 x_axis=NULL,
                                 xe_axis=NULL,
                                 y_axis=NULL,
                                 colorField=NULL,
                                 colorDomain=NULL,
                                 colorRange=NULL,
                                 tooltipField1=NULL,
                                 tooltipTitle=NULL,
                                 tooltipAlt1=NULL,
                                 tooltipField2=NULL,
                                 tooltipAlt2=NULL,
                                 tooltipField3=NULL,
                                 tooltipAlt3=NULL,
                                 tooltipField4=NULL,
                                 tooltipAlt4=NULL,
                                 tooltipField5=NULL,
                                 tooltipAlt5=NULL,
                                 legend=NULL,
                                 colorType=NULL,
                                 title=NULL) {
  return(
    add_single_track(
      data = track_data_gr(data, chromosomeField = 'seqnames', genomicFields = c('start','end'), value = dataValue),
      mark = 'bar',
      x = visual_channel_x(field = 'start', type = 'genomic', axis = x_axis),
      xe = visual_channel_x(field = 'end', type = 'genomic', axis = xe_axis),
      y = visual_channel_y(field = dataValue, type = 'quantitative', axis = y_axis),
      color = visual_channel_color(field = colorField, type = colorType, domain = colorDomain, range = colorRange, legend = legend),
      tooltip = visual_channel_tooltips(
        visual_channel_tooltip(field = "start", type = "genomic", alt = 'Start Position:'),
        visual_channel_tooltip(field = "end", type = "genomic", alt = "End Position:"),
        visual_channel_tooltip(field = tooltipField1, title = tooltipTitle, type = "quantitative", alt = paste(tooltipTitle, "Value:"), format = "0.2"),
        visual_channel_tooltip(field = tooltipField2, type = 'nominal', alt = tooltipAlt2),
        visual_channel_tooltip(field = tooltipField3, type = 'nominal', alt = tooltipAlt3),
        visual_channel_tooltip(field = tooltipField4, type = 'nominal', alt = tooltipAlt4),
        visual_channel_tooltip(field = tooltipField5, type = 'nominal', alt = tooltipAlt5)
      ),
      size = list(value = 1)
      , title=title)
  )
}
#######################################################################
########################################################################

.create_tracks <- function(data, gr){

  tracks <- list()
  track_cyto <- .create_cyto_track()
  if ("gene_genomic_res" %in% unique(data$omics)){
    gr$df_cnv$cnv_met2 <- "cnv/met"
    track_cnv <- .create_cnv_track(gr$df_cnv)
    track_met <- .create_met_track(gr$df_met)
    track_expr <- .create_expr_track(gr$df_cnv, genomic=T)
    track_coef_cnv <- .create_coef_track(gr$df_cnv)
    track_coef_met <- .create_coef_track(gr$df_met)
    tracks <- c(tracks, list(track_cnv=track_cnv,
                             track_met=track_met,
                             track_expr=track_expr,
                             track_coef_cnv=track_coef_cnv,
                             track_coef_met=track_coef_met))
  }

  if ("df_gene_cnv_res" %in% unique(names(gr))){
    gr$df_gene_cnv_res$cnv_met <- "cnv"
    track_cnv <- .create_cnv_track(gr$df_gene_cnv_res)
    track_expr <- .create_expr_track(gr$df_gene_cnv_res)
    track_coef_cnv <- .create_coef_track(gr$df_gene_cnv_res)
    tracks <- c(tracks,list(track_cnv=track_cnv,
                            track_expr=track_expr,
                            track_coef_cnv=track_coef_cnv))
  }

  if("df_gene_met_res"%in%unique(names(gr))){
    gr$df_gene_met_res$cnv_met <- "met"
    track_met <- .create_met_track(gr$df_gene_met_res)
    track_expr <- .create_expr_track(gr$df_gene_met_res)
    track_coef_met <- .create_coef_track(gr$df_gene_met_res)
    tracks <- c(tracks, list(track_met=track_met,
                             track_expr=track_expr,
                             track_coef_met=track_coef_met))
  }

  if ("df_mirna_cnv_res" %in% unique(names(gr))){
    gr$df_mirna_cnv_res$cnv_met <- "cnv"
    track_mirna_cnv <- .create_cnv_track(gr$df_mirna_cnv_res)
    track_mirna_expr <- .create_expr_track(gr$df_mirna_cnv_res)
    track_mirna_coef_cnv <- .create_coef_track(gr$df_mirna_cnv_res)
    tracks <- c(tracks,list(track_mirna_cnv=track_mirna_cnv,
                            track_mirna_expr=track_mirna_expr,
                            track_mirna_coef_cnv=track_mirna_coef_cnv))
  }
  tracks <- c(tracks,list(track_cyto=track_cyto))
  return(tracks)
}
#######################################################################
########################################################################

.create_cnv_track <- function(gr){

  gr$cov_value2 <- as.character(gr$cov_value_original)
  ccol <- .circos_track_cols(vvalues = gr$cov_value2)
  track_cnv <- .create_single_track(data=gr,
                                    dataValue='cov_value',
                                    x_axis="none",
                                    xe_axis="none",
                                    y_axis="none",
                                    colorField="cov_value2",
                                    colorDomain=unique(gr$cov_value2),
                                    colorRange=ccol,
                                    tooltipField1="cov_value",
                                    tooltipTitle="cnv",
                                    tooltipAlt1="CNV Value:",
                                    tooltipField2="gene",
                                    tooltipAlt2="Gene Name:",
                                    tooltipField3="class",
                                    tooltipAlt3="Class:",
                                    tooltipField4="cnv_met",
                                    tooltipAlt4="Integration Type:",
                                    tooltipField5="cov_value_original",
                                    tooltipAlt5="Original CNV Value:",
                                    legend=FALSE,
                                    colorType="nominal",
                                    title="CNV")
  return(track_cnv)
}
#######################################################################
########################################################################

.create_met_track <- function(gr){
  gr$cov_value2 <- as.character(gr$cov_value_original)
  ccol <- .circos_track_cols(vvalues = gr$cov_value2)
  track_met <- .create_single_track(data=gr,
                                    dataValue='cov_value',
                                    x_axis="none",
                                    xe_axis="none",
                                    y_axis="none",
                                    colorField="cov_value2",
                                    colorDomain=unique(gr$cov_value2),
                                    colorRange=ccol,
                                    tooltipField1="cov_value",
                                    tooltipTitle="met",
                                    tooltipAlt1="MET Value:",
                                    tooltipField2="gene",
                                    tooltipAlt2="Gene Name:",
                                    tooltipField3="class",
                                    tooltipAlt3="Class:",
                                    tooltipField4="cnv_met",
                                    tooltipAlt4="Integration Type:",
                                    tooltipField5="cov_value_original",
                                    tooltipAlt5="Original MET Value:",
                                    legend=FALSE,
                                    colorType="nominal",
                                    title="MET")
  return(track_met)

}
#######################################################################
########################################################################

.create_expr_track <- function(gr, genomic=F){
  cnv_met <- ifelse(genomic, "cnv_met2", "cnv_met")
  track_expr <- .create_single_track(data=gr,
                                     dataValue='response_value',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="response_value",
                                     colorDomain=NULL,
                                     colorRange=NULL,
                                     tooltipField1="response_value",
                                     tooltipTitle="expr",
                                     tooltipAlt1="Expression Value (log2):",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4=cnv_met,
                                     tooltipAlt4="Integration Type:",
                                     legend=T,
                                     colorType="quantitative",
                                     title="Expression")
  return(track_expr)
}
#######################################################################
########################################################################

.create_coef_track <- function(gr){
  gr$coef2 <- as.character(gr$coef_original)
  ccol <- .circos_track_cols(vvalues = gr$coef2)
  track_coef <- .create_single_track(data=gr,
                                     dataValue='coef',
                                     x_axis="none",
                                     xe_axis="none",
                                     y_axis="none",
                                     colorField="coef2",
                                     colorDomain=unique(gr$coef2),
                                     colorRange=ccol,
                                     tooltipField1="coef",
                                     tooltipTitle="coef",
                                     tooltipAlt1="Coef Value:",
                                     tooltipField2="gene",
                                     tooltipAlt2="Gene Name:",
                                     tooltipField3="class",
                                     tooltipAlt3="Class:",
                                     tooltipField4="cnv_met",
                                     tooltipAlt4="Integration Type:",
                                     tooltipField5="coef_original",
                                     tooltipAlt5="Original Coef Value:",
                                     legend=FALSE,
                                     colorType="nominal",
                                     title="Coefficients")
  return(track_coef)
}
#######################################################################
########################################################################

.create_cyto_track <- function(){
  track_cyto <- add_single_track(
    id = "track2",
    data = track_data(
      url = "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
      type = "csv",
      chromosomeField = "Chromosome",
      genomicFields = c("chromStart",
                        "chromEnd")
    ),
    mark = "rect",
    x = visual_channel_x(field = "chromStart",
                         type = "genomic"),
    xe = visual_channel_x(field = "chromEnd",
                          type = "genomic"),
    color = visual_channel_color(
      field = "Stain",
      type = "nominal",
      domain = c(
        "gneg",
        "gpos25",
        "gpos50",
        "gpos75",
        "gpos100",
        "gvar",
        "acen"           # acen: centromeric region (UCSC band files)
      ),
      range = c(
        "white",
        "#D9D9D9",
        "#979797",
        "#636363",
        "black",
        "#F0F0F0",
        "red"
      )
    ),
    stroke = visual_channel_stroke(
      value = "lightgray"
    ),
    strokeWidth = visual_channel_stroke_width(
      value = 0.5
    ),
  )
  return(track_cyto)
}
#######################################################################
########################################################################

.create_composed_view <- function(tracks,width, height) {

  composed_views <- list()

  if(sum(c("track_expr", "track_cnv", "track_met")%in%names(tracks))==3){

    composed_view_circos_genomic <- compose_view(
      multi = TRUE,
      tracks = add_multi_tracks(tracks$track_cyto,
                                tracks$track_expr,
                                tracks$track_cnv,
                                tracks$track_coef_cnv,
                                tracks$track_met,
                                tracks$track_coef_met),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail",
      width = width,
      height = height
    )
    composed_views <- c(composed_views, list(circos_genomic=composed_view_circos_genomic))
  }else{

    if("track_met"%in%names(tracks)){

      composed_view_circos_met_gene <- compose_view(
        multi = TRUE,
        tracks = add_multi_tracks(tracks$track_cyto,
                                  tracks$track_expr,
                                  tracks$track_met,
                                  tracks$track_coef_met),
        alignment = 'stack',
        spacing = 0.01,
        linkingId = "detail"
      )
      composed_views <- c(composed_views, list(circos_met_gene=composed_view_circos_met_gene))
    }
    if("track_cnv"%in%names(tracks)){

      composed_view_circos_cnv_gene <- compose_view(
        multi = TRUE,
        tracks = add_multi_tracks(tracks$track_cyto,
                                  tracks$track_expr,
                                  tracks$track_cnv,
                                  tracks$track_coef_cnv),
        alignment = 'stack',
        spacing = 0.01,
        linkingId = "detail"
      )
      composed_views <- c(composed_views, list(circos_cnv_gene=composed_view_circos_cnv_gene))
    }

  }
  if(sum(c("track_mirna_expr", "track_mirna_cnv")%in%names(tracks))==2){

    composed_view_circos_genomic_mirna <- compose_view(
      multi = TRUE,
      tracks = add_multi_tracks(tracks$track_cyto,
                                tracks$track_mirna_expr,
                                tracks$track_mirna_cnv,
                                tracks$track_mirna_coef_cnv),
      alignment = 'stack',
      spacing = 0.01,
      linkingId = "detail",
      width = width,
      height = height
    )
    composed_views <- c(composed_views,
                        list(circos_genomic_mirna=composed_view_circos_genomic_mirna))
  }
  return(composed_views)
}

#' @importFrom ComplexHeatmap Heatmap
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap rowAnnotation

.prepare_gen_heatmap <- function(data_table,
                                 df_heatmap,
                                 df_heatmap_t,
                                 significativityCriteria,
                                 pvalRange,
                                 fdrRange,
                                 numTopCNV,
                                 numTopMET,
                                 scale){

  tmp <- data.frame(cnv = data_table$coef[
    data_table$cnv_met == 'cnv'],
    pval_cnv = data_table$pval[
      data_table$cnv_met == 'cnv'],
    fdr_cnv = data_table$fdr[
      data_table$cnv_met == 'cnv'],
    row.names = data_table$response[
      data_table$cnv_met == 'cnv'])

  tmp2 <- data.frame(met = data_table$coef[
    data_table$cnv_met == 'met'],
    pval_met = data_table$pval[
      data_table$cnv_met == 'met'],
    fdr_met = data_table$fdr[
      data_table$cnv_met == 'met'],
    row.names = data_table$response[
      data_table$cnv_met == 'met'])
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
  expr_top <- rbind(top_cnv, top_met[
    !rownames(top_met)%in%rownames(top_cnv),])
  expr_top_subset <- expr_top[, -c((ncol(expr_top) - 5):ncol(expr_top))]
  set.seed(123)
  row_ha <- rowAnnotation(coef_cnv = expr_top$cnv, coef_met = expr_top$met)
  expr_top_subset <- as.matrix(expr_top_subset)
  if(scale=="row") expr_top_subset <- t(scale(t(log2(expr_top_subset+1))))
  if(scale=="col") expr_top_subset <- scale(log2(expr_top_subset+1))
  if(scale=="none") expr_top_subset <- log2(expr_top_subset+1)
  ht <- Heatmap(expr_top_subset,
                right_annotation = row_ha,
                heatmap_legend_param = list(title="log2 expression"))
  ht <- draw(ht)
  return(ht)
}
######################################################################
######################################################################
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap rowAnnotation

.prepare_cnv_heatmap <- function(data_table,
                                 df_heatmap,
                                 df_heatmap_t,
                                 significativityCriteria,
                                 pvalRange,
                                 fdrRange,
                                 numTopCNVonly,
                                 scale){

  tmp <-  data.frame(cnv=data_table$coef[
    data_table$omics == 'gene_cnv_res'],
    pval_cnv=data_table$pval[
      data_table$omics == 'gene_cnv_res'],
    fdr_cnv=data_table$fdr[
      data_table$omics == 'gene_cnv_res'],
    row.names = data_table$response[
      data_table$omics == 'gene_cnv_res'])
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
  df_heatmap_t <- as.data.frame(df_heatmap_t)
  top_cnv <- df_heatmap_t %>%
    arrange(desc(abs(cnv))) %>%
    head(numTopCNVonly)
  expr_top <- top_cnv
  expr_top_subset <- expr_top[, -c((ncol(expr_top) - 2):ncol(expr_top))]
  set.seed(123)
  row_ha <- rowAnnotation(coef_cnv=expr_top$cnv)
  expr_top_subset <- as.matrix(expr_top_subset)
  if(scale=="row") expr_top_subset <- t(scale(t(log2(expr_top_subset+1))))
  if(scale=="col") expr_top_subset <- scale(log2(expr_top_subset+1))
  if(scale=="none") expr_top_subset <- log2(expr_top_subset+1)
  ht <- Heatmap(expr_top_subset,
                right_annotation = row_ha,
                heatmap_legend_param = list(title="log2 expression"))
  ht <- draw(ht)
  return(ht)
}

#' @importFrom ComplexHeatmap Heatmap
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap rowAnnotation

.prepare_met_heatmap <- function(data_table,
                                 df_heatmap,
                                 df_heatmap_t,
                                 significativityCriteria,
                                 pvalRange,
                                 fdrRange,
                                 numTopMETonly,
                                 scale){

  tmp <-  data.frame(met=data_table$coef[
    data_table$omics == 'gene_met_res'],
    pval_met=data_table$pval[
      data_table$omics == 'gene_met_res'],
    fdr_met=data_table$fdr[
      data_table$omics == 'gene_met_res'],
    row.names = data_table$response[
      data_table$omics == 'gene_met_res'])
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
  df_heatmap_t <- as.data.frame(df_heatmap_t)
  top_met <- df_heatmap_t %>%
    arrange(desc(abs(met))) %>%
    head(numTopMETonly)
  expr_top <- top_met
  expr_top_subset <- expr_top[, -c((ncol(expr_top) - 2):ncol(expr_top))]
  set.seed(123)
  row_ha <- rowAnnotation(coef_met=expr_top$met)
  expr_top_subset <- as.matrix(expr_top_subset)
  if(scale=="row") expr_top_subset <- t(scale(t(log2(expr_top_subset+1))))
  if(scale=="col") expr_top_subset <- scale(log2(expr_top_subset+1))
  if(scale=="none") expr_top_subset <- log2(expr_top_subset+1)
  ht <- Heatmap(expr_top_subset,
                right_annotation = row_ha,
                heatmap_legend_param = list(title="log2 expression"))
  ht <- draw(ht)
  return(ht)
}

#' @importFrom ComplexHeatmap Heatmap
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap rowAnnotation

.prepare_mirna_heatmap <- function(data_table,
                                   df_heatmap,
                                   df_heatmap_t,
                                   significativityCriteria,
                                   pvalRange,
                                   fdrRange,
                                   numTopMiCNV,
                                   scale){

  tmp <-  data.frame(mirna_cnv=data_table$coef[
    data_table$omics == 'mirna_cnv_res'],
    pval_mirna_cnv=data_table$pval[
      data_table$omics == 'mirna_cnv_res'],
    fdr_mirna_cnv=data_table$fdr[
      data_table$omics == 'mirna_cnv_res'],
    row.names = data_table$response[
      data_table$omics == 'mirna_cnv_res'])
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
  expr_top_subset <- as.matrix(expr_top)
  if(scale=="row") expr_top_subset <- t(scale(t(log2(expr_top_subset+1))))
  if(scale=="col") expr_top_subset <- scale(log2(expr_top_subset+1))
  if(scale=="none") expr_top_subset <- log2(expr_top_subset+1)
  ht <- Heatmap(expr_top_subset,
                right_annotation = row_ha,
                heatmap_legend_param = list(title="log2 expression"))
  ht <- draw(ht)
  return(ht)
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
    reactive_network <- gINTomics:::.select_network(data_table = data_table,
                                                        input = input,
                                                        output = output,
                                                        network_data = nnet,
                                                        deg = FALSE)
    gINTomics:::.render_reactive_network(reactive_network = reactive_network,
                                         input = input,
                                         output = output,
                                         deg = FALSE)

    reactive_network_deg <- gINTomics:::.select_network(data_table = data_table,
                                                            input = input,
                                                            output = output,
                                                            network_data = nnet,
                                                            deg = TRUE)
    gINTomics:::.render_reactive_network(reactive_network = reactive_network_deg,
                                         input = input,
                                         output = output,
                                         deg = TRUE)
    ### ------------------------ VENN SERVER ----------------------
    reactive_venn <- gINTomics:::.prepare_reactive_venn(data_table = data_table,
                                                        input = input,
                                                        output = output,
                                                        deg = FALSE)
    output$venn_plot<- gINTomics:::.render_venn(reactive_venn)
    output$common_genes_table <- gINTomics:::.render_venn_table(reactive_venn)
    reactive_venn_deg <- gINTomics:::.prepare_reactive_venn(data_table = data_table,
                                                            input = input,
                                                            output = output,
                                                            deg = TRUE)
    output$venn_plotDEG<- gINTomics:::.render_venn(reactive_venn_deg)
    output$venn_tableDEG <- gINTomics:::.render_venn_table(reactive_venn_deg)
    output$download_csv_venn_gen <- download_csv(deg = FALSE,
                                             plotType = "venn",
                                             input=input,
                                             output=output,
                                             data_table=data_table)
    output$download_csv_venn_deg <- download_csv(deg = TRUE,
                                                 plotType = "venn",
                                                 input=input,
                                                 output=output,
                                                 data_table=data_table)
    ## -------------------------- VOLCANO SERVER ------------------------
    reactive_volcano <- gINTomics:::.prepare_reactive_volcano(data_table,
                                                              input = input,
                                                              output = output,
                                                              type = "genomic",
                                                              deg = FALSE)
    output$volcanoPlot <- gINTomics:::.render_volcano(reactive_volcano)
    reactive_volcano_transcript <- gINTomics:::.prepare_reactive_volcano(data_table,
                                                                         input = input,
                                                                         output = output,
                                                                         type = "transcript",
                                                                         deg = FALSE)
    output$transcriptVolcanoPlot <- gINTomics:::.render_volcano(reactive_volcano_transcript)
    reactive_volcano_deg <- gINTomics:::.prepare_reactive_volcano(data_table,
                                                                  input = input,
                                                                  output = output,
                                                                  type = "all",
                                                                  deg = TRUE)
    output$volcanoPlotDEG <- gINTomics:::.render_volcano(reactive_volcano_deg)
    ## -------------------------- HEATMAP SERVER ------------------------
    gINTomics:::.prepare_reactive_heatmap(data_table=data_table,
                                          multiomics_integration = multiomics_integration,
                                          input=input,
                                          output=output,
                                          session = session,
                                          deg = FALSE)
    gINTomics:::.prepare_reactive_heatmap(data_table=data_table,
                                          multiomics_integration = multiomics_integration,
                                          input=input,
                                          output=output,
                                          session = session,
                                          deg = TRUE)
    ## ---------------------- RIDGE SERVER ------------------------
    reactive_ridge <- gINTomics:::.prepare_reactive_ridge(data_table,
                                                          input = input,
                                                          output = output,
                                                          type = "genomic",
                                                          deg = FALSE)
    output$ridgelinePlot <- gINTomics:::.render_ridge(reactive_ridge)
    reactive_ridge_table <- gINTomics:::.prepare_reactive_ridge_table(data_table,
                                                                      input = input,
                                                                      output = output,
                                                                      type = "genomic",
                                                                      deg = FALSE)
    output$ridgelineTable <- gINTomics:::.render_ridge_table(reactive_ridge_table)
    reactive_ridge_transcript <- gINTomics:::.prepare_reactive_ridge(data_table,
                                                                     input = input,
                                                                     output = output,
                                                                     type = "transcript",
                                                                     deg = FALSE)
    output$ridgelinePlotTranscript <- gINTomics:::.render_ridge(reactive_ridge_transcript)
    reactive_ridge_tableTranscript <- gINTomics:::.prepare_reactive_ridge_table(data_table,
                                                                                input = input,
                                                                                output = output,
                                                                                type = "transcript",
                                                                                deg = FALSE)
    output$ridgelineTableTranscript <- gINTomics:::.render_ridge_table(reactive_ridge_tableTranscript)
    reactive_ridge_deg <- gINTomics:::.prepare_reactive_ridge(data_table,
                                                              input = input,
                                                              output = output,
                                                              type = "all",
                                                              deg = TRUE)
    output$ridgelinePlotDEG <- gINTomics:::.render_ridge(reactive_ridge_deg)
    reactive_ridge_tableDEG <- gINTomics:::.prepare_reactive_ridge_table(data_table,
                                                                         input = input,
                                                                         output = output,
                                                                         type = "all",
                                                                         deg = TRUE)
    output$ridgelineTableDEG <- gINTomics:::.render_ridge_table(reactive_ridge_tableDEG)
    output$download_csv_ridge_gen <- download_csv(type = "genomic",
                                                  deg = FALSE,
                                                  plotType = "ridge",
                                                  input=input,
                                                  output=output,
                                                  data_table=data_table)
    output$download_csv_ridge_transcr <- download_csv(type = "transcript",
                                                      deg = FALSE,
                                                      plotType = "ridge",
                                                      input=input,
                                                      output=output,
                                                      data_table=data_table)
    output$download_csv_ridge_deg <- download_csv(type = "all",
                                                  deg = TRUE,
                                                  plotType = "ridge",
                                                  input=input,
                                                  output=output,
                                                  data_table=data_table)
    ## ----------------------- HISTO SERVER --------------------------
    reactive_histo <- gINTomics:::.prepare_reactive_histo(data_table,
                                                          input = input,
                                                          output = output,
                                                          type = "genomic",
                                                          deg = FALSE)
    output$histogramPlot <- gINTomics:::.render_histo(reactive_histo)
    reactive_histo_table <- gINTomics:::.prepare_reactive_histo_table(data_table,
                                                                      input=input,
                                                                      output=output,
                                                                      type="genomic",
                                                                      deg=FALSE)
    output$histogramTable <- gINTomics:::.render_histo_table(reactive_histo_table)
    reactive_histo_transcript <- gINTomics:::.prepare_reactive_histo(data_table,
                                                                     input=input,
                                                                     output=output,
                                                                     type="transcript",
                                                                     deg=FALSE)
    output$histogramPlotTranscript <- gINTomics:::.render_histo(reactive_histo_transcript)
    reactive_histo_table_transcript <- gINTomics:::.prepare_reactive_histo_table(data_table,
                                                                                 input=input,
                                                                                 output=output,
                                                                                 type="transcript",
                                                                                 deg=FALSE)
    output$histogramTableTranscript <- gINTomics:::.render_histo_table(reactive_histo_table_transcript)
    reactive_histo_deg <- gINTomics:::.prepare_reactive_histo(data_table,
                                                              input=input,
                                                              output=output,
                                                              type="all",
                                                              deg=TRUE)
    output$histogramPlotDEG <- gINTomics:::.render_histo(reactive_histo_deg)
    reactive_histo_table_deg <- gINTomics:::.prepare_reactive_histo_table(data_table,
                                                                          input=input,
                                                                          output=output,
                                                                          type="all",
                                                                          deg=TRUE)
    output$histogramTableDEG <- gINTomics:::.render_histo_table(reactive_histo_table_deg)
    output$download_csv_histo_gen <- download_csv(type="genomic",
                                                  deg=FALSE,
                                                  plotType="histo",
                                                  input=input,
                                                  output=output,
                                                  data_table=data_table)
    output$download_csv_histo_transcr <- download_csv(type="transcript",
                                                      deg=FALSE,
                                                      plotType="histo",
                                                      input=input,
                                                      output=output,
                                                      data_table=data_table)
    output$download_csv_histo_deg <- download_csv(type="all",
                                                  deg=TRUE,
                                                  plotType="histo",
                                                  input=input,
                                                  output=output,
                                                  data_table=data_table)
    ## ----------------------- HISTO SERVER TF --------------------------
    reactive_histo_tf_transcript <- gINTomics:::.prepare_reactive_histo_tf(data_table,
                                                                input = input,
                                                                output = output,
                                                                deg = FALSE)
    output$histogramTfTranscript <- gINTomics:::.render_histo_TF(reactive_histo_tf_transcript)
    reactive_histo_tf_deg <- gINTomics:::.prepare_reactive_histo_tf(data_table,
                                                                input = input,
                                                                output = output,
                                                                deg = TRUE)
    output$histogramTfTranscriptDEG <- gINTomics:::.render_histo_TF(reactive_histo_tf_deg)





    reactive_tf_table <- gINTomics:::.prepare_reactive_histoTF_table(data_table,
                                                                      input=input,
                                                                      output=output,
                                                                      deg=FALSE)
    reactive_tf_table_deg <- gINTomics:::.prepare_reactive_histoTF_table(data_table,
                                                                      input=input,
                                                                      output=output,
                                                                      deg=TRUE)
    output$download_csv_histo_tf <- download_csv(deg=FALSE,
                                                 plotType="histo",
                                                 input=input,
                                                 output=output,
                                                 data_table=data_table)
    output$download_csv_histo_tf_deg <- download_csv(deg=TRUE,
                                                  plotType="histo",
                                                  input=input,
                                                  output=output,
                                                  data_table=data_table)
    output$histogramTfTable <- gINTomics:::.render_histo_tf_table(reactive_tf_table)
    output$histogramTfTableDEG <- gINTomics:::.render_histo_tf_table(reactive_tf_table_deg)
    #### ------------------- TABLE SERVER ----------------------------
    reactive_table <- gINTomics:::.prepare_reactive_table(data_table,
                                                          input=input,
                                                          output=output)
    output$res_table <- gINTomics:::.render_table(reactive_table)
    output$download_csv_table <- download_csv(plotType = "table",
                                              input=input,
                                              output=output,
                                              data_table=data_table)
    #### ------------------- ENRICHMENT SERVER ----------------------------
    data_gen_enrich <- data_table[data_table$omics=="gene_genomic_res",]
    data_tf_enrich <- data_table[data_table$omics=="tf_res",]
    callModule(gINTomics:::.background_srv, id = "prova", name="prova",
               data_gen_enrich=data_gen_enrich,
               data_tf_enrich=data_tf_enrich)
    #### ------------------- CIRCOS SERVER ----------------------------
    reactive_circos <- .prepare_reactive_circos(data = data$data,
                                                input = input,
                                                output = output)
    output$gosling_plot_circos <- renderGosling({
      gosling(component_id = "component_1",
              reactive_circos())
    })
    observe({
      export_pdf(component_id = "component_1")
    })%>%bindEvent(input$circosDownload_pdf)
    observe({
      export_png(component_id = "component_1")
    })%>%bindEvent(input$circosDownload_png)

  }


  shiny::shinyApp(ui, server)
}
