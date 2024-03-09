#' Server Genomic integration chr distribution
.server_histo <- function(input,
                          output,
                          session,
                          data_table){

reactive_histo <- .prepare_reactive_histo(data_table,
                                          input = input,
                                          output = output,
                                          deg = FALSE)
output$plotly <- .render_histo(reactive_histo)
reactive_histo_table <- .prepare_reactive_histo(data_table,
                                                input=input,
                                                output=output,
                                                deg=FALSE,
                                                table=TRUE)
output$table <- .render_histo_table(reactive_histo_table)
output$download_csv <- download_csv(deg=FALSE,
                                    plotType="histo",
                                    input=input,
                                    output=output,
                                    data_table=data_table)

}

#' Server transcriptional and deg integration
.server_histo2 <- function(input,
                           output,
                           session,
                           data_table,
                           deg=FALSE){

  reactive_histo_transcript <- .prepare_reactive_histo(data_table,
                                                       input=input,
                                                       output=output,
                                                       deg=deg)
  output$plotly <- .render_histo(reactive_histo_transcript)
  reactive_histo_table_transcript <- .prepare_reactive_histo(data_table,
                                                             input=input,
                                                             output=output,
                                                             deg=deg,
                                                             table=TRUE)
  output$table <- .render_histo_table(reactive_histo_table_transcript)
  output$download_csv <- download_csv(deg=deg,
                                      plotType="histo",
                                      input=input,
                                      output=output,
                                      data_table=data_table)

  reactive_histo_tf_transcript <- .prepare_reactive_histo_tf(data_table,
                                                             input = input,
                                                             output = output,
                                                             deg = deg)
  output$plotly_tf <- .render_histo_TF(reactive_histo_tf_transcript)
  reactive_tf_table <- .prepare_reactive_histo_tf(data_table,
                                                 input=input,
                                                 output=output,
                                                 deg=deg,
                                                 table=TRUE)
  output$download_csv_tf <- download_csv(deg=deg,
                                         tf=TRUE,
                                         plotType="histo",
                                         input=input,
                                         output=output,
                                         data_table=data_table)
  output$table_tf <- .render_histo_tf_table(reactive_tf_table)

}

#' Server network
.server_network <- function(input,
                            output,
                            session,
                            data_table,
                            deg=FALSE){
  nnet <- .prepare_network(data_table)
  reactive_network <- .select_network(data_table = data_table,
                                    input = input,
                                    output = output,
                                    network_data = nnet,
                                    deg = deg)
  .render_reactive_network(reactive_network = reactive_network,
                         input = input,
                         output = output)
}

#' Server genomic venn
.server_venn <- function(input,
                         output,
                         session,
                         data_table,
                         deg=FALSE){
  reactive_venn <- .prepare_reactive_venn(data_table = data_table,
                                          input = input,
                                          output = output,
                                          deg = deg)
  output$plotly<- .render_venn(reactive_venn)
  output$table <- .render_venn_table(reactive_venn)
  output$download_csv <- download_csv(deg = deg,
                                     plotType = "venn",
                                     input=input,
                                     output=output,
                                     data_table=data_table)
}


#' Server volcano
.server_volcano <- function(input,
                            output,
                            session,
                            data_table,
                            deg=FALSE){
  reactive_volcano <- .prepare_reactive_volcano(data_table,
                                                input = input,
                                                output = output,
                                                deg = deg)
  output$plotly <- .render_volcano(reactive_volcano)
}


#' Server ridge
.server_ridge <- function(input,
                          output,
                          session,
                          data_table,
                          deg=FALSE){
  reactive_ridge <- .prepare_reactive_ridge(data_table,
                                            input = input,
                                            output = output,
                                            deg = deg)
  output$plotly <- .render_ridge(reactive_ridge)
  reactive_ridge_table <- .prepare_reactive_ridge(data_table,
                                                  input = input,
                                                  output = output,
                                                  deg = deg,
                                                  table=TRUE)
  output$table <- .render_ridge_table(reactive_ridge_table)
  output$download_csv <- download_csv(deg = deg,
                                      plotType = "ridge",
                                      input=input,
                                      output=output,
                                      data_table=data_table)
}


#' Server table
.server_table <- function(input,
                          output,
                          session,
                          data_table,
                          deg=FALSE){
  reactive_table <- .prepare_reactive_table(data_table,
                                            input=input,
                                            output=output)
  output$table <- .render_table(reactive_table)
  output$download_csv <- download_csv(plotType = "table",
                                      input=input,
                                      output=output,
                                      data_table=data_table)

}


#' Enrichment server
#' @importFrom shiny renderText
#' @importFrom DT dataTableOutput

.server_enrich_bg <- function(input,
                              output,
                              session,
                              extracted_data,
                              name,
                              tf=FALSE){

  FFUN <- ifelse(tf, run_tf_enrich, run_genomic_enrich)
  if(is.null(extracted_data)) return(NULL)
  bg_enr <- .run_bg(FFUN = FFUN,
                   input = input,
                   output = output,
                   args = list(model_results = NULL,
                               qvalueCutoff = 1,
                               pvalueCutoff = 0.05,
                               extracted_data = extracted_data,
                               pAdjustMethod="none"))
  check <- .check_reactive_bg_enrich(bg_enrich = bg_enr,
                                     input = input,
                                     output = output,
                                     session = session)
  output$check <- renderText({check()})
  if(!tf){
      gen_plot <- .reactive_gen_enrich(bg_enrich = bg_enr,
                                       input = input,
                                       output = output,
                                       session = session)
      output$dotplot <- renderPlotly({
        gen_plot()[["plot"]]
      })
      output$table <- DT::renderDataTable({
        ans <- gen_plot()[["table"]]
        if(!bg_enr$is_alive()){
          if(!is.null(ans)){
            ans <- mutate_if(ans, is.numeric, ~ round(., 3))
          }
        }
      })
      output$download_csv <- download_csv(plotType = "dotPlot",
                                          input=input,
                                          output=output,
                                          bg_enr=bg_enr)

  }else{
      tf_plot <- .reactive_tf_enrich(bg_enrich = bg_enr,
                                     input = input,
                                     output = output,
                                     session = session)
      ns <- NS(name)
      output$dotplot <- renderUI({
        plots <- tf_plot()
        plot_list <- lapply(seq_along(plots), function(i) {
          list(plotlyOutput(ns(paste0(names(plots)[i], "_plot"))),
               HTML(paste0(rep("<br>", 20), collapse = "")),
               tags$div(
                 style = 'overflow-x: auto;',
                 dataTableOutput(ns(paste0(names(plots)[i], "_table")))),
               downloadButton(ns(paste0(names(plots)[i], "_download_csv"))),
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
            #output
            if(!bg_enr$is_alive()){
              if(!is.null(ans)){
                ans <- mutate_if(ans, is.numeric, ~ round(., 3))
                }
            }
          })
          output[[paste0(names(plots)[i], "_download_csv")]] <- download_csv(
             plotType = "dotPlot",
             input=input,
             output=output,
             bg_enr=bg_enr,
             tf = TRUE,
             i=i)
        }
      })
  }
  onStop(function(){
    cat(sprintf("Session was closed, killing background process...  "))
    bg_enr$kill_tree()
  })
}

#' Circos server
.server_circos <- function(input,
                          output,
                          session,
                          data){
  reactive_circos <- .prepare_reactive_circos(data = data,
                                              input = input,
                                              output = output)
  output$gosling_plot <- renderGosling({
    gosling(component_id = "component_1",
            reactive_circos())
  })
  observe({
    export_pdf(component_id = "component_1")
  })%>%bindEvent(input$Download_pdf)
  observe({
    export_png(component_id = "component_1")
  })%>%bindEvent(input$Download_png)

}

#
#
# ns('IntegrationSelect')
# paste0("input.", ns("IntegrationSelect"),"=='gene_genomic_res'")
# ns("genomicTypeSelect")
# ns('ClassSelect')
# ns('ChrSelect')
# ns('SignificativityCriteria')
# paste0("input.", ns("SignificativityCriteria"),"=='pval'")
# ns("PvalRange")
# paste0("input.", ns("SignificativityCriteria"),"=='FDR'")
# ns("FdrRange")
# ns("scaleHeatmap")
# ns("numTopGenesHeatmapmirna_cnv")
# ns("numTopGenesHeatmapMETonly")
# ns("numTopGenesHeatmapCNVonly")
# ns("numTopGenesHeatmapMET")
# ns("numTopGenesHeatmapCNV")
# ns('degSelect')
# ns("circosType")
# ns("Download_png")
# ns("Download_pdf")
