# Server Genomic integration chr distribution
.server_histo <- function(input,
                          output,
                          session,
                          data_table) {
    reactive_histo <- .prepare_reactive_histo(data_table,
        input = input,
        output = output,
        deg = FALSE
    )
    output$plotly <- .render_histo(reactive_histo)
    reactive_histo_table <- .prepare_reactive_histo(data_table,
        input = input,
        output = output,
        deg = FALSE,
        table = TRUE
    )
    output$table <- .render_histo_table(reactive_histo_table)
    output$download_csv <- .download_csv(
        deg = FALSE,
        plotType = "histo",
        input = input,
        output = output,
        data_table = data_table
    )
}

# Server transcriptional and deg integration
.server_histo2 <- function(input,
                           output,
                           session,
                           data_table,
                           deg = FALSE) {
    reactive_histo_transcript <- .prepare_reactive_histo(data_table,
        input = input,
        output = output,
        deg = deg
    )
    output$plotly <- .render_histo(reactive_histo_transcript)
    reactive_histo_table_transcript <- .prepare_reactive_histo(data_table,
        input = input,
        output = output,
        deg = deg,
        table = TRUE
    )
    output$table <- .render_histo_table(reactive_histo_table_transcript)
    output$download_csv <- .download_csv(
        deg = deg,
        plotType = "histo",
        input = input,
        output = output,
        data_table = data_table
    )

    reactive_histo_tf_transcript <- .prepare_reactive_histo_tf(data_table,
        input = input,
        output = output,
        deg = deg
    )
    output$plotly_tf <- .render_histo_TF(reactive_histo_tf_transcript)
    output$download_csv_tf <- .download_csv(
        deg = deg,
        tf = TRUE,
        plotType = "histo",
        input = input,
        output = output,
        data_table = data_table
    )
    output$table_tf <- .render_histo_tf_table(reactive_histo_tf_transcript)
}

# Server network
.server_network <- function(input,
                            output,
                            session,
                            data_table,
                            deg = FALSE,
                            degs = NULL) {
    nnet <- .prepare_network(data_table)
    if (is.null(nnet)) {
        return(NULL)
    }
    reactive_network <- .select_network(
        data_table = data_table,
        input = input,
        output = output,
        network_data = nnet,
        deg = deg,
        degs = degs
    )
    .render_reactive_network(
        reactive_network = reactive_network,
        input = input,
        output = output
    )
}

# Server genomic venn
.server_venn <- function(input,
                         output,
                         session,
                         data_table,
                         deg = FALSE) {
    reactive_venn <- .prepare_reactive_venn(
        data_table = data_table,
        input = input,
        output = output,
        deg = deg
    )
    output$plotly <- .render_venn(reactive_venn)
    output$table <- .render_venn_table(reactive_venn)
    output$download_csv <- .download_csv(
        deg = deg,
        plotType = "venn",
        input = input,
        output = output,
        data_table = data_table
    )
}


# Server volcano
.server_volcano <- function(input,
                            output,
                            session,
                            data_table,
                            deg = FALSE) {
    reactive_volcano <- .prepare_reactive_volcano(data_table,
        input = input,
        output = output,
        deg = deg
    )
    output$plotly <- .render_volcano(reactive_volcano)
}


# Server ridge
.server_ridge <- function(input,
                          output,
                          session,
                          data_table,
                          deg = FALSE) {
    reactive_ridge <- .prepare_reactive_ridge(data_table,
        input = input,
        output = output,
        deg = deg
    )
    output$plotly <- .render_ridge(reactive_ridge)
    reactive_ridge_table <- .prepare_reactive_ridge(data_table,
        input = input,
        output = output,
        deg = deg,
        table = TRUE
    )
    output$table <- .render_ridge_table(reactive_ridge_table)
    output$download_csv <- .download_csv(
        deg = deg,
        plotType = "ridge",
        input = input,
        output = output,
        data_table = data_table
    )
}


# Server table
.server_table <- function(input,
                          output,
                          session,
                          data_table,
                          deg = FALSE) {
    reactive_table <- .prepare_reactive_table(data_table,
        input = input,
        output = output
    )
    output$table <- .render_table(reactive_table)
    output$download_csv <- .download_csv(
        plotType = "table",
        input = input,
        output = output,
        data_table = data_table
    )
}


# Enrichment server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly renderPlotly plotlyOutput

.server_enrich_bg <- function(input,
                              output,
                              session,
                              extracted_data,
                              name,
                              tf = FALSE) {
    FFUN <- ifelse(tf, run_tf_enrich, run_genomic_enrich)
    if (is.null(extracted_data)) {
        return(NULL)
    }
    bg_enr <- .run_bg(
        FFUN = FFUN,
        input = input,
        output = output,
        args = list(
            model_results = NULL,
            qvalueCutoff = 1,
            pvalueCutoff = 0.05,
            extracted_data = extracted_data,
            pAdjustMethod = "none"
        )
    )
    check <- .check_reactive_bg_enrich(
        bg_enrich = bg_enr,
        input = input,
        output = output,
        session = session
    )
    output$check <- renderText({
        check()
    })
    if (!tf) {
        gen_plot <- .reactive_gen_enrich(
            bg_enrich = bg_enr,
            input = input,
            output = output,
            session = session
        )
        output$dotplot <- renderPlotly({
            gen_plot()[["plot"]]
        })
        output$table <- renderDataTable({
            ans <- gen_plot()[["table"]]
            if (!bg_enr$is_alive()) {
                if (!is.null(ans)) {
                    ans <- mutate_if(ans, is.numeric, ~ round(., 3))
                }
            }
        })
        output$download_csv <- .download_csv(
            plotType = "dotPlot",
            input = input,
            output = output,
            bg_enr = bg_enr
        )
    } else {
        tf_plot <- .reactive_tf_enrich(
            bg_enrich = bg_enr,
            input = input,
            output = output,
            session = session
        )
        ns <- NS(name)
        output$dotplot <- renderUI({
            plots <- tf_plot()
            plot_list <- lapply(seq_along(plots), function(i) {
                list(
                    plotlyOutput(ns(paste0("plot_", i))),
                    HTML(paste0(rep("<br>", 20), collapse = "")),
                    tags$div(
                        style = "overflow-x: auto;",
                        dataTableOutput(ns(paste0("table_", i)))
                    ),
                    downloadButton(ns(paste0("download_csv_", i))),
                    HTML(paste0(rep("<br>", 20), collapse = ""))
                )
            })
            plot_list <- as.list(unlist(plot_list, recursive = FALSE))
            do.call(tagList, plot_list)
        })

        observe({
            plots <- tf_plot()
            outplot_list <- lapply(seq_along(plots), function(i) {
                renderPlotly({
                    plots[[i]][["plot"]]
                })
            })
            outtable_list <- lapply(seq_along(plots), function(i) {
                renderDataTable({
                    ans <- plots[[i]][["table"]]
                    if (!is.null(ans)) {
                        ans <- mutate_if(ans, is.numeric, ~ round(., 3))
                    }
                    ans
                })
            })
            outdown_list <- lapply(seq_along(plots), function(i) {
                .download_csv(
                    plotType = "dotPlot",
                    input = input,
                    output = output,
                    bg_enr = bg_enr,
                    tf = TRUE,
                    i = i
                )
            })
            if (!bg_enr$is_alive()) {
                for (i in seq_along(outplot_list)) {
                    output[[paste0("plot_", i)]] <- outplot_list[[i]]
                    output[[paste0("table_", i)]] <- outtable_list[[i]]
                    output[[paste0("download_csv_", i)]] <- outdown_list[[i]]
                }
            }
        })
    }
    onStop(function() {
        message("Session was closed, killing background process...")
        bg_enr$kill_tree()
    })
}


# Differential Methylation server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom MethylMix MethylMix MethylMix_PlotModel
#' @importFrom shinyjs enable disable

.server_DiffMet_bg <- function(input,
                              output,
                              session,
                              multiomics) {
  
  tmp <- NULL
  if(!is.null(multiomics$gene_genomic_res)) tmp <- attr(multiomics$gene_genomic_res, "Class")
  if(!is.null(multiomics$gene_met_res)) tmp <- attr(multiomics$gene_met_res, "Class")
  if(is.null(tmp)) return(NULL)
  temp_dir <- tempdir() 
  logfile <- file.path(temp_dir, "logfile.txt")
  datafile <- file.path(temp_dir, "datafile.rds")
  pid_file <- file.path(temp_dir, "pid_file.txt")
  mmultiomics_file <- file.path(temp_dir, "mmultiomics_file.rds")
  status <- reactiveValues(text = 'In attesa...')
  child_pids <- reactiveValues(main_pid = NULL, child_pids = NULL)

  logData <- reactivePoll(1000, session,
                          checkFunc = function() {
                            if (file.exists(logfile)) file.info(logfile)$mtime[1] else NULL
                          },
                          valueFunc = function() {
                            if (file.exists(logfile)) readLines(logfile) else ""
                          }
  )
  if("gene_genomic_res"%in%names(multiomics)){
    .diffMet_start(input=input,
                   status=status,
                   session=session,
                   multiomics=multiomics,
                   mmultiomics_file=mmultiomics_file,
                   pid_file=pid_file,
                   datafile=datafile,
                   logfile=logfile,
                   child_pids=child_pids)
    
    .diffMet_stop(input=input,
                   status=status,
                   session=session,
                   mmultiomics_file=mmultiomics_file,
                   pid_file=pid_file,
                   datafile=datafile,
                   logfile=logfile,
                   child_pids=child_pids)
  }else{
    if("gene_met_res"%in%names(multiomics)){
      result <- list()
      result$res_expr_cnv <- t(multiomics$gene_met_res$data$response_var)
      attr(result, "Class") <- attr(multiomics$gene_met_res, "Class")
      saveRDS(result, datafile)
      fileConn <- file(logfile)
      writeLines("COMPLETED", fileConn)
      close(fileConn)
    }else{
    return(NULL)
  }}
  observeEvent(status, {
    output$status <- renderPrint({ status$text })
    })
  
  observe({
    if (grepl("COMPLETED", logData())) {
      status$text <- "Calcolo completato!"
      if (file.exists(datafile)) {
        result_data <- readRDS(datafile)
        cclass <- attr(result_data, "Class")
        unlink(logfile)
        unlink(datafile)
        unlink(mmultiomics_file)
        unlink(pid_file) 
        output$results <- renderPrint({ names(result_data) })
        message(paste("Computing"))
        updateSelectInput(session, "class1", choices = names(table(cclass)), selected = names(table(cclass))[1])
        updateSelectInput(session, "class2", choices = names(table(cclass)), selected = names(table(cclass))[2])
        data <- reactiveVal(NULL)
        .reactive_methylmix(input=input,
                            result_data=result_data,
                            multiomics=multiomics,
                            session=session,
                            data=data,
                            cclass=cclass)
        .reactive_methylmix_plots(input=input,
                                  session=session,
                                  data=data,
                                  output=output)
      }
    }
  })
  
  }
  

# Differential CNV server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom MethylMix MethylMix MethylMix_PlotModel
#' @importFrom shinyjs enable disable

.server_DiffCNV <- function(input,
                           output,
                           session,
                           multiomics) {
  
  tmp <- NULL
  if(!is.null(multiomics$gene_genomic_res)) tmp <- attr(multiomics$gene_genomic_res, "Class")
  if(!is.null(multiomics$gene_cnv_res)) tmp <- attr(multiomics$gene_met_res, "Class")
  if(is.null(tmp)) return(NULL)
  if("gene_genomic_res"%in%names(multiomics)){
    cclass <- attr(multiomics$gene_genomic_res, "Class")
    updateSelectInput(session, "class1", choices = names(table(cclass)), selected = names(table(cclass))[1])
    updateSelectInput(session, "class2", choices = names(table(cclass)), selected = names(table(cclass))[2])
    observeEvent(list(input$class1, input$class2), {
      req(cclass, input$class1, input$class2)
      disable("genes")
      cnv <- multiomics$gene_genomic_res$data$covariates
      cnv <- cnv[, grep("_cnv", colnames(cnv))]
      colnames(cnv) <- gsub("_cnv", "", colnames(cnv))
      print(cnv[1:5, 1:5])
      message(input$class1)
      ans <- lapply(colnames(cnv), function(x){
        t.test(cnv[names(cclass)[cclass==input$class1], x],
               cnv[names(cclass)[cclass==input$class2], x])
      })
      names(ans) <- colnames(cnv)
      ans <- lapply(ans, function(x){
        as.data.frame(t(c(x$statistic, pvalue=x$p.value, conf=x$conf.int, x$estimate)))
      })
      tmp <- plyr::rbind.fill(ans)
      rownames(tmp) <- names(ans)
      ans <- tmp
      ans$FDR <- p.adjust(ans$pvalue, method = "fdr")
      ans <- ans[ans$pvalue<=0.05,]
      output$table <- renderDataTable({round(ans, digits = 5)})
      
    })
    
      
    }else{
    if("gene_met_res"%in%names(multiomics)){
     
    }else{
      return(NULL)
    }}
  # observeEvent(status, {
  #   output$status <- renderPrint({ status$text })
  # })
  # 
  # observe({
  # 
  # })
  
}




# Circos server
#' @importFrom shiny.gosling renderGosling gosling export_png export_pdf
.server_circos <- function(input,
                           output,
                           session,
                           data) {
    reactive_circos <- .prepare_reactive_circos(
        data = data,
        input = input,
        output = output
    )
    output$gosling_plot <- renderGosling({
        gosling(
            component_id = "component_1",
            reactive_circos()
        )
    })
    observe({
        export_pdf(component_id = "component_1")
    }) %>% bindEvent(input$Download_pdf)
    observe({
        export_png(component_id = "component_1")
    }) %>% bindEvent(input$Download_png)
}
