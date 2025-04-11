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
    if(deg & length(grep("^deg_", colnames(data_table)))==0) return(NULL)
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
    if(deg & is.null(degs)) return(NULL)
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
    if(deg & length(grep("^deg_", colnames(data_table)))==0) return(NULL)
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
   if(deg & length(grep("^deg_", colnames(data_table)))==0) return(NULL)
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
    if(deg & length(grep("^deg_", colnames(data_table)))==0) return(NULL)
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

# Server table DEGs
.server_degsTable <- function(input,
                          output,
                          session,
                          multiomics) {
  if(is.null(.extract_deg_info(multiomics))) return(NULL)
  reactive_degsTable <- .prepare_reactive_degsTable(multiomics = multiomics,
                                                    input = input,
                                                    output = output)
  output$table <- .render_table(reactive_degsTable)
  output$download_csv <- .download_csv(
    plotType = "tableDegs",
    input = input,
    output = output,
    data_table = multiomics
  )
}


# Enrichment genomic server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly renderPlotly plotlyOutput

.server_enrich_bg <- function(input,
                              output,
                              session,
                              extracted_data) {

    bg_enr <- .run_bg(
        FFUN = run_genomic_enrich,
        args = list(
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
      onStop(function() {
        message("Session was closed, killing background process...")
        if (!is.null(bg_process)) {
          if (bg_enr$is_alive()) {
            bg_enr$kill_tree()
            message("Background process killed.")
          } else {
            message("No background process running.")
          }
        }
      })
    
}



# Enrichment TF server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' renderPlot updateSelectizeInput observe
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly renderPlotly plotlyOutput

.server_enrich_bgTF <- function(input,
                              output,
                              session,
                              extracted_data) {
  
      bg_process <<- NULL
      tf_data <- .run_reactive_tf_enrich(
      extracted_data = extracted_data,
      input = input,
      output = output,
      session = session
    )
    tf_plot <- .reactive_tf_enrich(
      bg_enrich = tf_data,
      input = input,
      output = output,
      session = session
    )
    observe({
      data <- extracted_data
      data <- data[data$cov != "(Intercept)", ]
      tmp <- unique(data$cov)
      tmp2 <- lapply(tmp, function(y) {
        check <- data$pval[data$cov == y]
        check <- sum(check <= 0.01)
        return(check)
      })
      names(tmp2) <- tmp
      tmp2 <- unlist(tmp2)
      tmp2 <- tmp2[tmp2 > 12]
      tmp2 <- sort(tmp2, decreasing = TRUE)
      
      updateSelectizeInput(session, "genes",
                           choices = names(tmp2),
                           selected = names(tmp2)[1],
                           server = TRUE)
    })
    observe({
      tf_plot()
    })
    check <- .check_reactive_bg_enrich(
      bg_enrich = tf_data,
      input = input,
      output = output,
      session = session
    )
    output$check <- renderText({
      check()
    })
    
    output$plot <- renderPlotly({
      plots <- tf_plot()
      plots[["plot"]]
    })
    
    output$table <- renderDataTable({
      plots <- tf_plot()
      ans <- plots[["table"]]
      if (!is.null(ans)) {
        ans <- mutate_if(ans, is.numeric, ~ round(., 3))
      }
      ans
    })
    
    onStop(function() {
      message("Session was closed, killing background process...")
      if (!is.null(bg_process)) {
        if (bg_process$is_alive()) {
          bg_process$kill_tree()
          message("Background process killed.")
        } else {
          message("No background process running.")
        }
      rm(bg_process, envir = .GlobalEnv)
      }
    })
  }


# Differential Methylation server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' updateSelectInput observeEvent observe reactiveVal reactiveValues 
#' reactivePoll renderPrint req
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom MethylMix MethylMix MethylMix_PlotModel
#' @importFrom shinyjs enable disable

.server_DiffMet_bg <- function(input,
                              output,
                              session,
                              multiomics,
                              data_table) {
  
  tmp <- NULL
  if(!is.null(multiomics$gene_genomic_res)){
    tmp <- attr(multiomics$gene_genomic_res, "Class")
    }
  if(!is.null(multiomics$gene_met_res)){
    tmp <- attr(multiomics$gene_met_res, "Class")
    }
  if(is.null(tmp)) return(NULL)
  temp_dir <- tempdir() 
  logfile <- file.path(temp_dir, "logfile.txt")
  datafile <- file.path(temp_dir, "datafile.rds")
  pid_file <- file.path(temp_dir, "pid_file.txt")
  mmultiomics_file <- file.path(temp_dir, "mmultiomics_file.rds")
  status <- reactiveValues(text = 'Waiting for start input...')
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
  
  observeEvent(list(input$class1, input$class2), {
    req(input$class1, input$class2)
    if(input$class1==input$class2){
      status$text <- paste("Please select a valid contrast")
    }else{
      status$text <- paste("Running MethylMix to detect",
                           "differentially methylated genes")
      }
    })
  observeEvent(status, {
    output$status <- renderPrint({ status$text })
    })
  
  observe({
    req(input$start)
    if (grepl("COMPLETED", logData())) {
      status$text <- paste("Running MethylMix to detect",
                           "differentially methylated genes")
      if (file.exists(datafile)) {
        result_data <- readRDS(datafile)
        cclass <- attr(result_data, "Class")
        unlink(logfile)
        unlink(datafile)
        unlink(mmultiomics_file)
        unlink(pid_file) 
        updateSelectInput(session, "class1",
                          choices = names(table(cclass)),
                          selected = names(table(cclass))[1])
        updateSelectInput(session, "class2",
                          choices = names(table(cclass)),
                          selected = names(table(cclass))[2])
        ddata <- reactiveVal(NULL)
        .reactive_methylmix(input=input,
                            result_data=result_data,
                            multiomics=multiomics,
                            session=session,
                            ddata=ddata,
                            cclass=cclass,
                            status = status)
        .reactive_methylmix_plots(input=input,
                                  session=session,
                                  ddata=ddata,
                                  output=output,
                                  data_table = data_table)
      }
    }
  })
  
  }
  

# Differential CNV server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' updateSelectInput
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom MethylMix MethylMix MethylMix_PlotModel
#' @importFrom shinyjs enable disable

.server_DiffCNV <- function(input,
                           output,
                           session,
                           multiomics) {
  
  tmp <- NULL
  if(!is.null(multiomics$gene_genomic_res)){
    tmp <- attr(multiomics$gene_genomic_res, "Class")
    }
  if(!is.null(multiomics$gene_cnv_res)){
    tmp <- attr(multiomics$gene_cnv_res, "Class")
    }
  if(is.null(tmp)) return(NULL)
  status <- reactiveValues(text = 'Waiting for start input...')
  observeEvent(input$start, {
      status$text <- paste("Running...")
  })
  observeEvent(status, {
    output$status <- renderPrint({ status$text })
  })
  if("gene_genomic_res"%in%names(multiomics)){
    cclass <- attr(multiomics$gene_genomic_res, "Class")
    cnv <- multiomics$gene_genomic_res$data$covariates
    cnv <- cnv[, grep("_cnv", colnames(cnv))]
    colnames(cnv) <- gsub("_cnv", "", colnames(cnv))
    eexpr <- multiomics$gene_genomic_res$data$response_var
    }else{
    if("gene_cnv_res"%in%names(multiomics)){
      cclass <- attr(multiomics$gene_cnv_res, "Class")
      cnv <- multiomics$gene_cnv_res$data$covariates
      colnames(cnv) <- gsub("_cov", "", colnames(cnv))
      eexpr <- multiomics$gene_cnv_res$data$response_var
    }else{
      return(NULL)
    }}
  updateSelectInput(session, "class1", choices = names(table(cclass)),
                    selected = names(table(cclass))[1])
  updateSelectInput(session, "class2", choices = names(table(cclass)),
                    selected = names(table(cclass))[2])
  ddata <- .reactive_cnv_test(input=input,
                     session=session,
                     cnv=cnv,
                     output=output,
                     cclass=cclass)
  .render_cnv_test_table(input=input,
                         session=session,
                         ddata=ddata,
                         output=output,
                         status=status)
  pplots <- .reactive_cnv_test_plots(input=input,
                           session=session,
                           cnv=cnv,
                           output=output,
                           eexpr=eexpr,
                           cclass=cclass)
  output$plot1 <- renderPlot({pplots()[[1]]})
  output$plot2 <- renderPlot({pplots()[[2]]})
}


# ExprIns genomic server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' updateSelectizeInput renderPlot
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom MethylMix MethylMix MethylMix_PlotModel
#' @importFrom shinyjs enable disable
#' @importFrom ggplot2 labs

.server_exprInsGen <- function(input,
                              output,
                              session,
                              multiomics) {
  if("gene_genomic_res"%in%names(multiomics)){
    ggenes <- colnames(multiomics$gene_genomic_res$data$response_var)
    cclass <- attr(multiomics$gene_genomic_res, "Class")
    cnv <- multiomics$gene_genomic_res$data$covariates
    met <- cnv[, grep("_met", colnames(cnv))]
    cnv <- cnv[, grep("_cnv", colnames(cnv))]
    colnames(cnv) <- gsub("_cnv", "", colnames(cnv))
    colnames(met) <- gsub("_met", "", colnames(met))
    eexpr <- multiomics$gene_genomic_res$data$response_var
  }else{
    if("gene_cnv_res"%in%names(multiomics)){
      ggenes <- colnames(multiomics$gene_cnv_res$data$response_var)
      cclass <- attr(multiomics$gene_cnv_res, "Class")
      cnv <- multiomics$gene_cnv_res$data$covariates
      colnames(cnv) <- gsub("_cov", "", colnames(cnv))
      eexpr <- multiomics$gene_cnv_res$data$response_var
    }else{
      if("gene_met_res"%in%names(multiomics)){
        ggenes <- colnames(multiomics$gene_met_res$data$response_var)
        cclass <- attr(multiomics$gene_met_res, "Class")
        met <- multiomics$gene_met_res$data$covariates
        colnames(met) <- gsub("_cov", "", colnames(met))
        eexpr <- multiomics$gene_met_res$data$response_var
      }else{
      ggenes <- NULL  
      cnv <- NULL
      eexpr <- NULL
      }}}
  if("mirna_cnv_res"%in%names(multiomics)){
    ggenes <- c(ggenes, colnames(multiomics$mirna_cnv_res$data$response_var))
    cclass <- attr(multiomics$mirna_cnv_res, "Class")
    if(!is.null(cnv)){
      cnv <- cbind(cnv,
                   multiomics$mirna_cnv_res$data$covariates[rownames(cnv),])
      eexpr <- cbind(eexpr,
                  multiomics$mirna_cnv_res$data$response_var[rownames(eexpr),])
      colnames(cnv) <- gsub("_cov", "", colnames(cnv))
    }else{
      cnv <- multiomics$mirna_cnv_res$data$covariates
      eexpr <- multiomics$mirna_cnv_res$data$response_var
      colnames(cnv) <- gsub("_cov", "", colnames(cnv))
      }
  }
  updateSelectizeInput(session, "genes",
                       choices = ggenes,
                       selected = ggenes[1],
                       server = TRUE)
  if(is.null(cclass)) cclass <- setNames(rep("", nrow(eexpr)), rownames(eexpr))
  pplots <- .reactive_cnv_test_plots(input=input,
                                     session=session,
                                     cnv=cnv,
                                     output=output,
                                     eexpr=eexpr,
                                     cclass=cclass)
  pplots2 <- .reactive_cnv_test_plots(input=input,
                                     session=session,
                                     cnv=met,
                                     output=output,
                                     eexpr=eexpr,
                                     cclass=cclass)
  output$plot1 <- renderPlot({pplots()[[1]]})
  output$plot2 <- renderPlot({pplots()[[2]]})
  output$plot3 <- renderPlot({
    pplots2()[[1]]+labs(
      title = paste(input$genes, "scatterplot"),
      subtitle = "Correlation between gene expression and Methylation",
      x = "Methylation"
    )
    })
  output$plot4 <- renderPlot({pplots2()[[2]]+labs(
    title = paste(input$genes, "Methylation values"),
    y = "Methylation"
  )
    })
}

# ExprIns genomic server
#' @importFrom shiny renderText downloadButton HTML tags renderUI tagList onStop
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom MethylMix MethylMix MethylMix_PlotModel
#' @importFrom shinyjs enable disable
#' @importFrom ggplot2 labs

.server_exprInsTransc <- function(input,
                                 output,
                                 session,
                                 multiomics) {
  
  observeEvent(input$IntegrationSelect, {
    data <- multiomics[[input$IntegrationSelect]]
    cov <- data$data$covariates
    res <- data$data$response
    ggenes1 <- colnames(res)
    ggenes2 <- colnames(cov)
    updateSelectizeInput(session, "gene1",
                         choices = ggenes1,
                         selected = ggenes1[1],
                         server = TRUE)
    updateSelectizeInput(session, "gene2",
                         choices = ggenes2,
                         selected = ggenes2[1],
                         server = TRUE)
  })
  
  pplots <- .reactive_transcExpr_plots(input=input,
                                      session=session,
                                      multiomics=multiomics,
                                      output=output)
  output$plot1 <- renderPlot({pplots()[[1]]})
  output$plot2 <- renderPlot({pplots()[[2]]})
  
  
  
  
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
