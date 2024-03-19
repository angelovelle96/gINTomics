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
                            deg = FALSE) {
    nnet <- .prepare_network(data_table)
    if (is.null(nnet)) {
        return(NULL)
    }
    reactive_network <- .select_network(
        data_table = data_table,
        input = input,
        output = output,
        network_data = nnet,
        deg = deg
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
