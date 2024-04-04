# Prepare network data for visualization
.prepare_network <- function(data_table) {
    tf_list <- c()
    mirna_list <- c()
    target_list <- c()
    tmp <- intersect(c("cov", "response", "pval", "fdr", "coef", "class",
                       "omics"), colnames(data_table))
    all <- subset(data_table, `omics` %in% c(
        "tf_res",
        "tf_mirna_res",
        "mirna_target_res"
    ),
    select = tmp
    )
    if (nrow(all) == 0) {
        return(NULL)
    }
    nodes <- data.frame(gene = c(all$cov, all$response))
    edges <- subset(all, select = tmp)
    nodes <- unique(nodes)
    names(nodes) <- "id"
    names(edges) <- gsub("cov", "from", names(edges))
    names(edges) <- gsub("response", "to", names(edges))
    tf_list <- all[all$omics == "tf_res", "cov"]
    tf_list2 <- all[all$omics == "tf_mirna_res", "cov"]
    mirna_list <- all[all$omics == "tf_mirna_res", "response"]
    mirna_list2 <- all[all$omics == "mirna_target_res", "cov"]
    target_list <- all[all$omics == "tf_res", "response"]
    target_list2 <- all[all$omics == "tf_mirna_res", "cov"]
    tf_list <- unique(c(tf_list, tf_list2))
    mirna_list <- unique(c(mirna_list, mirna_list2))
    target_list <- unique(c(target_list, target_list2))
    nodes$label <- nodes$id
    nodes$shape <- ifelse(nodes$id %in% tf_list,
        "diamond",
        ifelse(nodes$id %in% mirna_list,
            "triangle",
            ifelse(nodes$id %in% target_list,
                "circle",
                "diamond"
            )
        )
    )
    nodes$title <- nodes$id
    nodes$shadow <- TRUE
    nodes$color <- ifelse(nodes$id %in% tf_list,
        "#FFBA01",
        ifelse(nodes$id %in% mirna_list,
            "#9969C7",
            ifelse(nodes$id %in% target_list,
                "#CCE7C9",
                "#FFBA01"
            )
        )
    )
    nodes$width <- 10
    edges$width <- abs(edges$coef) * 10
    edges$width[edges$width > 40] <- 40
    edges$color <- ifelse(edges$coef > 0,
        "#4169E1",
        ifelse(edges$coef < 0,
            "#ED5564",
            "black"
        )
    )
    edges$length <- 500
    edges$pval <- round(edges$pval, 3)
    edges$fdr <- round(edges$fdr, 3)
    edges$coef <- round(edges$coef, 3)
    edges$title <- paste(
        "pval:",
        edges$pval,
        ",coef:", edges$coef,
        ",FDR:", edges$fdr
    )

    legend_nodes <- data.frame(
        label = c(
            "TF",
            "Target",
            "miRNA"
        ),
        shape = c(
            "diamond",
            "circle",
            "triangle"
        ),
        color = c(
            "#FFBA01",
            "#CCE7C9",
            "#9969C7"
        )
    )

    legend_edges <- data.frame(
        color = c(
            "#4169E1",
            "#ED5564",
            "black"
        ),
        label = c(
            "UPregulate",
            "DOWNregulate",
            "no-effect"
        ),
        arrows = c(
            "to",
            "to",
            "to"
        )
    )
    edges <- edges[order(abs(edges$coef), decreasing = TRUE), ]
    return(list(
        nodes = nodes,
        edges = edges,
        legend_edges = legend_edges,
        legend_nodes = legend_nodes
    ))
}


# Prepare network for visualization
#' @importFrom visNetwork visNetwork visGroups visLegend visEdges visOptions
#'   visLayout visExport
.build_network <- function(nodes,
                           edges,
                           legend_nodes,
                           legend_edges) {
    visNetwork(nodes,
        edges,
        width = "100%"
    ) %>%
        visGroups(
            groupname = "TF",
            color = "#FFBA01"
        ) %>%
        visGroups(
            groupname = "Target",
            color = "#CCE7C9"
        ) %>%
        visGroups(
            groupname = "miRNA",
            color = "#9969C7"
        ) %>%
        visLegend(
            addEdges = legend_edges,
            addNodes = legend_nodes,
            useGroups = FALSE,
            width = 0.2,
            position = "right",
            main = "Legend"
        ) %>%
        visEdges(
            shadow = TRUE,
            smooth = FALSE,
            arrows = list(to = list(enabled = TRUE)),
            color = list(color = "black", highlight = "red")
        ) %>%
        visOptions(
            highlightNearest = list(
                enabled = TRUE,
                degree = 2,
                hover = TRUE
            ),
            nodesIdSelection = TRUE,
            manipulation = TRUE
        ) %>%
        visLayout(randomSeed = 20, improvedLayout = TRUE) %>%
        visExport(
            type = "pdf", name = "export-network",
            float = "left", label = "Save network"
        )
}


# Build a Venn diagram
#' @importFrom ggvenn ggvenn
#' @importFrom RColorBrewer brewer.pal
.build_venn <- function(venn_data) {
    cnv_sign_genes <- venn_data$cnv_sign_genes
    met_sign_genes <- venn_data$met_sign_genes
    list_venn <- list(
        CNV = cnv_sign_genes$cov,
        MET = met_sign_genes$cov
    )
    ccol <- brewer.pal(name = "Set3", n = max(3, length(list_venn)))
    venn_diagram <- ggvenn(list_venn, fill_color = ccol)
    return(venn_diagram)
}

# Render a Venn diagram
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom ggplot2 ggplot theme_void labs
.render_venn <- function(reactive_venn) {
    renderPlotly({
        venn_data <- reactive_venn()
        if (!is.null(venn_data)) {
            venn_diagram <- .build_venn(venn_data)
            venn_diagram <- ggplotly(venn_diagram)
            return(venn_diagram)
        } else {
            ans <- ggplot() +
                theme_void() +
                labs(title = "No data available")
            return(ggplotly(ans))
        }
    })
}

# Render Venn diagram as a table
#' @importFrom DT renderDataTable
.render_venn_table <- function(reactive_venn) {
    renderDataTable({
        venn_data <- reactive_venn()

        if (!is.null(venn_data$cnv_sign_genes) && !is.null(
          venn_data$met_sign_genes)) {
            cnv_genes <- unlist(venn_data$cnv_sign_genes)
            met_genes <- unlist(venn_data$met_sign_genes)

            if (length(cnv_genes) > 0 & length(met_genes) > 0) {
                common_genes <- intersect(cnv_genes, met_genes)
                data.frame(Genes = common_genes)
            } else if (length(cnv_genes) == 1 & length(met_genes) == 1) {
                data.frame(CNV_Genes = cnv_genes, MET_Genes = met_genes,
                           stringsAsFactors = FALSE)
            } else {
                data.frame(Genes = character(), stringsAsFactors = FALSE)
            }
        } else if (!is.null(venn_data$cnv_sign_genes)) {
            data.frame(CNV_Genes = unlist(venn_data$cnv_sign_genes),
                       stringsAsFactors = FALSE)
        } else if (!is.null(venn_data$met_sign_genes)) {
            data.frame(MET_Genes = unlist(venn_data$met_sign_genes),
                       stringsAsFactors = FALSE)
        } else {
            data.frame(Genes = character(), stringsAsFactors = FALSE)
        }
    })
}

# Build a volcano plot
#' @importFrom plotly plot_ly
.build_volcano <- function(volcano_data) {
    plot_ly(volcano_data,
        width = 900,
        height = 700,
        x = ~coef,
        y = ~pval_fdr,
        mode = "markers",
        type = "scatter",
        color = ~group,
        colors = "Set2",
        symbol = ~class,
        text = ~ paste(
            "Group:", group, "<br>",
            "Class:", class, "<br>",
            "Name:", cov, "<br>",
            "Pval/FDR(-log10):", pval_fdr, "<br>",
            "coef", coef
        ),
        textposition = "top right"
    ) %>%
        layout(title = "Volcano Plot")
}

# Render a volcano plot
#' @importFrom plotly renderPlotly
#' @importFrom ggplot2 ggplot theme_void labs
.render_volcano <- function(reactive_volcano, annotations) {
    renderPlotly({
        volcano_data <- reactive_volcano()
        if (!is.null(volcano_data)) {
            volcano_plot <- .build_volcano(volcano_data)
            return(volcano_plot)
        } else {
            ans <- ggplot() +
                theme_void() +
                labs(title = "No data available")
            return(ggplotly(ans))
        }
    })
}

# Build a ridgeline plot
#' @importFrom ggplot2 ggplot labs theme_minimal scale_x_continuous theme
#' @importFrom ggridges geom_density_ridges theme_ridges position_raincloud
.build_ridge <- function(ridge_data,
                         quantiles) {
    ggplot(
        ridge_data,
        aes(
            x = `coef`,
            y = `significance`,
            fill = significance
        )
    ) +
        geom_density_ridges(
            jittered_points = TRUE,
            quantile_lines = TRUE,
            vline_width = 1,
            vline_color = "red",
            point_size = 0.4,
            point_alpha = 1,
            position = position_raincloud(adjust_vlines = TRUE)
        ) +
        labs(
            title = "Ridgeline Plot",
            x = "Value",
            y = "Significativity"
        ) +
        theme_minimal() +
        theme_ridges() +
        scale_x_continuous(limits = quantiles) +
        theme(legend.position = "bottom")
}

# Render a ridgeline plot
#' @importFrom shiny renderPlot
.render_ridge <- function(reactive_ridge) {
    renderPlot({
        tmp <- reactive_ridge()
        if (!is.null(tmp)) {
            ridge_data <- tmp$df
            quantiles <- tmp$quantiles
            ridge_plot <- .build_ridge(
                ridge_data = ridge_data,
                quantiles = quantiles
            )
            return(ridge_plot)
        } else {
            plot(1,
                type = "n", xlab = "", ylab = "", xlim = c(0, 1),
                ylim = c(0, 1),
                main = "No data available"
            )
        }
    })
}

# Build a histogram
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal
.build_histo <- function(histo_data) {
    ggplot(
        histo_data,
        aes(x = factor(`chr_cov`), fill = `significance`)
    ) +
        geom_bar() +
        labs(
            title = "Number of Genes with Significant Coefficients by
            Chromosome",
            x = "Chromosome",
            y = "Count"
        ) +
        theme_minimal()
}

# Render a histogram
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom ggplot2 ggplot theme_void labs
.render_histo <- function(reactive_histo) {
    renderPlotly({
        histo_data <- reactive_histo()
        if (!is.null(histo_data)) {
            histo_plot <- ggplotly(.build_histo(histo_data = histo_data))
            return(histo_plot)
        } else {
            ans <- ggplot() +
                theme_void() +
                labs(title = "No data available")
            return(ggplotly(ans))
        }
    })
}

# Render histogram as a table
#' @importFrom DT renderDataTable
.render_histo_table <- function(reactive_histo_table) {
    renderDataTable({
        table_data <- reactive_histo_table()
        if (!is.null(table_data)) {
            ttable <- .build_table(table_data)
            return(ttable)
        } else {
            return(NULL)
        }
    })
}

# Build histogram for TFs/miRNAs by chromosome
#' @importFrom plotly plot_ly layout
#' @importFrom stats reorder
#' @importFrom dplyr %>%
.build_histo_TFbyChr <- function(histo_data) {
    if (nrow(histo_data) == 0) {
        return(NULL)
    }
    plot_ly(histo_data,
        x = ~ reorder(TF, -Count), y = ~Count, type = "bar", color =
          ~Chromosome
    ) %>%
        layout(
            title = "Number of genes targeted by TFs/miRNAs",
            xaxis = list(title = "TF/miRNAs"),
            yaxis = list(title = "Number of targets"),
            barmode = "group"
        )
}

# Render histogram for TFs/miRNAs
#' @importFrom plotly renderPlotly
.render_histo_TF <- function(reactive_histo) {
    renderPlotly({
        histo_data <- reactive_histo()
        if (!is.null(histo_data)) {
            histo_plot <- .build_histo_TFbyChr(histo_data = histo_data)
            return(histo_plot)
        } else {
            ans <- ggplot() +
                theme_void() +
                labs(title = "No data available")
            return(ggplotly(ans))
        }
    })
}

# Render histogram for TFs/miRNAs as a table
#' @importFrom DT renderDataTable
.render_histo_tf_table <- function(reactive_histo_tf_table) {
    renderDataTable({
        table_data <- reactive_histo_tf_table()
        if (!is.null(table_data)) {
            ttable <- .build_table(table_data)
            return(ttable)
        } else {
            return(NULL)
        }
    })
}

# Render a ridgeline plot as a table
#' @importFrom DT renderDataTable
.render_ridge_table <- function(reactive_ridge_table) {
    renderDataTable({
        table_data <- reactive_ridge_table()
        if (!is.null(table_data)) {
            ttable <- .build_table(table_data)
            return(ttable)
        } else {
            return(NULL)
        }
    })
}

# Build a table
#' @importFrom DT datatable
.build_table <- function(table_data) {
    datatable(table_data,
        options = list(orderClasses = TRUE)
    )
}

# Render a table
#' @importFrom DT renderDataTable
.render_table <- function(reactive_table) {
    renderDataTable({
        table_data <- reactive_table()
        if (!is.null(table_data)) {
            ttable <- .build_table(table_data)
            return(ttable)
        } else {
            return(NULL)
        }
    })
}

# Handle CSV file download
#' @importFrom shiny downloadHandler
#' @importFrom utils write.table write.csv
.download_csv <- function(tf = FALSE,
                          deg = FALSE,
                          bg_enr = NULL,
                          plotType = "histo",
                          input,
                          output,
                          data_table,
                          i = NULL,
                          session = NULL) {
    if (plotType == "histo") {
        handler <- downloadHandler(
            filename = function() {
                "data_table_histo.csv"
            },
            content = function(file) {
                if (tf) {
                    reactive_histo_table <- .prepare_reactive_histo_tf(
                        deg = deg,
                        data_table = data_table,
                        input = input,
                        output = output
                    )
                } else {
                    reactive_histo_table <- .prepare_reactive_histo(
                        deg = deg,
                        data_table = data_table,
                        input = input,
                        output = output,
                        table = TRUE
                    )
                }
                write.csv(reactive_histo_table(), file, row.names = FALSE)
            }
        )
    }
    if (plotType == "ridge") {
        handler <- downloadHandler(
            filename = function() {
                "data_ridge.csv"
            },
            content = function(file) {
                reactive_ridge_table <- .prepare_reactive_ridge(
                    deg = deg,
                    data_table = data_table,
                    input = input,
                    output = output,
                    table = TRUE
                )
                write.csv(reactive_ridge_table(), file, row.names = FALSE)
            }
        )
    }
    if (plotType == "table") {
        handler <- downloadHandler(
            filename = function() {
                "data_table.csv"
            },
            content = function(file) {
                reactive_table <- .prepare_reactive_table(
                    data_table = data_table,
                    input = input,
                    output = output
                )
                write.csv(reactive_table(), file, row.names = FALSE)
            }
        )
    }
    if (plotType == "dotPlot") {
        handler <- downloadHandler(
            filename = function() {
                "enrich_table.csv"
            },
            content = function(file) {
                if (!tf) {
                    reactive_enr <- .reactive_gen_enrich(
                        bg_enrich = bg_enr,
                        input = input,
                        output = output,
                        session = session
                    )
                    write.csv(reactive_enr()[["table"]], file,
                              row.names = FALSE)
                } else {
                    reactive_enr <- .reactive_tf_enrich(
                        bg_enrich = bg_enr,
                        input = input,
                        output = output,
                        session = session
                    )
                    write.csv(reactive_enr()[[i]][["table"]], file,
                              row.names = FALSE)
                }
            }
        )
    }
    if (plotType == "venn") {
        handler <- downloadHandler(
            filename = function() {
                "venn_table.csv"
            },
            content = function(file) {
                reactive_venn <- .prepare_reactive_venn(
                    deg = deg,
                    data_table = data_table,
                    input = input,
                    output = output
                )
                venn_data <- reactive_venn()
                cnv_genes <- unlist(venn_data$cnv_sign_genes)
                met_genes <- unlist(venn_data$met_sign_genes)
                common_genes <- intersect(cnv_genes, met_genes)
                venn_data <- data.frame(Genes = common_genes)
                write.csv(venn_data, file, row.names = FALSE)
            }
        )
    }
    return(handler)
}


# Preprocess data for Circos plot
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom RColorBrewer brewer.pal
.circos_preprocess <- function(data) {
    dataframes <- lapply(unique(data$omics), function(x) {
        single_omic_df <- data[data$omics == x, ]
        if (max(single_omic_df$response_value, na.rm = TRUE) > 30) {
            single_omic_df$response_value <- log2(
              single_omic_df$response_value + 1)
        }
        return(single_omic_df)
    })
    names(dataframes) <- paste0("df_", unique(data$omics))
    if ("df_gene_genomic_res" %in% names(dataframes)) {
        dataframes$df_cnv <- filter(
            dataframes$df_gene_genomic_res,
            `cnv_met` == "cnv"
        )
        dataframes$df_met <- filter(
            dataframes$df_gene_genomic_res,
            `cnv_met` == "met"
        )
        tmp <- lapply(which(names(dataframes) != "df_gene_genomic_res"),
                      function(x) {
            ans <- dataframes[[x]]
            ans <- ans[ans$cov != "(Intercept)", ]
            ans$cov_value_original <- ans$cov_value
            ans$cov_value <- abs(ans$cov_value)
            ans$coef_original <- ans$coef
            ans$coef <- abs(ans$coef)
            return(ans)
        })
        names(tmp) <- names(dataframes)[
            which(names(dataframes) != "df_gene_genomic_res")
        ]
        dataframes <- tmp
    } else {
        tmp <- lapply(which(names(dataframes) != "df_gene_genomic_res"),
                      function(x) {
            ans <- dataframes[[x]]
            ans <- ans[ans$cov != "(Intercept)", ]
            ans$cov_value_original <- ans$cov_value
            ans$cov_value <- abs(ans$cov_value)
            ans$coef_original <- ans$coef
            ans$coef <- abs(ans$coef)
            return(ans)
        })
        names(tmp) <- names(dataframes)[
            which(names(dataframes) != "df_gene_genomic_res")
        ]
        dataframes <- tmp
    }
    gr <- lapply(dataframes, function(x) {
        ans <- makeGRangesFromDataFrame(x,
            seqnames.field = "chr",
            start.field = "start",
            end.field = "end",
            keep.extra.columns = TRUE,
            na.rm = TRUE
        )
        return(ans)
    })
    return(gr)
}

# Generate color palette for Circos track
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRamp
#' @importFrom grDevices colorRampPalette
.circos_track_cols <- function(vvalues,
                               n = 50,
                               colors = NULL) {
    if (is.null(colors)) {
        colors <- rev(brewer.pal(name = "RdBu", n = 11))
    }
    vvalues <- as.numeric(vvalues)
    if (min(vvalues) < 0) {
        breaks <- c(
            seq(min(vvalues, na.rm = TRUE), 0, len = n / 2),
            seq(0.00001, max(vvalues, na.rm = TRUE), len = n / 2)
        )
    } else {
        breaks <- c(
            seq(min(vvalues, na.rm = TRUE), mean(vvalues, na.rm = TRUE),
                len = n / 2),
            seq(mean(vvalues, na.rm = TRUE) + 0.00001,
                max(vvalues, na.rm = TRUE),
                len = n / 2
            )
        )
    }
    ii <- cut(unique(vvalues),
        breaks = breaks,
        include.lowest = TRUE
    )
    ccol <- colorRampPalette(colors)(n - 1)[ii]
    return(ccol)
}


# Create single track for Circos visualization
#' @importFrom shiny.gosling add_single_track visual_channel_x visual_channel_y
#' visual_channel_color visual_channel_tooltips visual_channel_tooltip
#' track_data_gr
.create_single_track <- function(data,
                                 dataValue = NULL,
                                 x_axis = NULL,
                                 xe_axis = NULL,
                                 y_axis = NULL,
                                 colorField = NULL,
                                 colorDomain = NULL,
                                 colorRange = NULL,
                                 tooltipField1 = NULL,
                                 tooltipTitle = NULL,
                                 tooltipAlt1 = NULL,
                                 tooltipField2 = NULL,
                                 tooltipAlt2 = NULL,
                                 tooltipField3 = NULL,
                                 tooltipAlt3 = NULL,
                                 tooltipField4 = NULL,
                                 tooltipAlt4 = NULL,
                                 tooltipField5 = NULL,
                                 tooltipAlt5 = NULL,
                                 legend = NULL,
                                 colorType = NULL,
                                 title = NULL) {
    return(
        add_single_track(
            data = track_data_gr(data, chromosomeField = "seqnames",
                                 genomicFields = c("start", "end"),
                                 value = dataValue),
            mark = "bar",
            x = visual_channel_x(field = "start", type = "genomic",
                                 axis = x_axis),
            xe = visual_channel_x(field = "end", type = "genomic",
                                  axis = xe_axis),
            y = visual_channel_y(field = dataValue, type = "quantitative",
                                 axis = y_axis),
            color = visual_channel_color(field = colorField, type = colorType,
                                         domain = colorDomain,
                                         range = colorRange, legend = legend),
            tooltip = visual_channel_tooltips(
                visual_channel_tooltip(field = "start", type = "genomic",
                                       alt = "Start Position:"),
                visual_channel_tooltip(field = "end", type = "genomic",
                                       alt = "End Position:"),
                visual_channel_tooltip(field = tooltipField1,
                                       title = tooltipTitle,
                                       type = "quantitative",
                                       alt = paste(tooltipTitle, "Value:"),
                                       format = "0.2"),
                visual_channel_tooltip(field = tooltipField2,
                                       type = "nominal", alt = tooltipAlt2),
                visual_channel_tooltip(field = tooltipField3,
                                       type = "nominal", alt = tooltipAlt3),
                visual_channel_tooltip(field = tooltipField4,
                                       type = "nominal", alt = tooltipAlt4),
                visual_channel_tooltip(field = tooltipField5, type = "nominal",
                                       alt = tooltipAlt5)
            ),
            size = list(value = 1),
            title = title
        )
    )
}


# Create tracks for Circos visualization
.create_tracks <- function(data, gr) {
    tracks <- list()
    track_cyto <- .create_cyto_track()
    if ("gene_genomic_res" %in% unique(data$omics)) {
        gr$df_cnv$cnv_met2 <- "cnv/met"
        track_cnv <- .create_cnv_track(gr$df_cnv)
        track_met <- .create_met_track(gr$df_met)
        track_expr <- .create_expr_track(gr$df_cnv, genomic = TRUE)
        track_coef_cnv <- .create_coef_track(gr$df_cnv)
        track_coef_met <- .create_coef_track(gr$df_met)
        tracks <- c(tracks, list(
            track_cnv = track_cnv,
            track_met = track_met,
            track_expr = track_expr,
            track_coef_cnv = track_coef_cnv,
            track_coef_met = track_coef_met
        ))
    }

    if ("df_gene_cnv_res" %in% unique(names(gr))) {
        gr$df_gene_cnv_res$cnv_met <- "cnv"
        track_cnv <- .create_cnv_track(gr$df_gene_cnv_res)
        track_expr <- .create_expr_track(gr$df_gene_cnv_res)
        track_coef_cnv <- .create_coef_track(gr$df_gene_cnv_res)
        tracks <- c(tracks, list(
            track_cnv = track_cnv,
            track_expr = track_expr,
            track_coef_cnv = track_coef_cnv
        ))
    }

    if ("df_gene_met_res" %in% unique(names(gr))) {
        gr$df_gene_met_res$cnv_met <- "met"
        track_met <- .create_met_track(gr$df_gene_met_res)
        track_expr <- .create_expr_track(gr$df_gene_met_res)
        track_coef_met <- .create_coef_track(gr$df_gene_met_res)
        tracks <- c(tracks, list(
            track_met = track_met,
            track_expr = track_expr,
            track_coef_met = track_coef_met
        ))
    }

    if ("df_mirna_cnv_res" %in% unique(names(gr))) {
        gr$df_mirna_cnv_res$cnv_met <- "cnv"
        track_mirna_cnv <- .create_cnv_track(gr$df_mirna_cnv_res)
        track_mirna_expr <- .create_expr_track(gr$df_mirna_cnv_res)
        track_mirna_coef_cnv <- .create_coef_track(gr$df_mirna_cnv_res)
        tracks <- c(tracks, list(
            track_mirna_cnv = track_mirna_cnv,
            track_mirna_expr = track_mirna_expr,
            track_mirna_coef_cnv = track_mirna_coef_cnv
        ))
    }
    tracks <- c(tracks, list(track_cyto = track_cyto))
    return(tracks)
}


# Create CNV track for Circos visualization
.create_cnv_track <- function(gr) {
    gr$cov_value2 <- as.character(gr$cov_value_original)
    ccol <- .circos_track_cols(vvalues = gr$cov_value2)
    track_cnv <- .create_single_track(
        data = gr,
        dataValue = "cov_value",
        x_axis = "none",
        xe_axis = "none",
        y_axis = "none",
        colorField = "cov_value2",
        colorDomain = unique(gr$cov_value2),
        colorRange = ccol,
        tooltipField1 = "cov_value",
        tooltipTitle = "cnv",
        tooltipAlt1 = "CNV Value:",
        tooltipField2 = "gene",
        tooltipAlt2 = "Gene Name:",
        tooltipField3 = "class",
        tooltipAlt3 = "Class:",
        tooltipField4 = "cnv_met",
        tooltipAlt4 = "Integration Type:",
        tooltipField5 = "cov_value_original",
        tooltipAlt5 = "Original CNV Value:",
        legend = FALSE,
        colorType = "nominal",
        title = "CNV"
    )
    return(track_cnv)
}


# Create Met track for Circos visualization
.create_met_track <- function(gr) {
    gr$cov_value2 <- as.character(gr$cov_value_original)
    ccol <- .circos_track_cols(vvalues = gr$cov_value2)
    track_met <- .create_single_track(
        data = gr,
        dataValue = "cov_value",
        x_axis = "none",
        xe_axis = "none",
        y_axis = "none",
        colorField = "cov_value2",
        colorDomain = unique(gr$cov_value2),
        colorRange = ccol,
        tooltipField1 = "cov_value",
        tooltipTitle = "met",
        tooltipAlt1 = "MET Value:",
        tooltipField2 = "gene",
        tooltipAlt2 = "Gene Name:",
        tooltipField3 = "class",
        tooltipAlt3 = "Class:",
        tooltipField4 = "cnv_met",
        tooltipAlt4 = "Integration Type:",
        tooltipField5 = "cov_value_original",
        tooltipAlt5 = "Original MET Value:",
        legend = FALSE,
        colorType = "nominal",
        title = "MET"
    )
    return(track_met)
}

# Create Expression track for Circos visualization
.create_expr_track <- function(gr, genomic = FALSE) {
    cnv_met <- ifelse(genomic, "cnv_met2", "cnv_met")
    track_expr <- .create_single_track(
        data = gr,
        dataValue = "response_value",
        x_axis = "none",
        xe_axis = "none",
        y_axis = "none",
        colorField = "response_value",
        colorDomain = NULL,
        colorRange = NULL,
        tooltipField1 = "response_value",
        tooltipTitle = "expr",
        tooltipAlt1 = "Expression Value (log2):",
        tooltipField2 = "gene",
        tooltipAlt2 = "Gene Name:",
        tooltipField3 = "class",
        tooltipAlt3 = "Class:",
        tooltipField4 = cnv_met,
        tooltipAlt4 = "Integration Type:",
        legend = TRUE,
        colorType = "quantitative",
        title = "Expression"
    )
    return(track_expr)
}

# Create Coefficient track for Circos visualization
.create_coef_track <- function(gr) {
    gr$coef2 <- as.character(gr$coef_original)
    ccol <- .circos_track_cols(vvalues = gr$coef2)
    track_coef <- .create_single_track(
        data = gr,
        dataValue = "coef",
        x_axis = "none",
        xe_axis = "none",
        y_axis = "none",
        colorField = "coef2",
        colorDomain = unique(gr$coef2),
        colorRange = ccol,
        tooltipField1 = "coef",
        tooltipTitle = "coef",
        tooltipAlt1 = "Coef Value:",
        tooltipField2 = "gene",
        tooltipAlt2 = "Gene Name:",
        tooltipField3 = "class",
        tooltipAlt3 = "Class:",
        tooltipField4 = "cnv_met",
        tooltipAlt4 = "Integration Type:",
        tooltipField5 = "coef_original",
        tooltipAlt5 = "Original Coef Value:",
        legend = FALSE,
        colorType = "nominal",
        title = "Coefficients"
    )
    return(track_coef)
}

# Create Cytoband track for Circos visualization
#' @importFrom shiny.gosling visual_channel_stroke_width visual_channel_color
#' visual_channel_x track_data_gr add_single_track visual_channel_stroke
#'
.create_cyto_track <- function() {
    track_cyto <- add_single_track(
        id = "track2",
        data = track_data_gr(
            granges = cyto_hs,
            chromosomeField = "seqnames",
            genomicFields = c("start", "end")
        ),
        mark = "rect",
        x = visual_channel_x(
            field = "start",
            type = "genomic"
        ),
        xe = visual_channel_x(
            field = "end",
            type = "genomic"
        ),
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
                "acen" # acen: centromeric region (UCSC band files)
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


# Create composed view for Circos visualization
#' @importFrom shiny.gosling compose_view add_multi_tracks
.create_composed_view <- function(tracks, width, height) {
    composed_views <- list()

    if (sum(c("track_expr", "track_cnv", "track_met") %in% names(tracks)) ==
        3) {
        composed_view_circos_genomic <- compose_view(
            multi = TRUE,
            tracks = add_multi_tracks(
                tracks$track_cyto,
                tracks$track_expr,
                tracks$track_cnv,
                tracks$track_coef_cnv,
                tracks$track_met,
                tracks$track_coef_met
            ),
            alignment = "stack",
            spacing = 0.01,
            linkingId = "detail",
            width = width,
            height = height
        )
        composed_views <- c(composed_views,
                            list(circos_genomic = composed_view_circos_genomic))
    } else {
        if ("track_met" %in% names(tracks)) {
            composed_view_circos_met_gene <- compose_view(
                multi = TRUE,
                tracks = add_multi_tracks(
                    tracks$track_cyto,
                    tracks$track_expr,
                    tracks$track_met,
                    tracks$track_coef_met
                ),
                alignment = "stack",
                spacing = 0.01,
                linkingId = "detail"
            )
            composed_views <- c(composed_views,
                                list(circos_met_gene =
                                       composed_view_circos_met_gene))
        }
        if ("track_cnv" %in% names(tracks)) {
            composed_view_circos_cnv_gene <- compose_view(
                multi = TRUE,
                tracks = add_multi_tracks(
                    tracks$track_cyto,
                    tracks$track_expr,
                    tracks$track_cnv,
                    tracks$track_coef_cnv
                ),
                alignment = "stack",
                spacing = 0.01,
                linkingId = "detail"
            )
            composed_views <- c(composed_views,
                                list(circos_cnv_gene =
                                       composed_view_circos_cnv_gene))
        }
    }
    if (sum(c("track_mirna_expr", "track_mirna_cnv") %in% names(tracks)) == 2) {
        composed_view_circos_genomic_mirna <- compose_view(
            multi = TRUE,
            tracks = add_multi_tracks(
                tracks$track_cyto,
                tracks$track_mirna_expr,
                tracks$track_mirna_cnv,
                tracks$track_mirna_coef_cnv
            ),
            alignment = "stack",
            spacing = 0.01,
            linkingId = "detail",
            width = width,
            height = height
        )
        composed_views <- c(
            composed_views,
            list(circos_genomic_mirna = composed_view_circos_genomic_mirna)
        )
    }
    return(composed_views)
}

# Prepare Genomic Heatmap
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom dplyr arrange desc
#' @importFrom utils head
.prepare_gen_heatmap <- function(data_table,
                                 df_heatmap,
                                 df_heatmap_t,
                                 significativityCriteria,
                                 pvalRange,
                                 fdrRange,
                                 numTopCNV,
                                 numTopMET,
                                 scale,
                                 numSamples) {
    tmp <- data.frame(
        `cnv` = data_table$coef[
            data_table$cnv_met == "cnv"
        ],
        `pval_cnv` = data_table$pval[
            data_table$cnv_met == "cnv"
        ],
        `fdr_cnv` = data_table$fdr[
            data_table$cnv_met == "cnv"
        ],
        row.names = data_table$response[
            data_table$cnv_met == "cnv"
        ]
    )

    tmp2 <- data.frame(
        `met` = data_table$coef[
            data_table$cnv_met == "met"
        ],
        `pval_met` = data_table$pval[
            data_table$cnv_met == "met"
        ],
        `fdr_met` = data_table$fdr[
            data_table$cnv_met == "met"
        ],
        row.names = data_table$response[
            data_table$cnv_met == "met"
        ]
    )
    tmp3 <- unique(c(rownames(tmp), rownames(tmp2)))
    data_table <- cbind(`cnv` = tmp[tmp3, ], `met` = tmp2[tmp3, ])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp3
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t), ])
    if (significativityCriteria == "pval") {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= pvalRange |
            df_heatmap_t$pval_met <= pvalRange, ]
    } else {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= fdrRange |
            df_heatmap_t$fdr_met <= fdrRange, ]
    }
    if (nrow(df_heatmap_t) == 0) {
        return(NULL)
    }
    top_met <- df_heatmap_t %>%
        arrange(desc(abs(met))) %>%
        head(numTopCNV)
    top_cnv <- df_heatmap_t %>%
        arrange(desc(abs(cnv))) %>%
        head(numTopMET)
    expr_top <- rbind(top_cnv, top_met[
        !rownames(top_met) %in% rownames(top_cnv),
    ])
    expr_top_subset <- expr_top[, -c((ncol(expr_top) - 5):ncol(expr_top))]
    my_colors <- colorRamp2(
        c(
            min(expr_top$cnv, na.rm = TRUE),
            0, max(expr_top$cnv, na.rm = TRUE)
        ),
        c("purple", "white", "green4")
    )
    my_colors2 <- colorRamp2(
        c(
            min(expr_top$met, na.rm = TRUE),
            0, max(expr_top$met, na.rm = TRUE)
        ),
        c("purple", "white", "green4")
    )
    row_ha <- rowAnnotation(
        `coef_cnv` = expr_top$cnv, `coef_met` = expr_top$met,
        col = list(
            `coef_cnv` = my_colors,
            `coef_met` = my_colors2
        )
    )
    expr_top_subset <- as.matrix(expr_top_subset)
    if (scale == "row") expr_top_subset <- t(scale(t(log2(
      expr_top_subset + 1))))
    if (scale == "col") expr_top_subset <- scale(log2(expr_top_subset + 1))
    if (scale == "none") expr_top_subset <- log2(expr_top_subset + 1)
    if (numSamples < ncol(expr_top_subset)) {
        tmp <- apply(expr_top_subset, 2, sd)
        names(tmp) <- colnames(expr_top_subset)
        sort(tmp, decreasing = TRUE)
        ans <- names(tmp)[seq_len(numSamples)]
        expr_top_subset <- expr_top_subset[, ans]
    }
    ht <- Heatmap(expr_top_subset,
        right_annotation = row_ha,
        heatmap_legend_param = list(title = "log2 expression")
    )
    return(ht)
}

# Prepare CNV Heatmap
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom dplyr desc
#' @importFrom utils head
.prepare_cnv_heatmap <- function(data_table,
                                 df_heatmap,
                                 df_heatmap_t,
                                 significativityCriteria,
                                 pvalRange,
                                 fdrRange,
                                 numTopCNVonly,
                                 scale,
                                 numSamples) {
    tmp <- data.frame(
        `cnv` = data_table$coef[
            data_table$omics == "gene_cnv_res"
        ],
        `pval_cnv` = data_table$pval[
            data_table$omics == "gene_cnv_res"
        ],
        `fdr_cnv` = data_table$fdr[
            data_table$omics == "gene_cnv_res"
        ],
        row.names = data_table$response[
            data_table$omics == "gene_cnv_res"
        ]
    )
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(`cnv` = tmp[tmp2, ])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t), ])
    if (significativityCriteria == "pval") {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_cnv <= pvalRange, ]
    } else {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_cnv <= fdrRange, ]
    }
    if (nrow(df_heatmap_t) == 0) {
        return(NULL)
    }
    df_heatmap_t <- as.data.frame(df_heatmap_t)
    top_cnv <- df_heatmap_t %>%
        arrange(desc(abs(cnv))) %>%
        head(numTopCNVonly)
    expr_top <- top_cnv
    expr_top_subset <- expr_top[, -c((ncol(expr_top) - 2):ncol(expr_top))]
    my_colors <- colorRamp2(
        c(
            min(expr_top$cnv, na.rm = TRUE),
            0, max(expr_top$cnv, na.rm = TRUE)
        ),
        c("purple", "white", "green4")
    )
    row_ha <- rowAnnotation(
        `coef_cnv` = expr_top$cnv,
        col = list(`coef_cnv` = my_colors)
    )
    expr_top_subset <- as.matrix(expr_top_subset)
    if (scale == "row") expr_top_subset <- t(scale(t(log2(expr_top_subset +
                                                            1))))
    if (scale == "col") expr_top_subset <- scale(log2(expr_top_subset + 1))
    if (scale == "none") expr_top_subset <- log2(expr_top_subset + 1)
    if (numSamples < ncol(expr_top_subset)) {
        tmp <- apply(expr_top_subset, 2, sd)
        names(tmp) <- colnames(expr_top_subset)
        sort(tmp, decreasing = TRUE)
        ans <- names(tmp)[seq_len(numSamples)]
        expr_top_subset <- expr_top_subset[, ans]
    }
    ht <- Heatmap(expr_top_subset,
        right_annotation = row_ha,
        heatmap_legend_param = list(title = "log2 expression")
    )
    return(ht)
}

# Prepare Met Heatmap
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom dplyr desc
#' @importFrom utils head
.prepare_met_heatmap <- function(data_table,
                                 df_heatmap,
                                 df_heatmap_t,
                                 significativityCriteria,
                                 pvalRange,
                                 fdrRange,
                                 numTopMETonly,
                                 scale,
                                 numSamples) {
    tmp <- data.frame(
        `met` = data_table$coef[
            data_table$omics == "gene_met_res"
        ],
        `pval_met` = data_table$pval[
            data_table$omics == "gene_met_res"
        ],
        `fdr_met` = data_table$fdr[
            data_table$omics == "gene_met_res"
        ],
        row.names = data_table$response[
            data_table$omics == "gene_met_res"
        ]
    )
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(`met` = tmp[tmp2, ])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t), ])
    if (significativityCriteria == "pval") {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_met <= pvalRange, ]
    } else {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_met <= fdrRange, ]
    }
    if (nrow(df_heatmap_t) == 0) {
        return(NULL)
    }
    df_heatmap_t <- as.data.frame(df_heatmap_t)
    top_met <- df_heatmap_t %>%
        arrange(desc(abs(met))) %>%
        head(numTopMETonly)
    expr_top <- top_met
    expr_top_subset <- expr_top[, -c((ncol(expr_top) - 2):ncol(expr_top))]
    my_colors <- colorRamp2(
        c(
            min(expr_top$met, na.rm = TRUE),
            0, max(expr_top$met, na.rm = TRUE)
        ),
        c("purple", "white", "green4")
    )
    row_ha <- rowAnnotation(
        `coef_met` = expr_top$met,
        col = list(`coef_met` = my_colors)
    )
    expr_top_subset <- as.matrix(expr_top_subset)
    if (scale == "row") expr_top_subset <- t(scale(t(log2(expr_top_subset +
                                                            1))))
    if (scale == "col") expr_top_subset <- scale(log2(expr_top_subset + 1))
    if (scale == "none") expr_top_subset <- log2(expr_top_subset + 1)
    if (numSamples < ncol(expr_top_subset)) {
        tmp <- apply(expr_top_subset, 2, sd)
        names(tmp) <- colnames(expr_top_subset)
        sort(tmp, decreasing = TRUE)
        ans <- names(tmp)[seq_len(numSamples)]
        expr_top_subset <- expr_top_subset[, ans]
    }
    ht <- Heatmap(expr_top_subset,
        right_annotation = row_ha,
        heatmap_legend_param = list(title = "log2 expression")
    )
    return(ht)
}

# Prepare miRNA Heatmap
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom dplyr desc
#' @importFrom utils head
.prepare_mirna_heatmap <- function(data_table,
                                   df_heatmap,
                                   df_heatmap_t,
                                   significativityCriteria,
                                   pvalRange,
                                   fdrRange,
                                   numTopMiCNV,
                                   scale,
                                   numSamples) {
    tmp <- data.frame(
        `mirna_cnv` = data_table$coef[
            data_table$omics == "mirna_cnv_res"
        ],
        `pval_mirna_cnv` = data_table$pval[
            data_table$omics == "mirna_cnv_res"
        ],
        `fdr_mirna_cnv` = data_table$fdr[
            data_table$omics == "mirna_cnv_res"
        ],
        row.names = data_table$response[
            data_table$omics == "mirna_cnv_res"
        ]
    )
    tmp2 <- unique(rownames(tmp))
    data_table <- cbind(`mirna_cnv` = tmp[tmp2, ])
    colnames(data_table) <- gsub("^.*\\.", "", colnames(data_table))
    rownames(data_table) <- tmp2
    df_heatmap_t <- cbind(df_heatmap_t, data_table[rownames(df_heatmap_t), ])
    if (significativityCriteria == "pval") {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$pval_mirna_cnv <= pvalRange, ]
    } else {
        df_heatmap_t <- df_heatmap_t[df_heatmap_t$fdr_mirna_cnv <= fdrRange, ]
    }
    if (nrow(df_heatmap_t) == 0) {
        return(NULL)
    }
    df_heatmap_t <- as.data.frame(df_heatmap_t)
    top_mirna_cnv <- df_heatmap_t %>%
        arrange(desc(abs(mirna_cnv))) %>%
        head(numTopMiCNV)
    expr_top <- top_mirna_cnv
    expr_top_subset <- expr_top[, -c((ncol(expr_top) - 2):ncol(expr_top))]
    my_colors <- colorRamp2(
        c(
            min(expr_top$mirna_cnv, na.rm = TRUE),
            0, max(expr_top$mirna_cnv, na.rm = TRUE)
        ),
        c("purple", "white", "green4")
    )
    row_ha <- rowAnnotation(
        `coef_mirna` = expr_top$mirna_cnv,
        col = list(`coef_mirna` = my_colors)
    )
    expr_top_subset <- as.matrix(expr_top_subset)
    if (scale == "row") expr_top_subset <- t(scale(t(log2(expr_top_subset +
                                                            1))))
    if (scale == "col") expr_top_subset <- scale(log2(expr_top_subset + 1))
    if (scale == "none") expr_top_subset <- log2(expr_top_subset + 1)
    if (numSamples < ncol(expr_top_subset)) {
        tmp <- apply(expr_top_subset, 2, sd)
        names(tmp) <- colnames(expr_top_subset)
        sort(tmp, decreasing = TRUE)
        ans <- names(tmp)[seq_len(numSamples)]
        expr_top_subset <- expr_top_subset[, ans]
    }
    ht <- Heatmap(expr_top_subset,
        right_annotation = row_ha,
        heatmap_legend_param = list(title = "log2 expression")
    )
    return(ht)
}

# Run function in background
#' @importFrom callr r_bg
#'
.run_bg <- function(FFUN,
                    args,
                    input,
                    output) {
    ans <- r_bg(
        func = FFUN,
        args = args,
        supervise = TRUE,
        package = "gINTomics"
    )
    return(ans)
}


#' Start a Shiny application for integrated multi-omics data analysis.
#'
#' The \code{run_shiny} function launches an interactive Shiny application that
#' allows users to explore and analyze integrated multi-omics data through
#' various visualizations and analyses.
#'
#' @param multiomics_integration An object representing the integration of
#' multi-omics data, compatible with the \code{\link{extract_model_res}}
#' function.
#'
#' @details The \code{run_shiny} function extracts model results from
#' \code{multiomics_integration}, performs preprocessing operations to prepare
#' the data for the Shiny user interface, creates the user interface and server
#' for the Shiny application.
#'
#' @return No return value. The function starts an interactive Shiny
#' application.
#'
#' @export
#' @importFrom shiny NS callModule shinyApp
#'
#' @examples
#' # Example usage:
#' library(MultiAssayExperiment)
#' data("mmultiassay_ov")
#' tmp <- lapply(experiments(mmultiassay_ov), function(x) x[1:20,])
#' mmultiassay_ov <- MultiAssayExperiment(experiments = tmp)
#' # multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' # app <- run_shiny(multiomics_integration)
#'
#' @seealso
#' \code{\link{extract_model_res}}
#'
#' @references
#' Description of the multi-omics data model and integrated analysis techniques
#' used.
#'
#' @keywords shiny multiomics integration visualization analysis
#' @family Shiny
#' @family Data analysis
#' @family Multi-omics
#' @family Integration
#' @family Visualization
#' @family Interactive
#' @family Function

run_shiny <- function(multiomics_integration) {
    data <- extract_model_res(multiomics_integration)
    data <- .shiny_preprocess(data)
    data_table <- data$data_table
    ui <- .create_ui(data_table)
    server <- function(input, output, session) {
        # ---------------------- NETWORK SERVERS -----------------------------
        callModule(.server_network,
            id = "network_trans",
            data_table = data_table
        )
        callModule(.server_network,
            id = "network_deg",
            data_table = data_table,
            deg = TRUE
        )

        ### ------------------------ VENN SERVERS ----------------------
        callModule(.server_venn,
            id = "venn_gen",
            data_table = data_table
        )
        callModule(.server_venn,
            id = "venn_deg",
            data_table = data_table,
            deg = TRUE
        )

        ## -------------------------- VOLCANO SERVERS ------------------------
        callModule(.server_volcano,
            id = "volcano_gen",
            data_table = data_table
        )
        callModule(.server_volcano,
            id = "volcano_trans",
            data_table = data_table
        )
        callModule(.server_volcano,
            id = "volcano_deg",
            data_table = data_table,
            deg = TRUE
        )

        ## -------------------------- HEATMAP SERVER ------------------------
        .prepare_reactive_heatmap(
            data_table = data_table,
            multiomics_integration = multiomics_integration,
            input = input,
            output = output,
            session = session,
            deg = FALSE,
            ns = NS("heat_gen")
        )
        .prepare_reactive_heatmap(
            data_table = data_table,
            multiomics_integration = multiomics_integration,
            input = input,
            output = output,
            session = session,
            deg = TRUE,
            ns = NS("heat_deg")
        )
        ## ---------------------- RIDGE SERVERS ------------------------
        callModule(.server_ridge,
            id = "ridge_gen",
            data_table = data_table
        )
        callModule(.server_ridge,
            id = "ridge_trans",
            data_table = data_table
        )
        callModule(.server_ridge,
            id = "ridge_deg",
            data_table = data_table,
            deg = TRUE
        )

        ## ----------------------- HISTO SERVERS --------------------------
        callModule(.server_histo,
            id = "histo_gen",
            data_table = data_table
        )
        callModule(.server_histo2,
            id = "histo_trans",
            data_table = data_table
        )
        callModule(.server_histo2,
            id = "histo_deg",
            data_table = data_table,
            deg = TRUE
        )

        #### ------------------- TABLE SERVER ----------------------------
        callModule(.server_table,
            id = "complete_table",
            data_table = data_table
        )

        #### ------------------- ENRICHMENT SERVER ----------------------------
        data_gen_enrich <- data_table[data_table$omics %in% c(
            "gene_genomic_res",
            "gene_cnv_res",
            "gene_met_res"
        ), ]
        data_tf_enrich <- data_table[data_table$omics == "tf_res", ]
        callModule(.server_enrich_bg,
            id = "enrich_gen", name = "enrich_gen",
            extracted_data = data_gen_enrich
        )
        callModule(.server_enrich_bg,
            id = "enrich_tf", name = "enrich_tf",
            extracted_data = data_tf_enrich, tf = TRUE
        )
        #### ------------------- CIRCOS SERVER ----------------------------
        callModule(.server_circos,
            id = "circos",
            data = data$data
        )
    }
    shinyApp(ui = ui, server = server)

}
