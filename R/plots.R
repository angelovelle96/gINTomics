#' Plotting network
#' @param data_table The data table containing network information.
#' @param num_interactions The number of interactions to display in the network
#' (default: 300).
#' @param class Optional. The class of interactions to include in the plot.
#' @param pval The p-value threshold for selecting interactions (default: 0.05).
#' @return A network plot.
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' # plot_network(data_table)
#' @export
plot_network <- function(data_table,
                         num_interactions = 300,
                         class = NULL,
                         pval = 0.05) {
    if (!is.null(class)) data_table <- data_table[data_table$class == class, ]
    nnet <- .prepare_network(data_table)
    nnet$edges <- nnet$edges[nnet$edges$pval <= pval, ]
    nnet$edges <- nnet$edges[seq_len(length.out = num_interactions), ]
    nodes_with_edges <- unique(c(nnet$edges$from, nnet$edges$to))
    nnet$nodes <- nnet$nodes[nnet$nodes$id %in% nodes_with_edges, ]
    .build_network(
        nodes = nnet$nodes,
        edges = nnet$edges,
        legend_nodes = nnet$legend_nodes,
        legend_edges = nnet$legend_edges
    )
}


#' plotting venn
#' @param data_table The data table containing information for the Venn diagram.
#' @param class Optional. The class of interactions to include in the Venn
#' diagram.
#' @return A Venn diagram plot.
#' @importFrom shiny isolate
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' plot_venn(data_table, omics = "gene_genomic_res", cnv_met = "cnv")
#' @export
plot_venn <- function(data_table,
                      class = NULL) {
    if (is.null(class) & "class" %in% colnames(data_table)) {
        stop(str_wrap("Please, specify class"))
    }
    input <- list()
    output <- list()
    input$ClassSelect <- class
    input$FdrRange <- c(0, 0.05)
    input$PvalRange <- c(0, 0.05)
    input$SignificativityCriteria <- "pval"
    reactive_venn <- .prepare_reactive_venn(
        data_table = data_table,
        input = input,
        output = output,
        deg = FALSE
    )
    tmp <- isolate(reactive_venn())
    ans <- .build_venn(tmp)
    ggplotly(ans)
}


#' plotting volcano
#' @param data_table The data table containing information for the volcano plot.
#' @param class Optional. The class of interactions to include in the volcano
#' plot.
#' @param omics Optional. The omics type for the volcano plot.
#' @param cnv_met Optional. Indicates whether the volcano plot is for CNV or
#' MET omics (only applicable if omics is specified).
#' @return A volcano plot.
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' plot_volcano(data_table, omics = "gene_genomic_res", cnv_met = "cnv")
#' @export
plot_volcano <- function(data_table,
                         class = NULL,
                         omics = NULL,
                         cnv_met = NULL) {
    if (!is.null(class) & !"class" %in% colnames(data_table)) {
        stop(str_wrap("Class not available"))
    }
    if (is.null(omics) & length(unique(data_table$omics)) > 0) {
        stop(str_wrap("Please, specify omics"))
    }

    if (!is.null(class)) data_table <- data_table[data_table$class == class, ]
    if (!is.null(omics)) data_table <- data_table[data_table$omics == omics, ]
    if (omics == "gene_genomic_res" & !is.null(cnv_met)) {
        if (!cnv_met %in% c("cnv", "met")) {
            stop(str_wrap("cnv_met should be one of cnv and met"))
        }
        data_table <- data_table[data_table$cnv_met == cnv_met, ]
        data_table <- data_table[!is.na(data_table$cnv_met), ]
    }
    data_table["group"] <- "Not Significant"
    data_table[data_table$pval <= 0.05, "group"] <- "Significant"
    data_table$pval_fdr <- -log10(data_table$pval)
    if (!"class" %in% colnames(data_table)) data_table$class <- data_table$omics
    .build_volcano(data_table)
}


#' plotting ridge
#' @param data_table The data table containing information for the ridge plot.
#' @param class Optional. The class of interactions to include in the ridge
#' plot.
#' @param omics Optional. The omics type for the ridge plot.
#' @param cnv_met Optional. Indicates whether the ridge plot is for CNV or MET
#' omics (only applicable if omics is specified).
#' @return A ridge plot.
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' plot_ridge(data_table, omics = "gene_genomic_res", cnv_met = "cnv")
#' @export
plot_ridge <- function(data_table,
                       class = NULL,
                       omics = NULL,
                       cnv_met = NULL) {
    if (is.null(class) & "class" %in% colnames(data_table)) {
        stop(str_wrap("Please, specify class"))
    }
    if (is.null(omics) & length(unique(data_table$omics)) > 0) {
        stop(str_wrap("Please, specify omics"))
    }

    if (!is.null(class)) data_table <- data_table[data_table$class == class, ]
    if (!is.null(omics)) data_table <- data_table[data_table$omics == omics, ]
    data_table$significance <- ifelse(data_table$pval <= 0.05,
        "Significant", "Not Significant"
    )
    if (omics == "gene_genomic_res" & !is.null(cnv_met)) {
        if (!cnv_met %in% c("cnv", "met")) {
            stop(str_wrap("cnv_met should be one of cnv and met"))
        }
        data_table <- data_table[data_table$cnv_met == cnv_met, ]
        data_table <- data_table[!is.na(data_table$cnv_met), ]
    }
    lower_quantile <- quantile(data_table$coef, 0.001)
    upper_quantile <- quantile(data_table$coef, 0.999)
    quantiles <- c(lower_quantile, upper_quantile)
    .build_ridge(ridge_data = data_table, quantiles = quantiles)
}


#' plotting heatmap
#' @param multiomics_integration The multiomics integration object.
#' @param data_table The data table containing information for the heatmap.
#' @param omics The type of omics data for the heatmap.
#' @param scale Optional. The scale type for the heatmap. Default is "none".
#' @param genes_number Optional. The number of genes to include in the heatmap.
#' Default is 50.
#' @param class Optional. The class of interactions to include in the heatmap.
#' @param pval Optional. The p-value threshold for significance in the heatmap.
#' Default is 0.05.
#' @param samples_number Number of samples to include in the heatmap. If this
#' number is inferior to the total number of samples, the n most variable
#' samples will be selected
#' @return A heatmap plot.
#' @importFrom methods is
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' plot_heatmap(data_table, omics = "gene_genomic_res")
#' @export
plot_heatmap <- function(multiomics_integration,
                         data_table,
                         omics,
                         scale = "none",
                         genes_number = 50,
                         samples_number = 50,
                         class = NULL,
                         pval = 0.05) {
    if (!"omics" %in% colnames(data_table)) data_table$omics <- omics
    if (is.null(class) & "class" %in% colnames(data_table)) {
        stop(str_wrap("Please, specify class"))
    }
    tmp <- c("gene_genomic_res", "gene_cnv_res", "gene_met_res",
             "mirna_cnv_res")
    if (!omics %in% tmp) {
        stop(paste(
            "Omics should be one of",
            paste(tmp, collapse = ", ")
        ))
    }
    if (!is(multiomics_integration, "MultiOmics")) {
        tmp <- list()
        tmp[[omics]] <- multiomics_integration
        multiomics_integration <- tmp
    }

    if ("class" %in% colnames(data_table)) {
        df_heatmap <- multiomics_integration[[omics]][[
            class
        ]]$data$response_var
        data_table <- data_table[data_table$omics == omics, ]
        data_table <- data_table[data_table$class == class, ]
    } else {
        df_heatmap <- multiomics_integration[[
            omics
        ]]$data$response_var
        data_table <- data_table[data_table$omics == omics, ]
    }

    df_heatmap_t <- t(as.matrix(df_heatmap))
    if (omics == "gene_genomic_res") {
        ans <- .prepare_gen_heatmap(
            data_table = data_table,
            df_heatmap = df_heatmap,
            df_heatmap_t = df_heatmap_t,
            significativityCriteria = "pval",
            pvalRange = pval,
            fdrRange = pval,
            numTopCNV = round(genes_number / 2),
            numTopMET = round(genes_number / 2),
            scale = scale,
            numSamples = samples_number
        )
    }
    if (omics == "gene_cnv_res") {
        ans <- .prepare_cnv_heatmap(
            data_table = data_table,
            df_heatmap = df_heatmap,
            df_heatmap_t = df_heatmap_t,
            significativityCriteria = "pval",
            pvalRange = pval,
            fdrRange = pval,
            numTopCNVonly = genes_number,
            scale = scale,
            numSamples = samples_number
        )
    }
    if (omics == "gene_met_res") {
        ans <- .prepare_met_heatmap(
            data_table = data_table,
            df_heatmap = df_heatmap,
            df_heatmap_t = df_heatmap_t,
            significativityCriteria = "pval",
            pvalRange = pval,
            fdrRange = pval,
            numTopMETonly = genes_number,
            scale = scale,
            numSamples = samples_number
        )
    }
    if (omics == "mirna_cnv_res") {
        ans <- .prepare_mirna_heatmap(
            data_table = data_table,
            df_heatmap = df_heatmap,
            df_heatmap_t = df_heatmap_t,
            significativityCriteria = "pval",
            pvalRange = pval,
            fdrRange = pval,
            numTopMiCNV = genes_number,
            scale = scale,
            numSamples = samples_number
        )
    }
    return(ans)
}


#' plotting chr distribution
#' @param data_table The data table containing information for plotting
#' chromosome distribution.
#' @param class Optional. The class of interactions to include in the plot.
#' @param omics Optional. The type of omics data for the plot.
#' @param cnv_met Optional. The type of copy number variation or methylation
#' data.
#' @param pval Optional. The p-value threshold for significance. Default is
#' 0.05.
#' @return A histogram plot showing chromosome distribution.
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' plot_chr_distribution(data_table, omics = "gene_genomic_res")
#' @export
plot_chr_distribution <- function(data_table,
                                  class = NULL,
                                  omics = NULL,
                                  cnv_met = NULL,
                                  pval = 0.05) {
    if (is.null(class) & "class" %in% colnames(data_table)) {
        stop(str_wrap("Please, specify class"))
    }
    if (is.null(omics) & length(unique(data_table$omics)) > 0) {
        stop(str_wrap("Please, specify omics"))
    }
    data_table <- data_table[!is.na(data_table$chr_cov), ]
    data_table <- data_table[!is.na(data_table$chr_response), ]
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)

    if (!is.null(class) & "class" %in% colnames(data_table)) {
        data_table <- data_table[data_table$class == class, ]
    }
    if (!is.null(omics) & length(unique(data_table$omics)) > 0) {
        data_table <- data_table[data_table$omics == omics, ]
    }
    if (omics == "gene_genomic_res" & !is.null(cnv_met)) {
        if (!cnv_met %in% c("cnv", "met")) {
            stop(str_wrap("cnv_met should be one of cnv and met"))
        }
        data_table <- data_table[data_table$cnv_met == cnv_met, ]
        data_table <- data_table[!is.na(data_table$cnv_met), ]
    }
    data_table$significance <- ifelse(data_table$pval <= pval,
        "Significant",
        "Not Significant"
    )
    if (nrow(data_table) == 0) {
        return(NULL)
    }
    .build_histo(data_table)
}

#' plotting TF distribution
#' @param data_table The data table containing TF information.
#' @param class Optional. The class of interactions to include in the
#' distribution plot.
#' @param pval Optional. The p-value threshold for significance in the
#' distribution plot. Default is 0.05.
#' @return A TF distribution plot.
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' data_table <- extract_model_res(multiomics_integration)
#' plot_tf_distribution(data_table)
#' @export
plot_tf_distribution <- function(data_table,
                                 class = NULL,
                                 pval = 0.05) {
    omics <- "tf_res"
    if (is.null(class) & "class" %in% colnames(data_table)) {
        stop(str_wrap("Please, specify class"))
    }
    chr_order <- mixedsort(unique(data_table$chr_cov))
    chr_order <- chr_order[!is.na(chr_order)]
    data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
    if ("omics" %in% colnames(data_table)) {
        data_table <- data_table[data_table$omics == omics, ]
    }
    data_table <- data_table[
        ,
        colnames(data_table) %in% c(
            "response", "cov",
            "pval", "fdr", "chr_cov",
            "deg", "class"
        )
    ]
    if (!is.null(class) & "class" %in% colnames(data_table)) {
        data_table <- data_table[data_table$class == class, ]
    }
    df_filtered_histo_tf <- data_table
    df_filtered_histo_tf <- df_filtered_histo_tf[
        df_filtered_histo_tf$pval <= pval,
    ]
    genes_count <- table(df_filtered_histo_tf$cov, df_filtered_histo_tf$chr_cov)
    genes_count_df <- as.data.frame.table(genes_count)
    genes_count_df <- subset(genes_count_df, `Freq` != 0)
    colnames(genes_count_df) <- c("TF", "Chromosome", "Count")
    genes_count_df <- genes_count_df[order(-genes_count_df$Count), ]
    .build_histo_TFbyChr(genes_count_df)
}


#' plotting enrichment
#' @param enrich_result Enrichment analysis results.
#' @param title Title of the plot.
#' @param showCategory Number of categories to display.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @return A plotly object containing the dot plot.
#' @importFrom ggtree fortify
#' @importFrom plotly add_markers subplot plot_ly
#' @examples
#' # Example usage:
#' data("ov_test_tcga_omics")
#' # multiomics_integration <- run_multiomics(data = mmultiassay_ov)
#' # gen_enr <- run_genomic_enrich(multiomics_integration, qvalueCutoff = 1,
#'  pvalueCutoff = 0.05, pAdjustMethod = "none")
#' # dot_plotly(gen_enr, title = "Enrichment Analysis",showCategory = 10)
#' @export
dot_plotly <- function(enrich_result,
                       title = NULL,
                       showCategory = 10,
                       width = 800,
                       height = 700) {
    df <- fortify(enrich_result, showCategory = showCategory)
    if (nrow(df) == 0) {
        return(NULL)
    }
    df$Description <- as.character(df$Description)
    df <- df[order(df$GeneRatio, decreasing = TRUE), ]
    df$Description <- unlist(lapply(df$Description, function(label) {
        words <- strsplit(label, " ")[[1]]
        split_words <- lapply(seq(1, length(words), by = 2), function(i) {
            paste(words[i:min(i + 1, length(words))], collapse = " ")
        })
        paste(split_words, collapse = "<br>")
    }))
    if (length(df$Count[!is.na(df$Count)]) > 0) {
        legend.sizes <- seq(
            min(df$Count, na.rm = TRUE),
            max(df$Count, na.rm = TRUE),
            max(c(1, round(((max(df$Count, na.rm = TRUE) -
                min(df$Count, na.rm = TRUE)) / 4))))
        )
    } else {
        return(NULL)
    }
    lprop <- c(0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55)
    lprop <- lprop[length(legend.sizes)]
    ax <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = TRUE,
        showgrid = FALSE,
        side = "top"
    )
    mk <- list(sizeref = 0.1, sizemode = "area")

    pplot <- plot_ly(df, width = width, height = height) %>%
        add_markers(
            x = ~GeneRatio,
            y = ~Description,
            name = "DotPlot",
            text = ~ paste(
                "Count:", Count, "<br>",
                "pValue:", round(pvalue, digits = 4), "<br>",
                "qValue:", round(qvalue, digits = 4)
            ),
            size = ~Count,
            color = ~pvalue,
            marker = mk,
            type = "scatter",
            mode = "markers"
        ) %>%
        layout(
            yaxis = list(
                automargin = TRUE,
                tickfont = list(size = 7),
                categoryorder = "array",
                categoryarray = rev(df$Description)
            ),
            xaxis = list(title = "GeneRatio"),
            title = title
        )

    llegend <- plot_ly() %>%
        add_markers(
            x = "Count",
            y = legend.sizes,
            size = legend.sizes,
            showlegend = FALSE,
            fill = ~"",
            marker = c(mk, color = "black"),
            text = legend.sizes,
            hoverinfo = "text"
        ) %>%
        layout(
            xaxis = ax,
            yaxis = list(showgrid = FALSE, tickvals = legend.sizes)
        )

    empty_trace <- plot_ly(
        x = numeric(0),
        y = numeric(0),
        type = "scatter",
        mode = "markers",
        showlegend = FALSE
    ) %>%
        layout(
            xaxis = list(
                showgrid = FALSE,
                zeroline = FALSE,
                showticklabels = FALSE
            ),
            yaxis = list(
                showgrid = FALSE,
                zeroline = FALSE,
                showticklabels = FALSE
            )
        )

    ans <- subplot(empty_trace, llegend, empty_trace,
        heights = c(0.02, lprop, (0.98 - lprop)), nrows = 3
    )
    ans <- subplot(pplot, ans,
        widths = c(0.9, 0.1), titleX = TRUE, shareX = FALSE
    )
    return(ans)
}
