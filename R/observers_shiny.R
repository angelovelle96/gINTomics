# Render Reactive Network
#' @importFrom visNetwork renderVisNetwork
#' @importFrom visNetwork visHierarchicalLayout
#' @importFrom visNetwork visIgraphLayout
#' @importFrom dplyr %>%
#' @importFrom shiny observe bindEvent
.render_reactive_network <- function(reactive_network,
                                     input,
                                     output) {
    observe({
        output$networkPlot <- renderVisNetwork({
            network_data <- reactive_network()
            if (nrow(network_data$nodes) == 0) {
                return(NULL)
            }
            network <- .build_network(
                nodes = network_data$nodes,
                edges = network_data$edges,
                legend_nodes = network_data$legend_nodes,
                legend_edges = network_data$legend_edges
            )
            network <- network %>%
                visHierarchicalLayout(enabled = input$layout) %>%
                visIgraphLayout(
                    physics = input$physics,
                    layout = "layout_with_fr",
                    randomSeed = 20
                )

            return(network)
        })
    }) %>% bindEvent(
        input$layout,
        input$physics
    )
}

# Select Network Data
#' @importFrom shiny reactive bindEvent
#' @importFrom dplyr %>%
.select_network <- function(data_table,
                            network_data,
                            input,
                            output,
                            deg = FALSE) {
    reactive({
        data_table <- data_table[data_table$omics %in% c(
            "tf_res",
            "tf_mirna_res",
            "mirna_target_res"
        ), ]

        numInteractions <- input$numInteractions
        significativityCriteria <- input$SignificativityCriteria
        pval <- input$PvalRange
        fdr <- input$FdrRange
        class <- input$ClassSelect
        ans <- network_data
        if ("class" %in% colnames(ans$edges)) {
            ans$edges <- ans$edges[ans$edges$class == class, ]
        }
        nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
        ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges, ]
        if (significativityCriteria == "pval") {
            ans$edges <- ans$edges[ans$edges$pval <= pval, ]
        } else {
            ans$edges <- ans$edges[ans$edges$fdr <= fdr, ]
        }
        if (deg) {
            tto <- unique(data_table$response[data_table$deg])
            ans$edges <- ans$edges[ans$edges$to %in% tto, ]
        }
        ans$edges <- ans$edges[seq_len(length.out = numInteractions), ]
        nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
        ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges, ]
        return(ans)
    }) %>% bindEvent(
        input$layout,
        input$numInteractions,
        input$SignificativityCriteria,
        input$PvalRange,
        input$fdrNetwork,
        input$FdrRange,
        input$ClassSelect
    )
}

# Prepare Reactive Venn
#' @importFrom dplyr %>%
#' @importFrom shiny reactive bindEvent
.prepare_reactive_venn <- function(data_table,
                                   input,
                                   output,
                                   deg = FALSE) {
    reactive_venn <- reactive({
        data_venn <- NULL
        classSelect <- input$ClassSelect
        pvalRange <- input$PvalRange
        fdrRange <- input$FdrRange
        significativityCriteria <- input$SignificativityCriteria
        if (deg == FALSE & "gene_genomic_res" %in% unique(data_table$omics)) {
            data_table <- data_table[data_table$omics == "gene_genomic_res", ]
            if ("class" %in% colnames(data_table)) {
                data_table <- data_table[data_table$class == classSelect, ]
            }
            if (significativityCriteria == "pval") {
                cnv_sign_genes <- subset(data_table,
                    `cnv_met` == "cnv" & `pval` >= pvalRange[1] &
                        `pval` <= pvalRange[2],
                    select = c("cov")
                )
                met_sign_genes <- subset(data_table,
                    `cnv_met` == "met" & `pval` >= pvalRange[1] &
                        `pval` <= pvalRange[2],
                    select = c("cov")
                )
            } else {
                cnv_sign_genes <- subset(data_table,
                    `cnv_met` == "cnv" & `fdr` >= fdrRange[1] &
                        `fdr` <= fdrRange[2],
                    select = c("cov")
                )
                met_sign_genes <- subset(data_table,
                    `cnv_met` == "met" & `fdr` >= fdrRange[1] &
                        `fdr` <= fdrRange[2],
                    select = c("cov")
                )
            }
            data_venn <- list(
                `cnv_sign_genes` = cnv_sign_genes,
                `met_sign_genes` = met_sign_genes
            )
        }
        if (deg == FALSE & !"gene_genomic_res" %in% unique(data_table$omics)) {
            if ("gene_cnv_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_cnv_res", ]
            }
            if ("gene_met_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_met_res", ]
            }
        }
        if (deg == TRUE & "gene_genomic_res" %in% unique(data_table$omics)) {
            if ("class" %in% colnames(data_table)) {
                data_table <- data_table[data_table$class == classSelect, ]
                data_table <- data_table[data_table$deg, ]
                if (significativityCriteria == "pval") {
                    cnv_sign_genes <- subset(data_table,
                        cnv_met == "cnv" & pval >= pvalRange[1] &
                            pval <= pvalRange[2],
                        select = c("cov")
                    )
                    met_sign_genes <- subset(data_table,
                        cnv_met == "met" & pval >= pvalRange[1] &
                            pval <= pvalRange[2],
                        select = c("cov")
                    )
                } else {
                    cnv_sign_genes <- subset(data_table,
                        cnv_met == "cnv" & fdr >= fdrRange[1] &
                            fdr <= fdrRange[2],
                        select = c("cov")
                    )
                    met_sign_genes <- subset(data_table,
                        cnv_met == "met" & fdr >= fdrRange[1] &
                            fdr <= fdrRange[2],
                        select = c("cov")
                    )
                }
                data_venn <- list(
                    cnv_sign_genes = cnv_sign_genes,
                    met_sign_genes = met_sign_genes
                )
            } else {
                (return(NULL))
            }
        }
        if (deg == TRUE & !"gene_genomic_res" %in% unique(data_table$omics)) {
            if ("gene_cnv_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_cnv_res", ]
            }
            if ("gene_met_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_met_res", ]
            }
        }
        return(data_venn)
    }) %>% bindEvent(
        input$ClassSelect,
        input$PvalRange,
        input$FdrRange,
        input$SignificativityCriteria
    )
}

# Prepare Reactive Volcano
#' @importFrom dplyr %>%
#' @importFrom shiny reactive bindEvent
.prepare_reactive_volcano <- function(data_table,
                                      input,
                                      output,
                                      deg = FALSE) {
    reactive({
        integrationSelect <- input$IntegrationSelect
        typeSelect <- input$genomicTypeSelect
        significativityCriteria <- input$SignificativityCriteria
        pvalRange <- input$PvalRange
        fdrRange <- input$FdrRange
        if (deg) data_table <- data_table[data_table$deg, ]
        if (!"class" %in% colnames(data_table)) data_table$class <- data_table$
          omics
        data_table <- data_table[data_table$omics == integrationSelect, ]
        if (integrationSelect == "gene_genomic_res") {
            data_table <- data_table[data_table$cnv_met == typeSelect, ]
        }
        if (nrow(data_table) == 0) {
            return(NULL)
        }
        if (significativityCriteria == "pval") {
            data_table["group"] <- "Not Significant"
            data_table[data_table$pval <= pvalRange, "group"] <- "Significant"
            data_table$pval_fdr <- -log10(data_table$pval)
        } else {
            data_table["group"] <- "Not Significant"
            data_table[data_table$fdr <= fdrRange, "group"] <- "Significant"
            data_table$pval_fdr <- -log10(data_table$fdr)
        }
        return(data_table)
    }) %>% bindEvent(
        input$IntegrationSelect,
        input$genomicTypeSelect,
        input$SignificativityCriteria,
        input$PvalRange,
        input$FdrRange
    )
}

# Prepare Reactive Heatmap
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom ComplexHeatmap ht_opt
#' @importFrom shiny observe bindEvent
#' @importFrom dplyr %>%
.prepare_reactive_heatmap <- function(data_table,
                                      multiomics_integration,
                                      input,
                                      output,
                                      session,
                                      ns,
                                      deg = FALSE) {
    observe({
        ht_opt$message <- FALSE
        integrationSelect <- input[[ns("IntegrationSelect")]]
        numTopCNV <- input[[ns("numTopGenesHeatmapCNV")]]
        numTopMET <- input[[ns("numTopGenesHeatmapMET")]]
        numTopCNVonly <- input[[ns("numTopGenesHeatmapCNVonly")]]
        numTopMETonly <- input[[ns("numTopGenesHeatmapMETonly")]]
        numTopMiCNV <- input[[ns("numTopGenesHeatmapmirna_cnv")]]
        numSamples <- input[[ns("numSamples")]]
        classSelect <- input[[ns("ClassSelect")]]
        significativityCriteria <- input[[ns("SignificativityCriteria")]]
        pvalRange <- input[[ns("PvalRange")]]
        fdrRange <- input[[ns("FdrRange")]]
        scale <- input[[ns("scaleHeatmap")]]
        if ("class" %in% colnames(data_table) & deg == FALSE) {
            df_heatmap <- multiomics_integration[[integrationSelect]][[
                classSelect
            ]]$data$response_var
            data_table <- data_table[data_table$omics == integrationSelect, ]
            data_table <- data_table[data_table$class == classSelect, ]
        }
        if ("class" %in% colnames(data_table) & deg == TRUE) {
            df_heatmap <- multiomics_integration[[integrationSelect]][[
                classSelect
            ]]$data$response_var
            data_table <- data_table[data_table$omics == integrationSelect, ]
            data_table <- data_table[data_table$class == classSelect, ]
            data_table <- data_table[data_table$deg, ]
            df_heatmap <- df_heatmap[, unique(data_table$response)]
        }
        if (!"class" %in% colnames(data_table)) {
            df_heatmap <- multiomics_integration[[
                integrationSelect
            ]]$data$response_var
            data_table <- data_table[data_table$omics == integrationSelect, ]
        }
        if (is.null(df_heatmap)) {
            return(NULL)
        }
        if (!"class" %in% colnames(data_table) & deg == TRUE) {
            return(NULL)
        }
        df_heatmap_t <- t(as.matrix(df_heatmap))
        if (integrationSelect == "gene_genomic_res") {
            ans <- .prepare_gen_heatmap(
                data_table = data_table,
                df_heatmap = df_heatmap,
                df_heatmap_t = df_heatmap_t,
                significativityCriteria = significativityCriteria,
                pvalRange = pvalRange,
                fdrRange = fdrRange,
                numTopCNV = numTopCNV,
                numTopMET = numTopMET,
                scale = scale,
                numSamples = numSamples
            )
        }
        if (integrationSelect == "gene_cnv_res") {
            ans <- .prepare_cnv_heatmap(
                data_table = data_table,
                df_heatmap = df_heatmap,
                df_heatmap_t = df_heatmap_t,
                significativityCriteria = significativityCriteria,
                pvalRange = pvalRange,
                fdrRange = fdrRange,
                numTopCNVonly = numTopCNVonly,
                scale = scale,
                numSamples = numSamples
            )
        }
        if (integrationSelect == "gene_met_res") {
            ans <- .prepare_met_heatmap(
                data_table = data_table,
                df_heatmap = df_heatmap,
                df_heatmap_t = df_heatmap_t,
                significativityCriteria = significativityCriteria,
                pvalRange = pvalRange,
                fdrRange = fdrRange,
                numTopMETonly = numTopMETonly,
                scale = scale,
                numSamples = numSamples
            )
        }
        if (integrationSelect == "mirna_cnv_res") {
            ans <- .prepare_mirna_heatmap(
                data_table = data_table,
                df_heatmap = df_heatmap,
                df_heatmap_t = df_heatmap_t,
                significativityCriteria = significativityCriteria,
                pvalRange = pvalRange,
                fdrRange = fdrRange,
                numTopMiCNV = numTopMiCNV,
                scale = scale,
                numSamples = numSamples
            )
        }
        if (is.null(ans)) {
            return(NULL)
        }
        ht <- makeInteractiveComplexHeatmap(input, output, session, ans,
                                            ns("heatmap"))
    }) %>% bindEvent(
        input[[ns("IntegrationSelect")]],
        input[[ns("numTopGenesHeatmapCNV")]],
        input[[ns("numTopGenesHeatmapMET")]],
        input[[ns("numTopGenesHeatmapCNVonly")]],
        input[[ns("numTopGenesHeatmapMETonly")]],
        input[[ns("numTopGenesHeatmapmirna_cnv")]],
        input[[ns("ClassSelect")]],
        input[[ns("SignificativityCriteria")]],
        input[[ns("PvalRange")]],
        input[[ns("FdrRange")]],
        input[[ns("scaleHeatmap")]],
        input[[ns("numSamples")]]
    )
}

# Prepare Reactive Ridge Plot
#' @importFrom dplyr mutate_if %>%
#' @importFrom shiny reactive bindEvent
.prepare_reactive_ridge <- function(data_table,
                                    input,
                                    output,
                                    deg = FALSE,
                                    table = FALSE) {
    reactive({
        df <- data_table
        if (table) df <- mutate_if(df, is.numeric, ~ round(., 3))
        integrationSelect <- input$IntegrationSelect
        classSelect <- input$ClassSelect
        significativityCriteria <- input$SignificativityCriteria
        pvalRange <- input$PvalRange
        fdrRange <- input$FdrRange
        typeSelect <- input$genomicTypeSelect
        if (deg) df <- df[df$deg, ]
        if ("class" %in% colnames(data_table)) df <- df[df$class ==
                                                          classSelect, ]
        if (integrationSelect == "gene_genomic_res") {
            df <- df[df$cnv_met == typeSelect, ]
            df <- df[!is.na(df$cnv_met), ]
        }
        df <- df[df$omics == integrationSelect, ]
        if (table) {
            if (significativityCriteria == "pval") {
                df <- df[df$pval >= pvalRange[1] & df$pval <= pvalRange[2], ]
            } else {
                df <- df[df$fdr >= fdrRange[1] & df$fdr <= fdrRange[2], ]
            }
        } else {
            if (significativityCriteria == "pval") {
                df$significance <- ifelse(df$pval >= pvalRange[1] & df$pval <=
                                            pvalRange[2],
                    "Significant", "Not Significant"
                )
            } else {
                df$significance <- ifelse(df$fdr >= fdrRange[1] & df$fdr <=
                                            fdrRange[2],
                    "Significant", "Not Significant"
                )
            }
        }
        if (nrow(df) == 0) {
            return(NULL)
        }
        lower_quantile <- quantile(df$coef, 0.001)
        upper_quantile <- quantile(df$coef, 0.999)
        ans <- list(df = df, quantiles = c(lower_quantile, upper_quantile))
        if (table) {
            return(df)
        } else {
            return(ans)
        }
    }) %>% bindEvent(
        input$IntegrationSelect,
        input$ClassSelect,
        input$SignificativityCriteria,
        input$PvalRange,
        input$FdrRange,
        input$genomicTypeSelect
    )
}


# Prepare Reactive Histogram
#' @importFrom dplyr mutate_if %>%
#' @importFrom gtools mixedsort
#' @importFrom shiny reactive bindEvent
.prepare_reactive_histo <- function(data_table,
                                    input,
                                    output,
                                    deg = FALSE,
                                    table = FALSE) {
    reactive({
        if (table) data_table <- mutate_if(data_table, is.numeric, ~
                                             round(., 3))
        integrationSelect <- input$IntegrationSelect
        typeSelect <- input$genomicTypeSelect
        classSelect <- input$ClassSelect
        chrSelect <- input$ChrSelect
        significativityCriteria <- input$SignificativityCriteria
        pvalRange <- input$PvalRange
        fdrRange <- input$FdrRange
        if (!table) {
            data_table <- data_table[!is.na(data_table$chr_cov), ]
        }
        chr_order <- mixedsort(unique(data_table$chr_cov))
        chr_order <- chr_order[!is.na(chr_order)]
        data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
        data_table <- data_table[data_table$omics == integrationSelect, ]
        if (integrationSelect == "gene_genomic_res") {
            data_table <- data_table[data_table$cnv_met == typeSelect, ]
        }
        if (deg == TRUE) data_table <- data_table[data_table$deg, ]
        if (chrSelect != "All") {
            data_table <- data_table[data_table$chr_cov == chrSelect, ]
        }
        if ("class" %in% colnames(data_table)) {
            data_table <- data_table[data_table$class == classSelect, ]
        }
        if (table) {
            if (significativityCriteria == "pval") {
                data_table <- data_table[data_table$pval >= pvalRange[1] &
                    data_table$pval <= pvalRange[2], ]
            } else {
                data_table <- data_table[data_table$fdr >= fdrRange[1] &
                    data_table$fdr <= fdrRange[2], ]
            }
        } else {
            if (significativityCriteria == "pval") {
                data_table$significance <- ifelse(data_table$pval >=
                                                    pvalRange[1] &
                    data_table$pval <= pvalRange[2],
                "Significant",
                "Not Significant"
                )
            } else {
                data_table$significance <- ifelse(data_table$fdr >=
                                                    fdrRange[1] &
                    data_table$fdr <= fdrRange[2],
                "Significant",
                "Not Significant"
                )
            }
        }
        if (nrow(data_table) == 0) {
            return(NULL)
        } else {
            return(data_table)
        }
    }) %>% bindEvent(
        input$IntegrationSelect,
        input$genomicTypeSelect,
        input$ClassSelect,
        input$ChrSelect,
        input$SignificativityCriteria,
        input$PvalRange,
        input$FdrRange
    )
}

# Prepare Reactive Histogram for TF
#' @importFrom gtools mixedsort
#' @importFrom shiny reactive bindEvent
#' @importFrom dplyr %>%
.prepare_reactive_histo_tf <- function(data_table,
                                       input,
                                       output,
                                       deg = FALSE) {
    reactive({
        genes_count_df <- NULL
        integrationSelect <- input$IntegrationSelect
        typeSelect <- input$genomicTypeSelect
        classSelect <- input$ClassSelect
        chrSelect <- input$ChrSelect
        significativityCriteria <- input$SignificativityCriteria
        pvalRange <- input$PvalRange
        fdrRange <- input$FdrRange
        chr_order <- mixedsort(unique(data_table$chr_cov))
        chr_order <- chr_order[!is.na(chr_order)]
        data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)
        if (integrationSelect %in% c("tf_res", "mirna_target_res",
                                     "tf_mirna_res")) {
            data_table <- data_table[data_table$omics == integrationSelect, ]
            data_table <- data_table[
                ,
                colnames(data_table) %in% c(
                    "response", "cov",
                    "pval", "fdr", "chr_cov",
                    "deg", "class"
                )
            ]
            if ("class" %in% colnames(data_table)) {
                data_table <- data_table[data_table$class == classSelect, ]
            }
            if (nrow(data_table) == 0) {
                return(NULL)
            }
            if (deg) data_table <- data_table[data_table$deg, ]
            if (significativityCriteria == "pval") {
                data_table <- data_table[data_table$pval >= pvalRange[1] &
                    data_table$pval <= pvalRange[2], ]
            } else {
                data_table <- data_table[data_table$fdr >= fdrRange[1] &
                    data_table$fdr <= fdrRange[2], ]
            }
            genes_count <- table(data_table$cov, data_table$chr_cov)
            genes_count_df <- as.data.frame.table(genes_count)
            genes_count_df <- subset(genes_count_df, `Freq` != 0)
            colnames(genes_count_df) <- c("TF", "Chromosome", "Count")
            genes_count_df <- genes_count_df[order(-genes_count_df$Count), ]
        }
        return(genes_count_df)
    }) %>% bindEvent(
        input$IntegrationSelect,
        input$genomicTypeSelect,
        input$ClassSelect,
        input$ChrSelect,
        input$SignificativityCriteria,
        input$PvalRange,
        input$FdrRange
    )
}

# Prepare Reactive Table
#' @importFrom gtools mixedsort
#' @importFrom shiny reactive bindEvent
#' @importFrom dplyr %>% mutate_if
.prepare_reactive_table <- function(data_table,
                                    input,
                                    output) {
    reactive({
        data_table <- data_table[!is.na(data_table$chr_cov), ]
        data_table <- data_table[!is.na(data_table$response), ]
        data_table <- mutate_if(data_table, is.numeric, ~ round(., 3))
        chr_order <- mixedsort(unique(data_table$chr_cov))
        chr_order <- chr_order[!is.na(chr_order)]
        data_table$chr_cov <- factor(data_table$chr_cov, levels = chr_order)

        data_table <- data_table[data_table$omics == input$IntegrationSelect, ]
        if ("class" %in% colnames(data_table)) {
            data_table <- data_table[data_table$class == input$ClassSelect, ]
        }

        if (input$SignificativityCriteria == "pval") {
            data_table <- data_table[data_table$pval >= input$PvalRange[1] &
                data_table$pval <= input$PvalRange[2], ]
        } else {
            data_table <- data_table[data_table$fdr >= input$FdrRange[1] &
                data_table$fdr <= input$FdrRange[2], ]
        }
        if (input$degSelect == "Only DEGs") {
            data_table <- data_table[data_table$deg, ]
        }
        if (input$ChrSelect != "All") {
            data_table <- data_table[data_table$chr_cov == input$ChrSelect, ]
        }
        return(data_table)
    }) %>% bindEvent(
        input$IntegrationSelect,
        input$ClassSelect,
        input$ChrSelect,
        input$PvalRange,
        input$FdrRange,
        input$degSelect,
        input$SignificativityCriteria
    )
}


# Check Reactive Background Enrichment
#' @importFrom shiny invalidateLater reactive
.check_reactive_bg_enrich <- function(bg_enrich,
                                      input,
                                      output,
                                      session) {
    reactive({
        invalidateLater(millis = 1000, session = session)
        if (bg_enrich$is_alive()) {
            x <- "Enrichment running in background, this may take several
            minutes"
        } else {
            x <- "Enrichment completed"
        }
        return(x)
    })
}

# Reactive Gene Enrichment
#' @importFrom shiny reactive invalidateLater
.reactive_gen_enrich <- function(bg_enrich,
                                 input,
                                 output,
                                 session) {
    reactive({
        if (bg_enrich$is_alive()) {
            invalidateLater(millis = 1000, session = session)
        }
        class <- input$ClassSelect
        db <- input$DBSelectEnrich
        type <- input$genomicTypeSelect
        if (!bg_enrich$is_alive()) {
            data <- bg_enrich$get_result()
            if (sum(c("cnv", "met") %in% names(data)) == 2) {
                if (!class %in% names(data[[type]])) class <- 1
                ans <- data[[type]][[class]][[db]]
                if (is.null(ans)) {
                    return(NULL)
                }
                ans2 <- dot_plotly(ans)
                ans <- list(plot = ans2, table = ans@result)
            } else {
                if (!class %in% names(data)) class <- 1
                ans <- data[[class]][[db]]
                if (is.null(ans)) {
                    return(NULL)
                }
                ans2 <- dot_plotly(ans)
                ans <- list(plot = ans2, table = ans@result)
            }
            return(ans)
        }
    })
}

# Reactive TF Enrichment
#' @importFrom shiny reactive invalidateLater
.reactive_tf_enrich <- function(bg_enrich,
                                input,
                                output,
                                session) {
    reactive({
        if (bg_enrich$is_alive()) {
            invalidateLater(millis = 1000, session = session)
        }
        class <- input$ClassSelect
        db <- input$DBSelectEnrich
        if (!bg_enrich$is_alive()) {
            data <- bg_enrich$get_result()
            if (!class %in% names(data)) class <- 1
            tf <- names(data[[class]])
            ans <- lapply(tf, function(x) {
                ans <- data[[class]][[x]][[db]]
                if (!is.null(ans)) {
                    ans2 <- dot_plotly(ans, title = x)
                    ans <- list(plot = ans2, table = ans@result)
                    return(ans)
                } else {
                    return(NULL)
                }
            })
            names(ans) <- tf
            return(ans)
        }
    })
}

# Prepare Reactive Circos Plot
#' @importFrom shiny.gosling arrange_views
#' @importFrom shiny reactive bindEvent
#' @importFrom dplyr %>%
.prepare_reactive_circos <- function(data, input, output) {
    reactive({
        if ("class" %in% colnames(data)) data <- data[data$class ==
                                                        input$ClassSelect, ]
        gr <- .circos_preprocess(data = data)
        tracks <- .create_tracks(data = data, gr = gr)
        width <- 800
        height <- 800
        if (input$layout == "linear") {
            width <- 800
            height <- 100
        }
        composed_view <- .create_composed_view(tracks, height = height,
                                               width = width)
        if (input$circosType == "Gene") {
            ssel <- intersect(
                c("circos_genomic", "circos_cnv_gene", "circos_met_gene"),
                names(composed_view)
            )
        }
        if (input$circosType == "miRNA") {
            ssel <- intersect(c("circos_genomic_mirna"), names(composed_view))
        }
        composed_view <- composed_view[[ssel]]
        arranged_view <- arrange_views(
            views = composed_view,
            layout = input$layout,
            xDomain = list(
                chromosome = input$ChrSelect
            ),
            assembly = "hg38"
        )

        return(arranged_view)
    }) %>% bindEvent(
        input$ClassSelect,
        input$layout,
        input$circosType,
        input$ChrSelect
    )
}
