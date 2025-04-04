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
            if (nrow(network_data$nodes) == 0 | is.null(network_data)) {
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
                            deg = FALSE,
                            degs = NULL) {
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
        nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
        ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges, ]
        if (significativityCriteria == "pval") {
            ans$edges <- ans$edges[ans$edges$pval <= pval, ]
        } else {
            ans$edges <- ans$edges[ans$edges$fdr <= fdr, ]
        }
        if (deg) {
            if(is.null(degs)) return(NULL)
            tto <- unique(data_table$response[data_table[, class]])
            ans$edges <- ans$edges[ans$edges$to %in% tto, ]
            degs <- degs[[gsub("deg_", "", class)]]
            up <- degs[degs$logFC>0,]
            down <- degs[degs$logFC<0,]
            ans$nodes$color[ans$nodes$id%in%up$id] <- "#ED5564"
            ans$nodes$color[ans$nodes$id%in%down$id] <- "#4169E1"
            ans$legend_nodes$color <- "lightgrey"
            tmp <- data.frame(
              label = c("Up", "Down"),
              shape = "dot",
              color = c("#ED5564", "#4169E1"))
            ans$legend_nodes <- rbind(ans$legend_nodes, tmp)
        }
        ans$edges <- ans$edges[order(abs(ans$edges$coef),
                                     decreasing = TRUE), ]
        ans$edges <- ans$edges[seq_len(length.out = numInteractions), ]
        nodes_with_edges <- unique(c(ans$edges$from, ans$edges$to))
        ans$nodes <- ans$nodes[ans$nodes$id %in% nodes_with_edges, ]
        return(ans)
    }) %>% bindEvent(
        input$layout,
        input$numInteractions,
        input$SignificativityCriteria,
        input$PvalRange,
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
        if (!deg & "gene_genomic_res" %in% unique(data_table$omics)) {
            data_table <- data_table[data_table$omics == "gene_genomic_res", ]
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
        if (!deg & !"gene_genomic_res" %in% unique(data_table$omics)) {
            if ("gene_cnv_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_cnv_res", ]
            }
            if ("gene_met_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_met_res", ]
            }
        }
        if (deg & "gene_genomic_res" %in% unique(data_table$omics)) {
            if (!is.null(data_table[, classSelect])) {
                data_table <- data_table[data_table[,classSelect], ]
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
        if (deg & !"gene_genomic_res" %in% unique(data_table$omics)) {
          if (!is.null(data_table[, classSelect])) {
            data_table <- data_table[data_table[,classSelect], ]
            if ("gene_cnv_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_cnv_res", ]
            }
            if ("gene_met_res" %in% unique(data_table$omics)) {
                data_venn <- data_table[data_table$omics == "gene_met_res", ]
            }
          }else {
              (return(NULL))
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
        classSelect <- input$ClassSelect
        pvalRange <- input$PvalRange
        fdrRange <- input$FdrRange
        if (deg) data_table <- data_table[data_table[,classSelect], ]
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
        input$FdrRange,
        input$ClassSelect
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
        if (length(.get_deg_col(data_table))==0 & deg) {
            return(NULL)
        }
        if (!deg) {
            df_heatmap <- multiomics_integration[[
              integrationSelect]]$data$response_var
            data_table <- data_table[data_table$omics == integrationSelect, ]
        }
        if (deg) {
            df_heatmap <- multiomics_integration[[
              integrationSelect]]$data$response_var
            data_table <- data_table[data_table$omics == integrationSelect, ]
            data_table <- data_table[data_table[,classSelect], ]
            df_heatmap <- df_heatmap[, unique(data_table$response)]
        }
        if (is.null(df_heatmap)) {
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
        if (deg) df <- df[df[,classSelect], ]
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
        if (deg) data_table <- data_table[data_table[, classSelect], ]
        if (chrSelect != "All") {
            data_table <- data_table[data_table$chr_cov == chrSelect, ]
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
                    "pval", "fdr", "chr_cov", classSelect
                )
            ]
            if (deg) data_table <- data_table[data_table[, classSelect], ]
            if (nrow(data_table) == 0) {
                return(NULL)
            }
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
        if (input$SignificativityCriteria == "pval") {
            data_table <- data_table[data_table$pval >= input$PvalRange[1] &
                data_table$pval <= input$PvalRange[2], ]
        } else {
            data_table <- data_table[data_table$fdr >= input$FdrRange[1] &
                data_table$fdr <= input$FdrRange[2], ]
        }
        if (input$degSelect!="All_genes") {
            data_table <- data_table[data_table[, input$degSelect], ]
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
            x <- paste("Enrichment running in background, this may take",
                       "several minutes")
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
        db <- input$DBSelectEnrich
        type <- input$genomicTypeSelect
        if (!bg_enrich$is_alive()) {
            data <- bg_enrich$get_result()
            class <- 1
            if (sum(c("cnv", "met") %in% names(data)) == 2) {
                ans <- data[[type]][[class]][[db]]
                if (is.null(ans)) {
                    return(NULL)
                }
                ans2 <- dot_plotly(ans)
                ans <- list(plot = ans2, table = ans@result)
            } else {
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
        db <- input$DBSelectEnrich
        if (!bg_enrich$is_alive()) {
            data <- bg_enrich$get_result()
            class <- 1
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
        input$layout,
        input$circosType,
        input$ChrSelect
    )
}


# start Diffmet 
.diffMet_start <- function(input,
                           status,
                           session,
                           multiomics,
                           mmultiomics_file,
                           pid_file,
                           datafile,
                           logfile,
                           child_pids
                           ){

  observeEvent(input$start, {
    status$text <- "Calcolo in corso..."
    session$sendCustomMessage("disableButton", "DiffMet_deg-start")
    session$sendCustomMessage("enableButton", "DiffMet_deg-stop")
    saveRDS(multiomics, mmultiomics_file)
    file.create(pid_file)
    num_cores <- input$num_cores
    vfile <- c(
      'pid <- Sys.getpid()',
      paste0('fileConn <- file(\"', pid_file, '\")'),
      'writeLines(as.character(pid), fileConn)',
      'close(fileConn)',
      'library(BiocParallel)',
      'library(gINTomics)',
      paste0('mmultiomics <- readRDS(\"', mmultiomics_file, '\")'),
      paste0('bpparam <- MulticoreParam(', num_cores, ')'),
      paste0('result <- gINTomics:::.get_mo_filtered_genexp(mmultiomics, BPPARAM = bpparam, run_met=F)'),
      paste0('saveRDS(result, \"', datafile, '\")'),
      'on.exit({',
      '  rm(mmultiomics, result, bpparam)',
      '  gc()',
      '})',
      'rm(list = ls())',
      'gc()',
      paste0('fileConn <- file(\"', logfile, '\")'),
      'writeLines(\"COMPLETED\", fileConn)',
      'close(fileConn)'
    )
    scriptfile <- tempfile(fileext = ".R")
    writeLines(vfile, scriptfile)
    if (Sys.info()['sysname'] == "Windows") {
      system_cmd <- paste0("start /B ", Sys.getenv("R_HOME"), "/bin/Rscript ", scriptfile)
    } else {
      system_cmd <- paste0(Sys.getenv("R_HOME"), "/bin/Rscript ", scriptfile, " &")
    }
    system(system_cmd)
    Sys.sleep(2)
    main_pid <- as.numeric(readLines(pid_file))
    child_pids$main_pid <- main_pid
    print(paste("PID principale:", main_pid))
  })
  
  }
  
  

# stop Diffmet 
.diffMet_stop <- function(input,
                           status,
                           session,
                           mmultiomics_file,
                           pid_file,
                           datafile,
                           logfile,
                           child_pids
){
  
  observeEvent(input$stop, {
    status$text <- "Interruzione in corso..."
    
    if (!is.null(child_pids$main_pid)) {
      if (.Platform$OS.type == "unix") {
        child_pids_list <- system(paste("pgrep -P", child_pids$main_pid), intern = TRUE)
        child_pids$child_pids <- as.numeric(child_pids_list)
        
        if (length(child_pids$child_pids) > 0) {
          print(paste("PID figli trovati:", paste(child_pids$child_pids, collapse = ", ")))
        } else {
          print("Nessun PID figlio trovato.")
        }
        
        if (length(child_pids$child_pids) > 0) {
          child_pids2 <- isolate(child_pids$child_pids)  
          print(paste("Terminando processo figlio PID:", child_pids2[1]))
          system(paste("kill", child_pids2[1]))
        }
        print(paste("Terminando processo principale PID:", child_pids$main_pid))
        system(paste("kill", child_pids$main_pid))
        Sys.sleep(2)
        
        for (pid in c(child_pids$main_pid, child_pids$child_pids)) {
          if (system(paste("ps -p", pid), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
            print(paste("Forzando la chiusura del processo PID:", pid))
            system(paste("kill -9", pid))
          }
        }
      } else if (Sys.info()['sysname'] == "Windows") {
        print(paste("Terminando processo principale PID su Windows:", child_pids$main_pid))
        system(paste("taskkill /PID", child_pids$main_pid, "/F"))
        Sys.sleep(2)
      }
      session$sendCustomMessage("enableButton", "DiffMet_deg-start")
      session$sendCustomMessage("disableButton", "DiffMet_deg-stop")
      status$text <- "Calcolo interrotto!"
    } else {
      print("Errore: Nessun PID salvato, impossibile terminare il processo!")
    }
      gc()
      unlink(pid_file)
      unlink(logfile)
      unlink(datafile)
      unlink(mmultiomics_file)
      unlink(pid_file)
  })
  
}

# MethylMix observer
.reactive_methylmix <- function(input,
                               result_data,
                               multiomics,
                               session,
                               data,
                               cclass){
  
  observeEvent(list(input$class1, input$class2), {
    disable("genes")
    req(input$class1, input$class2)
    if(input$class1==input$class2) return(NULL)
    cnv_residuals <- result_data$res_expr_cnv
    if("gene_genomic_res"%in%names(multiomics)){
      tmp <- multiomics$gene_genomic_res$data$covariates
      meth_data <- tmp[, grep("_met$", colnames(tmp))]
      colnames(meth_data) <- gsub("_met$", "", colnames(meth_data))
      }
    if("gene_met_res"%in%names(multiomics)){
      meth_data <- multiomics$gene_met_res$data$covariates
      colnames(meth_data) <- gsub("_cov$", "", colnames(meth_data))
      }
    message(paste("Computing", as.character(input$class1)))
    MethylMixResults <- MethylMix(t(meth_data[names(cclass)[cclass==input$class1],]),
                                  as.matrix(cnv_residuals[,names(cclass)[cclass==input$class1]]),
                                  t(meth_data[names(cclass)[cclass==input$class2],]))
    message(paste("Done", as.character(input$class1)))
    tmp <- list(MethylMixResults=MethylMixResults,
                cnv_residuals=cnv_residuals,
                meth_data=meth_data,
                cclass=cclass
    )
    data(tmp)
    enable("genes")
    updateSelectInput(session, "genes", choices = names(MethylMixResults$Models), selected = names(MethylMixResults$Models)[1])
  })
}

# MethylMix plots observer
.reactive_methylmix_plots <- function(input,
                                      session,
                                      data,
                                      output){
  
  observeEvent(list(input$genes, data()), {
    req(input$genes, data())
    MethylMixResults <- data()$MethylMixResults
    meth_data <- data()$meth_data
    cnv_residuals <- data()$cnv_residuals
    cclass <- data()$cclass
    message(paste("Plotting", as.character(input$genes)))
    plots <- MethylMix_PlotModel(input$genes, MethylMixResults,
                                 t(meth_data[names(cclass)[cclass==input$class1],]),
                                 as.matrix(cnv_residuals[,names(cclass)[cclass==input$class1]]),
                                 t(meth_data[names(cclass)[cclass==input$class2],]))
    message(paste("Plotting2", as.character(input$genes)))
    output$plot1 <- renderPlot({plots$MixtureModelPlot})
    output$plot2 <- renderPlot({plots$CorrelationPlot})
    output$table <- renderDataTable({as.data.frame(MethylMixResults$MethylationDrivers)})
  })
  
  
}

  