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
#' @importFrom shiny observe bindEvent updateSliderInput
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
        mmax <- nrow(
          multiomics_integration[[integrationSelect]]$data$response_var)
        updateSliderInput(session = session, ns("numSamples"),
                          value = min(10, mmax),
                          min = 1,
                          max = mmax,
                          step = 5
        )
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


# Prepare Reactive Table
#' @importFrom gtools mixedsort
#' @importFrom shiny reactive bindEvent
#' @importFrom dplyr %>% mutate_if
#' @importFrom plyr rbind.fill
.prepare_reactive_degsTable <- function(multiomics,
                                    input,
                                    output) {
  
  reactive({
    contrast <- input$ClassSelect
    contrast <- gsub("deg_", "", contrast)
    ans <- lapply(multiomics, function(x) {
      if (!is.null(x$deg)) {
        return(x$deg)
      } else {
        return(NULL)
      }
    })
    ans <- Filter(Negate(is.null), ans)
    ans <- lapply(ans, function(x) cbind(genes=rownames(x[[contrast]]),
                                         x[[contrast]]))
    ans <- rbind.fill(ans)
    if(input$SignificativityCriteria=="pval"){
      ans <- ans[ans$PValue <= input$PvalRange, ]
    }else{
      ans <- ans[ans$FDR <= input$FdrRange, ]
    }
    ans[,2:ncol(ans)] <- round(ans[,2:ncol(ans)], digits = 5)
    return(ans)
  })%>% bindEvent(
    input$ClassSelect,
    input$PvalRange,
    input$FdrRange,
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
        if(is(bg_enrich, 'reactiveExpr')) bg_enrich <- bg_enrich()
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
        ont <- input$ont
        if(is.null(ont)) ont <- "none"
        if (!bg_enrich$is_alive()) {
            data <- bg_enrich$get_result()
            class <- 1
            if (sum(c("cnv", "met") %in% names(data)) == 2) {
                ans <- data[[type]][[class]][[db]]
                if (is.null(ans)) {
                    return(NULL)
                }
                if(ont%in%c("BP", "MF", "CC") & db=="go"){
                  ans@result <- ans@result[ans@result$ONTOLOGY==ont,]}
                ans2 <- dot_plotly(ans)
                ans <- list(plot = ans2, table = ans@result)
            } else {
                ans <- data[[class]][[db]]
                if (is.null(ans)) {
                    return(NULL)
                }
                if(ont%in%c("BP", "MF", "CC") & db=="go"){
                  ans@result <- ans@result[ans@result$ONTOLOGY==ont,]}
                ans2 <- dot_plotly(ans)
                ans <- list(plot = ans2, table = ans@result)
            }
            return(ans)
        }
    })
}

# Reactive TF Enrichment
#' @importFrom shiny reactive invalidateLater req
.reactive_tf_enrich <- function(bg_enrich,
                                input,
                                output,
                                session) {
    reactive({
      proc <- bg_enrich()
      req(proc)
      if (proc$is_alive()) {
        invalidateLater(1000, session)
        return(NULL)
      }
      if (proc$get_exit_status() != 0) {
        stop("process: ", proc$read_error())
      }
      data <- proc$get_result()
      ans2 <- dot_plotly(data[[1]], title = input$genes)
      ans <- list(
        plot = ans2,
        table = data[[1]]@result
      )
      return(ans)
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
#' @importFrom shiny observeEvent
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
    status$text <- "Filtering out CNV from expression data..."
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
      paste0('result <- gINTomics:::.get_mo_filtered_genexp(mmultiomics,',
             ' BPPARAM = bpparam, run_met=F)'),
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
    rscript_path <- file.path(Sys.getenv("R_HOME"), "bin", "Rscript")
    if (Sys.info()[["sysname"]] == "Windows") {
      # Windows: use "cmd.exe" with "start /B"
      system2("cmd.exe",
              args = c("/c", "start", "/B",
                       shQuote(rscript_path),
                       shQuote(scriptfile)),
              wait = FALSE)
    } else {
      # Unix-like (Linux/macOS)
      system2(rscript_path, args = shQuote(scriptfile), wait = FALSE)
    }
    Sys.sleep(2)
    main_pid <- as.numeric(readLines(pid_file))
    child_pids$main_pid <- main_pid
    print(paste("Main PID:", main_pid))
  })
  
  }
  
  

# stop Diffmet 
#' @importFrom shiny observeEvent
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
    status$text <- "Stopping the process..."
    
    if (!is.null(child_pids$main_pid)) {
      if (.Platform$OS.type == "unix") {
        child_pids_list <- system2("pgrep",
                                   args = c("-P",
                                            as.character(child_pids$main_pid)),
                                   stdout = TRUE)
        
        child_pids$child_pids <- as.numeric(child_pids_list)
        
        if (length(child_pids$child_pids) > 0) {
          print(paste("PIDs:", paste(child_pids$child_pids, collapse = ", ")))
        } else {
          print("No PIDs found.")
        }
        
        if (length(child_pids$child_pids) > 0) {
          child_pids2 <- isolate(child_pids$child_pids)  
          print(paste("Stopping PID:", child_pids2[1]))
          system2("kill", args = as.character(child_pids2[1]))
        }
        print(paste("Stopping main PID:", child_pids$main_pid))
        system2("kill", args = as.character(child_pids$main_pid))
        Sys.sleep(2)
        
        for (pid in c(child_pids$main_pid, child_pids$child_pids)) {
          check <- system2("ps",args = c("-p", pid),stdout = NULL,stderr = NULL)
          if (check == 0) {
            print(paste("Forcing PID to stop:", pid))
            system2("kill", args = c("-9", as.character(pid)))
          }
        }
      } else if (Sys.info()['sysname'] == "Windows") {
        print(paste("Stopping main PID:", child_pids$main_pid))
        system2("taskkill", args = c("/PID",
                                     as.character(child_pids$main_pid), "/F"))
        Sys.sleep(2)
      }
      session$sendCustomMessage("enableButton", "DiffMet_deg-start")
      session$sendCustomMessage("disableButton", "DiffMet_deg-stop")
      status$text <- "Stopped!"
    } else {
      print("Errore: No PID saved, impossible to stop the process!")
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
#' @importFrom shiny observeEvent req updateSelectizeInput
.reactive_methylmix <- function(input,
                               result_data,
                               multiomics,
                               session,
                               ddata,
                               cclass,
                               status){
  
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
    MethylMixResults <- MethylMix(
      t(meth_data[names(cclass)[cclass==input$class1],]),
      as.matrix(cnv_residuals[,names(cclass)[cclass==input$class1]]),
      t(meth_data[names(cclass)[cclass==input$class2],]))
    tmp <- list(MethylMixResults=MethylMixResults,
                cnv_residuals=cnv_residuals,
                meth_data=meth_data,
                cclass=cclass
    )
    ddata(tmp)
    enable("genes")
    status$text <- "Analysis Completed !"
    updateSelectizeInput(session, "genes",
                         choices = names(MethylMixResults$Models),
                         selected = names(MethylMixResults$Models)[1],
                         server = TRUE)
  })
}

# MethylMix plots observer
#' @importFrom shiny observeEvent req
.reactive_methylmix_plots <- function(input,
                                      session,
                                      ddata,
                                      data_table,
                                      output){
  
  observeEvent(list(input$genes, ddata()), {
    req(input$genes, ddata())
    MethylMixResults <- ddata()$MethylMixResults
    meth_data <- ddata()$meth_data
    cnv_residuals <- ddata()$cnv_residuals
    cclass <- ddata()$cclass
    plots <- MethylMix_PlotModel(
      input$genes, MethylMixResults,
      t(meth_data[names(cclass)[cclass==input$class1],]),
      as.matrix(cnv_residuals[,names(cclass)[cclass==input$class1]]),
      t(meth_data[names(cclass)[cclass==input$class2],]))
    ttable <- MethylMixResults$MethylationDrivers
    data_table <- data_table[data_table$omics%in%c("gene_genomic_res",
                                                   "gene_met_res"),]
    if("gene_genomic_res"%in%data_table$omics){
      data_table <- data_table[data_table$cnv_met=="met",]
    }
    ttable <- data_table[data_table$response%in%ttable,]
    output$plot1 <- renderPlot({plots$MixtureModelPlot})
    output$plot2 <- renderPlot({plots$CorrelationPlot})
    output$table <- renderDataTable({
      as.data.frame(ttable)
      })
  })
  
  
}


# CNV ttest
#' @importFrom plyr rbind.fill
#' @importFrom shiny eventReactive req
#' @importFrom stats p.adjust t.test
.reactive_cnv_test <- function(input,
                               session,
                               cnv,
                               output,
                               cclass){
  
  eventReactive(input$start, {
    req(cclass, input$class1, input$class2)
    disable("genes")
    ans <- lapply(colnames(cnv), function(x){
      t.test(cnv[names(cclass)[cclass==input$class1], x],
             cnv[names(cclass)[cclass==input$class2], x])
    })
    names(ans) <- colnames(cnv)
    ans <- lapply(ans, function(x){
      as.data.frame(t(c(x$statistic, pvalue=x$p.value,
                        conf=x$conf.int, x$estimate)))
    })
    tmp <- rbind.fill(ans)
    rownames(tmp) <- names(ans)
    ans <- tmp
    ans$FDR <- p.adjust(ans$pvalue, method = "fdr")
    return(ans)
  })
  
}

# CNV ttest table
#' @importFrom shiny observeEvent req updateSelectizeInput
.render_cnv_test_table <- function(input,
                                   session,
                                   ddata,
                                   status,
                                   output){
  observeEvent(list(input$SignificativityCriteria,
                    input$PvalRange,
                    input$FdrRange,
                    ddata()), {
    req(ddata())
    status$text <- paste("Analysis Completed !")
    ssign <- input$SignificativityCriteria
    pvalRange <- input$PvalRange
    fdrRange <- input$FdrRange
    ans <- ddata()
    if(ssign=="pval") ans <- ans[ans$pvalue>=pvalRange[1] &
                                   ans$pvalue<=pvalRange[2],]
    if(ssign=="FDR") ans <- ans[ans$FDR>=fdrRange[1] &
                                  ans$FDR<=fdrRange[2],]
    updateSelectizeInput(session, "genes",
                         choices = rownames(ans),
                         selected = rownames(ans)[1],
                         server = TRUE)
    enable("genes")
    output$table <- renderDataTable({round(ans, digits = 5)})
  })
}


# CNV ttest plots
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme_minimal labs
#' geom_boxplot scale_fill_brewer theme
#' @importFrom shiny req reactive
.reactive_cnv_test_plots <- function(input,
                                     session,
                                     cnv,
                                     eexpr,
                                     output,
                                     cclass=NULL){
  reactive({
    req(input$genes, cnv)
    tmp <- data.frame(x=cnv[names(cclass), input$genes],
                      y=eexpr[names(cclass), input$genes],
                      class=cclass)
    ans <- ggplot(tmp, aes(x = x, y = log2(y+1), colour = class)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, color = "#404080",
                  linetype = "dashed") +
      theme_minimal(base_size = 14) +
      labs(
        title = paste(input$genes, "scatterplot"),
        subtitle = "Correlation between gene expression and CNV",
        x = "CNV",
        y = "Log2 Gene expression"
      )
    
    ans2 <- ggplot(tmp, aes(x = class, y = x, fill = class)) +
        geom_boxplot(alpha = 0.7, outlier.color = "red", outlier.shape = 16) +
        theme_minimal(base_size = 14) +
        labs(
          title = paste(input$genes, "CNV values"),
          x = "Class",
          y = "CNV"
        ) +
        theme(legend.position = "none")
    return(list(ans, ans2))
    
  })%>% bindEvent(input$genes)
}


# Transcriptional integration expression plots
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme_minimal labs
#' geom_boxplot scale_fill_brewer theme
#' @importFrom shiny req reactive
.reactive_transcExpr_plots <- function(input,
                                       session,
                                       multiomics,
                                       output){
  reactive({
    
    req(input$gene1, input$gene2)
    data <- multiomics[[input$IntegrationSelect]]
    cclass <- attr(data, "Class")
    cov <- data$data$covariates
    res <- data$data$response
    if(is.null(cclass)) cclass <- setNames(rep("", nrow(res)), rownames(res))
    tmp <- data.frame(x=cov[names(cclass), input$gene2],
                      y=res[names(cclass), input$gene1],
                      class=cclass)
    ans <- ggplot(tmp, aes(x = log2(x+1), y = y, colour = class)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, color = "#404080",
                  linetype = "dashed") +
      theme_minimal(base_size = 14) +
      labs(
        title = paste(input$gene1, gsub("_cov", "", input$gene2),
                      "scatterplot"),
        subtitle = "Correlation between regulator-target expression",
        x = paste(gsub("_cov", "", input$gene2),"Log2 expression"),
        y = paste(input$gene1,"genomic-filtered expression")
      )

    ans2 <- ggplot(tmp, aes(x = class, y = log2(x+1), fill = class)) +
      geom_boxplot(alpha = 0.7, outlier.color = "red", outlier.shape = 16) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste(gsub("_cov", "", input$gene2), "Log2 Expression"),
        x = "Class",
        y = "Log2 Expression"
      ) +
      theme(legend.position = "none")
    return(list(ans, ans2))
    
  })%>% bindEvent(input$gene1,
                  input$gene2,
                  input$IntegrationSelect)
}


# Running  reactive TF enrichment analysis
#' @importFrom shiny req reactive
.run_reactive_tf_enrich <- function(input,
                                    output,
                                    session,
                                    extracted_data,
                                    qvalueCutoff = 1,
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "none",
                                    species="hsa") {
  reactive({
    req(input$genes, input$DBSelectEnrich, input$ont)
    data <- extracted_data
    data <- data[data$cov != "(Intercept)", ]
    tmp <- data[data$cov == input$genes, ]
    run_go <- run_kegg <- run_reactome <- FALSE
    if (input$DBSelectEnrich == "KEGG") run_kegg <- TRUE
    if (input$DBSelectEnrich == "Reactome") run_reactome <- TRUE
    if (input$DBSelectEnrich == "GO") run_go <- TRUE
    bg_enr <- isolate({
      .run_bg(
        FFUN = .def_enrich,
        args = list(
          data = tmp,
          species = species,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          ont = input$ont,
          run_kegg = run_kegg,
          run_go = run_go,
          run_reactome = run_reactome
        )
      )
    })
    
    
    bg_process <<- bg_enr
    
    return(bg_enr)
  })
}



  