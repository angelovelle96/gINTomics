# Enrichment
#' @importFrom clusterProfiler enrichKEGG enrichGO
#' @importFrom ReactomePA enrichPathway
#' @importFrom stats setNames
#' @import org.Hs.eg.db org.Mm.eg.db
.def_enrich <- function(data, species, pvalueCutoff,
    pAdjustMethod, qvalueCutoff, ont, run_go = TRUE,
    run_kegg = TRUE, run_reactome = FALSE, ...) {
    kegg <- go <- reactome <- NULL
    orgdb <- setNames(c("org.Hs.eg.db", "org.Mm.eg.db"),
        c("hsa", "mmu"))
    organism <- setNames(c("human", "mouse"), c("hsa",
        "mmu"))
    orgdb <- orgdb[species]
    organism <- organism[species]
    ssel <- setNames(data$coef[data$pval <= 0.01],
        data$entrez_response[data$pval <= 0.01])
    if (length(ssel) > 500) {
        tmp <- data[data$pval <= 0.01, ]
        tmp <- tmp[order(abs(tmp$coef), decreasing = TRUE),
            ]
        tmp <- tmp$entrez_response[seq_len(500)]
        ssel <- ssel[names(ssel) %in% tmp]
    }
    ssel <- ssel[!is.na(names(ssel))]
    ssel <- ssel[!is.na(ssel)]
    ssel <- ssel[!duplicated(names(ssel))]
    universe <- data$entrez_response
    universe <- universe[!is.na(universe)]
    universe <- universe[!duplicated(universe)]
    if (run_kegg) {
        kegg <- enrichKEGG(gene = as.character(names(ssel)),
            universe = as.character(universe), organism = species,
            pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
            qvalueCutoff = qvalueCutoff, ...)
    }
    if (run_go) {
        go <- enrichGO(gene = as.character(names(ssel)),
            universe = as.character(universe), OrgDb = orgdb,
            pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
            qvalueCutoff = qvalueCutoff, ont = ont,
            readable = TRUE, ...)
    }
    if (run_reactome) {
        reactome <- enrichPathway(gene = as.character(names(ssel)),
            universe = as.character(universe), organism = organism,
            pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
            qvalueCutoff = qvalueCutoff, readable = TRUE,
            ...)
    }
    enrichment <- Filter(Negate(is.null), list(kegg = kegg,
        go = go, reactome = reactome))
    return(enrichment)
}
#' Running genomic enrichment analysis
#' @param model_results Model integration results, typically a list containing
#' different types of genomic results
#' @param species Species to select for the enrichment analysis. Default is
#' 'hsa' (Homo sapiens).
#' @param pvalueCutoff P-value cutoff for significant enrichment. Default is
#' 0.1.
#' @param pAdjustMethod Method for adjusting p-values. Default is 'BH'
#' (Benjamini & Hochberg).
#' @param qvalueCutoff Q-value cutoff for significant enrichment. Default is
#' 0.1.
#' @param ont Ontology to use for the enrichment analysis. Default is 'all'.
#' @param BPPARAM A BiocParallelParam object specifying parallelization
#' options. Default is BiocParallel::SerialParam().
#' @param extracted_data Pre-extracted data for enrichment analysis. If NULL,
#' function will extract relevant data from model_results.
#' @param ... Additional arguments to be passed to the internal enrichment
#' function.
#' @return A list containing enrichment results. If CNV and methylation data
#' are available, it returns a nested list with results for each data type.
#' @examples
#' # Example usage:
#' data(ov_test_tcga_omics)
#' #multiomics_integration <- run_multiomics(mmultiassay_ov)
#' #gen_enr <- run_genomic_enrich(multiomics_integration, qvalueCutoff = 1,
#'  pvalueCutoff = 0.05, pAdjustMethod = 'none')
#' @export
run_genomic_enrich <- function(model_results, species = "hsa",
    pvalueCutoff = 0.1, pAdjustMethod = "BH", qvalueCutoff = 0.1,
    ont = "all", BPPARAM = BiocParallel::SerialParam(), extracted_data = NULL,
    ...) {
    data <- extracted_data
    if (is.null(data)) {
        if ("gene_genomic_res" %in% names(model_results)) {
            model_results <- model_results[["gene_genomic_res"]]
        }
        if ("gene_cnv_res" %in% names(model_results)) {
            model_results <- model_results[["gene_cnv_res"]]
        }
        if ("gene_met_res" %in% names(model_results)) {
            model_results <- model_results[["gene_met_res"]]
        }
        data <- extract_model_res(model_results = model_results)
    }
    data <- data[data$cov != "(Intercept)", ]
    if ("class" %in% colnames(data)) {
        tmp <- lapply(unique(data$class), function(x) data[data$class ==
            x, ])
        names(tmp) <- unique(data$class)
        data <- tmp
    }
    if (is.data.frame(data)) {
        tmp <- list()
        tmp[[1]] <- data
        data <- tmp
    }
    if (length(unique(data[[1]]$cnv_met[!is.na(data[[1]]$cnv_met)])) >
        0) {
        tmp <- lapply(data, function(x) x[x$cnv_met == "cnv",
            ])
        enrichment_cnv <- BiocParallel::bplapply(tmp, gINTomics:::.def_enrich,
            species = species, pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
            ont = ont, BPPARAM = BPPARAM, ...)
        tmp <- lapply(data, function(x) x[x$cnv_met == "met",
            ])
        enrichment_met <- BiocParallel::bplapply(tmp, gINTomics:::.def_enrich,
            species = species, pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
            ont = ont, BPPARAM = BPPARAM, ...)
        enrichment <- list(cnv = enrichment_cnv, met = enrichment_met)
    } else {
        enrichment <- BiocParallel::bplapply(data, gINTomics:::.def_enrich,
            species = species, pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
            ont = ont, BPPARAM = BPPARAM, ...)
    }
    return(enrichment)
}
#' Running TF enrichment analysis
#' @param model_results Model integration results, typically a list containing
#' TF data.
#' @param species Species to select for the enrichment analysis. Default is
#' 'hsa' (Homo sapiens).
#' @param pvalueCutoff P-value cutoff for significant enrichment. Default is
#' 0.1.
#' @param qvalueCutoff Q-value cutoff for significant enrichment. Default is
#' 0.1.
#' @param pAdjustMethod Method for adjusting p-values. Default is 'BH'
#' (Benjamini & Hochberg).
#' @param ont Ontology to use for the enrichment analysis. Default is 'all'.
#' @param BPPARAM A BiocParallelParam object specifying parallelization
#' options. Default is BiocParallel::SerialParam().
#' @param extracted_data Pre-extracted data for enrichment analysis. If NULL,
#' function will extract relevant data from model_results.
#' @param ... Additional arguments to be passed to the internal enrichment
#' function.
#' @return A list containing TF enrichment results.
#' @examples
#' # Example usage:
#' data(ov_test_tcga_omics)
#' #multiomics_integration <- run_multiomics(mmultiassay_ov)
#' #run_tf_enrich(multiomics_integration, qvalueCutoff = 1, pvalueCutoff = 0.05,
#'  pAdjustMethod = 'none')
#' @export
run_tf_enrich <- function(model_results, species = "hsa", pvalueCutoff = 0.1,
    qvalueCutoff = 0.1, pAdjustMethod = "BH", ont = "all",
    BPPARAM = BiocParallel::SerialParam(), extracted_data = NULL,
    ...) {
    data <- extracted_data
    if (is.null(data)) {
        if ("tf_res" %in% names(model_results)) {
            model_results <- model_results[["tf_res"]]
        }
        data <- extract_model_res(model_results = model_results)
    }
    data <- data[data$cov != "(Intercept)", ]
    if ("class" %in% colnames(data)) {
        tmp <- lapply(unique(data$class), function(x) data[data$class ==
            x, ])
        names(tmp) <- unique(data$class)
        data <- tmp
    }
    if (is.data.frame(data)) {
        tmp <- list()
        tmp[[1]] <- data
        data <- tmp
    }
    enrichment <- lapply(data, function(x) {
        tmp <- unique(x$cov)
        tmp2 <- lapply(tmp, function(y) {
            check <- x$pval[x$cov == y]
            check <- sum(check <= 0.01)
            return(check)
        })
        names(tmp2) <- tmp
        tmp2 <- unlist(tmp2)
        tmp2 <- tmp2[tmp2 > 12]
        tmp2 <- sort(tmp2, decreasing = TRUE)
        if (length(tmp2) > 10)
            tmp2 <- tmp2[seq_len(10)]
        tmp <- lapply(names(tmp2), function(y) x[x$cov == y,
            ])
        names(tmp) <- names(tmp2)
        enrichment <- BiocParallel::bplapply(tmp, gINTomics:::.def_enrich,
            species = species, pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
            ont = ont, BPPARAM = BPPARAM)
        return(enrichment)
    })
    return(enrichment)
}
