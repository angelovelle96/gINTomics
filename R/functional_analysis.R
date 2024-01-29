

#' @importFrom clusterProfiler enrichKEGG enrichGO
#' @importFrom ReactomePA enrichPathway
#' @import org.Hs.eg.db org.Mm.eg.db

.def_enrich <- function(data,
                        species,
                        pvalueCutoff,
                        pAdjustMethod,
                        ont,
                        run_go=T,
                        run_kegg=T,
                        run_reactome=T,
                        ...){

    kegg <- go <- reactome <- NULL
    orgdb <- setNames(c("org.Hs.eg.db", "org.Mm.eg.db"),
                      c("hsa", "mmu"))
    organism <- setNames(c("human", "mouse"),
                      c("hsa", "mmu"))
    orgdb <- orgdb[species]
    organism <- organism[species]
    ssel <- setNames(data$coef[data$pval<=0.05],
                     data$entrez_response[data$pval<=0.05])
    ssel <- ssel[!is.na(names(ssel))]
    ssel <- ssel[!is.na(ssel)]
    ssel <- ssel[!duplicated(names(ssel))]
    universe <- data$entrez_response
    universe <- universe[!is.na(universe)]
    universe <- universe[!duplicated(universe)]
    if(run_kegg){
        kegg <- enrichKEGG(gene = as.character(names(ssel)),
                           universe = as.character(universe),
                           organism = species,
                           pvalueCutoff = pvalueCutoff,
                           pAdjustMethod=pAdjustMethod,
                           ...
        )
    }
    if(run_go){
        go <- enrichGO(gene = as.character(names(ssel)),
                       universe = as.character(universe),
                       OrgDb = orgdb,
                       pvalueCutoff = pvalueCutoff,
                       pAdjustMethod=pAdjustMethod,
                       ont = ont,
                       readable = T,
                       ...
        )
    }
    if(run_reactome){
      reactome <- enrichPathway(gene = as.character(names(ssel)),
                               universe = as.character(universe),
                               organism = organism,
                               pvalueCutoff = pvalueCutoff,
                               pAdjustMethod=pAdjustMethod,
                               readable = T,
                               ...
      )
    }

    enrichment <- Filter(Negate(is.null), list(kegg=kegg,
                                               go=go,
                                               reactome=reactome))


    return(enrichment)
}


#' Title
#'
#' @param model_results
#' @param species
#' @param pvalueCutoff
#' @param pAdjustMethod
#' @param ont
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

run_genomic_enrich <- function(model_results,
                       species="hsa",
                       pvalueCutoff = 0.1,
                       pAdjustMethod="BH",
                       ont = "all",
                       ...
){

  if("gene_genomic_res"%in%names(model_results)){
    model_results <- model_results[["gene_genomic_res"]]
  }

  data <- extract_model_res(model_results = model_results)
  if("class"%in%colnames(data)){
    tmp <- lapply(unique(data$class), function(x) data[data$class==x,])
    names(tmp) <- unique(data$class)
    data <- tmp
  }

  if(is.data.frame(data)){
    tmp <- list()
    tmp[[1]] <- data
    data <- tmp
  }

  if("cnv_met"%in%colnames(data[[1]])){
        enrichment_cnv <- lapply(data, function(x)
        .def_enrich(data = x[x$cnv_met=="cnv",],
                    species=species,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    ...))
      enrichment_met <- lapply(data, function(x)
        .def_enrich(data = x[x$cnv_met=="met",],
                    species=species,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    ...))

      enrichment <- list(cnv=enrichment_cnv,
                         met=enrichment_met)
  }else{

      enrichment <- lapply(data, function(x)
        .def_enrich(data = x,
                    species=species,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    ...))

  }


return(enrichment)

}


#' Title
#'
#' @param model_results
#' @param species
#' @param pvalueCutoff
#' @param pAdjustMethod
#' @param ont
#' @param BPPARAM
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

run_tf_enrich <- function(model_results,
                          species="hsa",
                          pvalueCutoff = 0.1,
                          pAdjustMethod="BH",
                          ont = "all",
                          BPPARAM = BiocParallel::SerialParam(),
                          ...
){

  if("tf_res"%in%names(model_results)){
    model_results <- model_results[["tf_res"]]
  }

  data <- extract_model_res(model_results = model_results)
  if("class"%in%colnames(data)){
    tmp <- lapply(unique(data$class), function(x) data[data$class==x,])
    names(tmp) <- unique(data$class)
    data <- tmp
  }

  if(is.data.frame(data)){
    tmp <- list()
    tmp[[1]] <- data
    data <- tmp
  }


  enrichment <- lapply(data, function(x){
    x <- x[x$cov!="(Intercept)",]
    tmp <- unique(x$cov)
    tmp2 <- lapply(tmp, function(y) {
      check <- x$pval[x$cov==y]
      check <- sum(check<=0.05)
      return(check)
      })
    names(tmp2) <- tmp
    tmp2 <- unlist(tmp2)
    tmp2 <- tmp2[tmp2>12]
    tmp2 <- sort(tmp2, decreasing = T)
    if(length(tmp2)>10) tmp2 <- tmp2[1:10]
    tmp <- lapply(names(tmp2), function(y) x[x$cov==y,])
    names(tmp) <- names(tmp2)
    ans <- bplapply(tmp, .def_enrich,
                    species=species,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    run_kegg = T,
                    BPPARAM = BPPARAM)
    return(ans)
    })
  return(enrichment)
}






