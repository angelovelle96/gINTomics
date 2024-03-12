
#' Enrichment
#' @importFrom clusterProfiler enrichKEGG enrichGO
#' @importFrom ReactomePA enrichPathway
#' @importFrom stats setNames
#' @import org.Hs.eg.db org.Mm.eg.db

.def_enrich <- function(data,
                        species,
                        pvalueCutoff,
                        pAdjustMethod,
                        qvalueCutoff,
                        ont,
                        run_go=TRUE,
                        run_kegg=TRUE,
                        run_reactome=FALSE,
                        ...){

    kegg <- go <- reactome <- NULL
    orgdb <- setNames(c("org.Hs.eg.db", "org.Mm.eg.db"),
                      c("hsa", "mmu"))
    organism <- setNames(c("human", "mouse"),
                      c("hsa", "mmu"))
    orgdb <- orgdb[species]
    organism <- organism[species]
    ssel <- setNames(data$coef[data$pval<=0.01],
                     data$entrez_response[data$pval<=0.01])
    if(length(ssel)>500){
      tmp <- data[data$pval<=0.01,]
      tmp <- tmp[order(abs(tmp$coef), decreasing = TRUE),]
      tmp <- tmp$entrez_response[1:500]
      ssel <- ssel[names(ssel)%in%tmp]
    }
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
                           qvalueCutoff=qvalueCutoff,
                           ...
        )
    }
    if(run_go){
        go <- enrichGO(gene = as.character(names(ssel)),
                       universe = as.character(universe),
                       OrgDb = orgdb,
                       pvalueCutoff = pvalueCutoff,
                       pAdjustMethod=pAdjustMethod,
                       qvalueCutoff=qvalueCutoff,
                       ont = ont,
                       readable = TRUE,
                       ...
        )
    }
    if(run_reactome){
      reactome <- enrichPathway(gene = as.character(names(ssel)),
                               universe = as.character(universe),
                               organism = organism,
                               pvalueCutoff = pvalueCutoff,
                               pAdjustMethod=pAdjustMethod,
                               qvalueCutoff=qvalueCutoff,
                               readable = TRUE,
                               ...
      )
    }

    enrichment <- Filter(Negate(is.null), list(kegg=kegg,
                                               go=go,
                                               reactome=reactome))


    return(enrichment)
}


#' Running genomic enrichment analysis
#' @param model_results
#' @param species
#' @param pvalueCutoff
#' @param pAdjustMethod
#' @param qvalueCutoff description
#' @param ont
#' @param BPPARAM description
#' @param extracted_data description
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
                       qvalueCutoff=0.1,
                       ont = "all",
                       BPPARAM = SerialParam(),
                       extracted_data=NULL,
                       ...
){
  data <- extracted_data
  if(is.null(data)){
    if("gene_genomic_res"%in%names(model_results)){
    model_results <- model_results[["gene_genomic_res"]]
    }
    if("gene_cnv_res"%in%names(model_results)){
      model_results <- model_results[["gene_cnv_res"]]
    }
    if("gene_met_res"%in%names(model_results)){
      model_results <- model_results[["gene_met_res"]]
      }
    data <- extract_model_res(model_results = model_results)
  }
  data <- data[data$cov!="(Intercept)",]
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

  if(length(unique(data[[1]]$cnv_met[!is.na(data[[1]]$cnv_met)]))>0){
    tmp <- lapply(data, function(x) x[x$cnv_met=="cnv",])
    enrichment_cnv <- BiocParallel::bplapply(tmp, gINTomics:::.def_enrich,
                              species=species,
                              pvalueCutoff = pvalueCutoff,
                              pAdjustMethod=pAdjustMethod,
                              qvalueCutoff=qvalueCutoff,
                              ont = ont,
                              BPPARAM = BPPARAM,
                              ...)
    tmp <- lapply(data, function(x) x[x$cnv_met=="met",])
    enrichment_met <- BiocParallel::bplapply(tmp, gINTomics:::.def_enrich,
                               species=species,
                               pvalueCutoff = pvalueCutoff,
                               pAdjustMethod=pAdjustMethod,
                               qvalueCutoff=qvalueCutoff,
                               ont = ont,
                               BPPARAM = BPPARAM,
                               ...)

      enrichment <- list(cnv=enrichment_cnv,
                         met=enrichment_met)
  }else{

      enrichment <- BiocParallel::bplapply(data, gINTomics:::.def_enrich,
                             species=species,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             qvalueCutoff=qvalueCutoff,
                             ont = ont,
                             BPPARAM = BPPARAM,
                             ...)

  }


return(enrichment)

}


#' Running TF enrichment analysis
#' @param model_results
#' @param species
#' @param pvalueCutoff
#' @param qvalueCutoff description
#' @param pAdjustMethod
#' @param ont
#' @param BPPARAM
#' @param extracted_data description
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

run_tf_enrich <- function(model_results,
                          species="hsa",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1,
                          pAdjustMethod="BH",
                          ont = "all",
                          BPPARAM = BiocParallel::SerialParam(),
                          extracted_data=NULL,
                          ...
){

  data <- extracted_data
  if(is.null(data)){
    if("tf_res"%in%names(model_results)){
    model_results <- model_results[["tf_res"]]
    }
    data <- extract_model_res(model_results = model_results)
  }

  data <- data[data$cov!="(Intercept)",]
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
    tmp <- unique(x$cov)
    tmp2 <- lapply(tmp, function(y) {
      check <- x$pval[x$cov==y]
      check <- sum(check<=0.01)
      return(check)
      })
    names(tmp2) <- tmp
    tmp2 <- unlist(tmp2)
    tmp2 <- tmp2[tmp2>12]
    tmp2 <- sort(tmp2, decreasing = T)
    if(length(tmp2)>10) tmp2 <- tmp2[1:10]
    tmp <- lapply(names(tmp2), function(y) x[x$cov==y,])
    names(tmp) <- names(tmp2)
    enrichment <- BiocParallel::bplapply(tmp, gINTomics:::.def_enrich,
                    species=species,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    qvalueCutoff = qvalueCutoff,
                    ont = ont,
                    BPPARAM = BPPARAM)
    return(enrichment)
    })
  return(enrichment)
}


######################################
#' Making groups
#' @importFrom ComplexHeatmap pheatmap
.make_groups <- function(model_results){

  data <- model_results$data$covariates
  tmp <- extract_model_res(model_results)
  if("cnv_met"%in%colnames(tmp)){
    tmp <- tmp[tmp$cov!="(Intercept)",]
    tmp <- paste(tmp$cov[tmp$pval<=0.05],
                 tmp$cnv_met[tmp$pval<=0.05], sep = "_")
  }else{
    tmp <- tmp[tmp$cov!="(Intercept)",]
    tmp <- tmp$cov[tmp$pval<=0.05]
  }

 data1 <- data[, tmp]
 data1 <- scale(data1)
 tmp <- unique(gsub("_.*$", "", tmp))
 data2 <- log2(model_results$data$response_var[, tmp]+1)

 pheatmap(data1, scale = "none")
 pheatmap(data2, scale = "none")

gmax <- round(nrow(data1)/5)
mclust <- Mclust(data1, G = 1:gmax)
mclust2 <- Mclust(data2, G = 1:gmax)
}


######################################
#' Survival
#' @importFrom survival survfit
survival <- function(surv_time,
                     surv_event,
                     data){



  fit <- survfit(Surv(pfs, status) ~ treatment , data = tmp)
}
