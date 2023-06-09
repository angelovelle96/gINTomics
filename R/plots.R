
setGeneric("extract_model_res", function(model_results, ...) standardGeneric("extract_model_res"))

setMethod("extract_model_res", "list",
          function(model_results,
                   outliers=F,
                   species="hsa"){

  data <- cbind(response=rownames(model_results$coef_data),
                model_results$coef_data)
  data <- reshape2::melt(data, id.vars="response", variable.name = "cov")
  rownames(data) <- paste0(data$response, "_", data$cov)
  data <- data[!is.na(data$value),]
  pval <- cbind(response=rownames(model_results$pval_data),
                model_results$pval_data)
  pval <- reshape2::melt(pval, id.vars="response", variable.name = "cov")
  rownames(pval) <- paste0(pval$response, "_", pval$cov)
  pval <- pval[rownames(data),]
  tmp <- rep("not_significant", nrow(pval))
  tmp[pval$value<=0.05] <- "significant"
  names(tmp) <- rownames(pval)
  data$pval <- pval[rownames(data), "value"]
  data$significativity <- tmp[rownames(data)]
  data$sign <- rep("negative", nrow(data))
  data$sign[data$value>0]="positive"
  genes_info <- integrazione::genes_annotation[[species]]
  tmp <- grep("^CHR_", genes_info$chromosome_name)
  genes_info <- genes_info[-tmp,]
  tmp <- grep("^K", genes_info$chromosome_name)
  genes_info <- genes_info[-tmp,]
  tmp <- grep("^G", genes_info$chromosome_name)
  genes_info <- genes_info[-tmp,]
  data$cov <- as.character(data$cov)
  for(i in 1:nrow(data)){
    data$cov[i] <- gsub("cov", data$response[i], as.character(data$cov[i]))
  }
  tmp <- c(intersect(data$cov, genes_info$hgnc_symbol))
  tmp <- rbind(genes_info[genes_info$hgnc_symbol%in%tmp,])
  tmp <- tmp[!duplicated(tmp$hgnc_symbol),]
  rownames(tmp) <- tmp$hgnc_symbol
  tmp2 <- intersect(data$cov, genes_info$ensembl_gene_id)
  tmp2 <- genes_info[genes_info$ensembl_gene_id%in%tmp2,]
  tmp2 <- tmp2[!duplicated(tmp2$ensembl_gene_id),]
  rownames(tmp2) <- tmp2$hgnc_symbol
  genes_info <- rbind(tmp, tmp2)
  data$chr_cov <- genes_info[as.character(data$cov), "chromosome_name"]
  data$cytoband_cov <- genes_info[as.character(data$cov), "band"]

  mmin <- quantile(data$value, 0.25) - 1.5*IQR(data$value)
  mmax <- quantile(data$value, 0.75) + 1.5*IQR(data$value)

  if(outliers==F){
    data <- data[data$value>mmin,]
    data <- data[data$value<mmax,]
  }

  return(data)

      }
)


#' @import ggplot2 ggridges

ridgeline_plot <- function(data,
                           outliers=F){

    data <- extract_model_res(data)
    ggplot(data, aes(x = value, y = significativity, fill=significativity))+
      geom_density_ridges() +
      theme_ridges()

}


#' @import ggplot2 ggridges
#' @importFrom gtools mixedsort

chr_distribution_plot <- function(data,
                           outliers=F,
                           show_sign=F){

  data <- extract_model_res(data)
  data <- data[!is.na(data$chr_cov),]
  if(show_sign){
    ggplot(data, aes(x = factor(chr_cov, level=mixedsort(unique(data$chr_cov))),
                     fill=sign))+
      geom_bar(position="dodge", stat="count")+
      xlab("Chromosome")

  }else{
  ggplot(data, aes(x = factor(chr_cov, level=mixedsort(unique(data$chr_cov))),
                   fill=significativity))+
    geom_bar(position="dodge", stat="count")+
      xlab("Chromosome")
        }
}

#' @import ComplexHeatmap

heatmap_sign <- function(data,
                         outliers=T,
                         number=50){

  data <- extract_model_res(tf_mirna_res, outliers=outliers)
  data <- data[data$cov!="(Intercept)",]
  data <- data[data$significativity=="significant",]
  tmp <- sort(table(data$response), decreasing = T)
  if(length(tmp)>number) tmp <- tmp[1:number]
  pval <- lapply(names(tmp), function(x) data[data$response%in%x, "pval"])
  names(pval) <- names(tmp)
  pval <- sapply(pval, min)

}







