
#' @import ggplot2 ggridges
#' @export
#'

ridgeline_plot <- function(data,
                           outliers=F){

    data <- extract_model_res(data)
    ggplot(data, aes(x = value, y = significativity, fill=significativity))+
      geom_density_ridges() +
      theme_ridges()

}


#' @import ggplot2 ggridges
#' @importFrom gtools mixedsort
#' @export

chr_distribution_plot <- function(data,
                           outliers=F,
                           show_sign=F){

  data <- extract_model_res(data)
  data <- data[!is.na(data$chr_cov),]
  if(show_sign){
    data <- data[data$significativity=="significant",]
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
#' @export

heatmap_sign <- function(data,
                         outliers=T,
                         number=50){

  data <- extract_model_res(mmultiomics, outliers=outliers)
  data <- data[data$cov!="(Intercept)",]
  data <- data[data$pval<=0.1,]
  if(!"omics"%in%colnames(data)) data$omics <-rep("omic", nrow(data))
  tmp <- unique(data$omics)
  data <- lapply(tmp, function(x) data[data$omics==x,])
  names(data) <- tmp
  tmp <- lapply(data, function(x) sort(table(x$response), decreasing = T))
  tmp2 <- lapply(data, function(x) sort(setNames(x$pval, x$response)))
  tmp2 <- lapply(tmp2, function(x) x[!duplicated(names(x))])
  tmp2 <- lapply(tmp2, function(x) 1-x)
  tmp3 <- lapply(names(tmp), function(x) tmp[[x]]+tmp2[[x]][names(tmp[[x]])])
  names(tmp3) <- names(tmp)
  tmp <- lapply(tmp3, function(x) sort(x, decreasing = T))
  tmp <- lapply(tmp, function(x){
    if(length(x)>number) x <- x[1:number]
    return(x)
  })
  tmp2 <- lapply(tmp2, function(x) 1-x)
  tmp2 <- lapply(names(tmp), function(x) tmp2[[x]][names(tmp[[x]])])
  names(tmp2) <- names(tmp)
  tmp <- lapply(names(tmp2), function(x)
    as.data.frame(t(data.frame(row.names = names(tmp[[x]]),
                               freq = as.numeric(tmp[[x]]-(1-tmp2[[x]]))))))
  names(tmp) <- names(tmp2)
  tmp2 <- lapply(tmp2, function(x)
    as.data.frame(t(data.frame(row.names = names(x), pval=x))))


  hheatmap <- plyr::rbind.fill(tmp)
  hheatmap <- t(hheatmap)
  colnames(hheatmap) <- names(tmp)
  pval <- plyr::rbind.fill(tmp2)
  pval <- t(pval)
  colnames(pval) <- names(tmp2)

  ComplexHeatmap::pheatmap(hheatmap, cluster_rows = F, cluster_cols = F)


}


#' @import circlize
#' @export

ccircos_genLines_chr <- function(chr=c(1:22, "X", "Y", "MT"),
                                    rregion,
                                    vvalue,
                                    ttrack,
                                    ttype,
                                    aarea,
                                    ccol,
                                    bbaseline,
                                    ...){
  for (i in chr) {

    circos.genomicLines(region = rregion[which(rregion[,1] ==i) ,2:3],
                        value = vvalue[which(rregion[,1] ==i)],
                        sector.index = i,
                        track.index = ttrack,
                        type = ttype,
                        area = aarea,
                        col = ccol[which(rregion[,1] ==i)],
                        baseline = bbaseline)


  }

}


#' @import circlize
#' @export

ccircos_genPoints_chr <- function(chr=c(1:22, "X", "Y", "MT"),
                                  rregion,
                                  vvalue,
                                  ppch,
                                  ccex,
                                  ttrack,
                                  ccol,
                                  bbaseline){
  for (i in chr) {

    circos.genomicPoints(region = rregion[which(rregion[,1] ==i) ,2:3],
                         value = vvalue[which(rregion[,1] ==i)],
                         sector.index = i,
                         track.index = ttrack,
                         pch = ppch ,
                         col = ccol,
                         baseline = bbaseline,
                         cex = ccex)


  }

}

#' @import circlize
#' @importFrom gtools mixedsort
#' @export

circos_plot <- function(model_results,
                        species="hg38"){

    data <- extract_data(model_results)
    coef <- "coef_layer"%in%names(data)
    if(coef) {
      colnames(data$coef_layer) <- gsub("cov", "mean",
                                        colnames(data$coef_layer))
      data$coef_layer <- data$coef_layer[data$coef_layer$pval<=0.05,]
    }

    tmp <- lapply(names(data), function(x) {
      ans <- data[[x]]
      ans$chromosome_name <- paste0("chr", ans$chromosome_name)
      ans <- ans[ans$chromosome_name!="chrMT",]
      ans$ccol <- rep("red3", nrow(ans))
      ans$ccol[ans$mean<0] <- "royalblue4"
      return(ans)
    })
    names(tmp) <- names(data)
    data <- tmp
    wwhich <- data$res_layer$mean>quantile(data$res_layer$mean,0.99)
    data$res_layer$mean[wwhich] <- quantile(data$res_layer$mean, 0.99)
    wwhich <- data$cov_layer$mean>quantile(data$cov_layer$mean,0.99)
    data$cov_layer$mean[wwhich] <- quantile(data$cov_layer$mean, 0.99)
    tmp <- unique(gtools::mixedsort(data$res_layer$chromosome_name))
    circos.initializeWithIdeogram(species = "hg38", chromosome.index = tmp)+
    circos.genomicTrack(data, ylim=c(min(data$cov_layer$mean),
                                     max(data$cov_layer$mean)))+
    circos.genomicTrack(data, ylim=c(min(data$res_layer$mean),
                                     max(data$res_layer$mean)))+
    ccircos_genLines_chr(chr = tmp,
                         rregion = data$cov_layer[,4:6],
                         vvalue = data$cov_layer$mean,
                         ttrack = 3,
                         ttype = "h",
                         aarea = TRUE,
                         ccol = data$cov_layer$ccol,
                         bbaseline = 0,
                         lwd = 1)+
    ccircos_genLines_chr(chr = tmp,
                         rregion = data$res_layer[,4:6],
                         vvalue = data$res_layer$mean,
                         ttrack = 4,
                         ttype = "h",
                         aarea = TRUE,
                         ccol = data$res_layer$ccol,
                         bbaseline = 0,
                         lwd = 1)+
    if(coef) circos.genomicTrack(data, ylim=c(min(data$coef_layer$mean),
                                              max(data$coef_layer$mean)))+
      ccircos_genLines_chr(chr = tmp,
                           rregion = data$coef_layer[,5:7],
                           vvalue = data$coef_layer$mean,
                           ttrack = 5,
                           ttype = "h",
                           aarea = TRUE,
                           ccol = data$coef_layer$ccol,
                           bbaseline = 0,
                           lwd = 1)

}












