
extract_model_res <-  function(object, ...) UseMethod("extract_model_res")


extract_model_res.default <- function(object,
                                      outliers=F){

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
  tmp <- rep("not_sign", nrow(pval))
  tmp[pval$value<=0.05] <- "sign"
  names(tmp) <- rownames(pval)
  data$pval <- tmp[rownames(data)]
  mmin <- quantile(data$value, 0.25) - 1.5*IQR(data$value)
  mmax <- quantile(data$value, 0.75) + 1.5*IQR(data$value)

  if(outliers==F){
    data <- data[data$value>mmin,]
    data <- data[data$value<mmax,]
  }

  return(data)

}





#' @import ggplot2 ggridges
#' @importFrom reshape melt

ridgeline_plot <- function(model_results,
                           outliers=F){


    ggplot(data, aes(x = value, y = "pval", fill="pval"))+
      geom_density_ridges() +
      theme_ridges()


}


