allgene_edgeR_model <- function(response_var,covariates, design_mat_allgene=NULL, offset_allgene=NULL){


if(is.null(design_mat_allgene)){
  design_mat_allgene <- model.matrix(~1, data=covariates)
  }
y_all <- DGEList(counts=response_var)
if (!is.null(offset_allgene)) {
  y_all$offset <- offset_allgene
}
y_all <- estimateGLMCommonDisp(y_all, design_mat_allgene)
y_all <- estimateGLMTagwiseDisp(y_all, design_mat_allgene)
return(y_all)

}
