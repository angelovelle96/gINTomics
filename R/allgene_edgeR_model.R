allgene_edgeR_model <- function(response_var,
    covariates = NULL, design_mat_allgene = NULL,
    offset_allgene = NULL, norm_method_all = "TMM") {

    if (is.null(covariates) & is.null(design_mat_allgene))
        stop("Design matrix should be supplied if covariates are not specified")

    if (is.null(design_mat_allgene)) {
        design_mat_allgene <- model.matrix(~1,
            data = covariates)
    }
    y_all <- DGEList(counts = response_var)
    if (!is.null(offset_allgene)) {
        y_all$offset <- offset_allgene
    } else {
        y_all <- calcNormFactors(y_all, method = norm_method_all)
    }
    y_all <- estimateGLMCommonDisp(y_all,
        design_mat_allgene)
    y_all <- estimateGLMTagwiseDisp(y_all,
        design_mat_allgene)
    return(y_all)

}
