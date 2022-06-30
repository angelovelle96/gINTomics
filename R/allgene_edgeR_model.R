allgene_edgeR_model <- function(response_var,
                                design_mat_allgene = NULL,
                                offset_allgene = NULL,
                                norm_method = "TMM") {

    response_var <- t(response_var+1)
    if(!norm_method %in% c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
        stop(str_wrap("norm_method should be one of TMM, TMMwsp, RLE,
                        upperquartile, none"))
    }

    if(is.null(design_mat_allgene)) {
        design_mat_allgene <- model.matrix(~1,
                data = as.data.frame(t(response_var)))
    }
    if(ncol(response_var) != nrow(design_mat_allgene)) {
        stop(str_wrap("Number of samples differ between respone_var and
                        design_mat_allgene"))
    }
    y_all <- DGEList(counts = response_var)
    y_all <- calcNormFactors(y_all, method = norm_method)
    if (!is.null(offset_allgene)) y_all$offset <- offset_allgene
    y_all <- estimateGLMCommonDisp(y_all, design_mat_allgene)
    y_all <- estimateGLMTagwiseDisp(y_all, design_mat_allgene)
    fit_all <- glmFit(y_all, design_mat_allgene)
    fit_all$common.dispersion <- y_all$common.dispersion
    fit_all$tagwise.dispersion <- y_all$tagwise.dispersion
    return(fit_all)
}
