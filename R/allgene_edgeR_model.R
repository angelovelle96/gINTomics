allgene_edgeR_model <- function(response_var,
                                covariates = NULL,
                                design_mat_allgene = NULL,
                                offset_allgene = NULL,
                                norm_method = "TMM") {

    if(is.null(covariates)&is.null(design_mat_allgene)) {
        stop("Design matrix should be supplied if covariates are not specified")
    }

    if(!norm_method %in% c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
        stop(str_wrap("norm_method should be one of TMM, TMMwsp, RLE,
                        upperquartile, none"))
    }
    if(is.null(design_mat_allgene)) {
        design_mat_allgene <- model.matrix(~1, data = covariates)
    }
    if(ncol(response_var) != nrow(design_mat_allgene)) {
        stop("Number of samples differ between respone_var and design_mat")
    }
    y_all <- DGEList(counts = response_var)
    if (!is.null(offset_allgene)) {
        y_all$offset <- offset_allgene
    } else {
        y_all <- calcNormFactors(y_all, method = norm_method)
    }
    y_all <- estimateGLMCommonDisp(y_all, design_mat_allgene)
    y_all <- estimateGLMTagwiseDisp(y_all, design_mat_allgene)
    return(y_all)
}
