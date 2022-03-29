#' Single gene omics integration
#'
#' @description Compute single gene edgeR's models for each of the genes
#' provided in **response_var**
#' @param response_var Matrix or data.frame containing the expression values
#' for each model. Rows represent samples, while each column represents
#' the different response variables of the models. The number of rows should
#' be identical to that of **covariates** or of design matrices provided in
#' **design_mat_singlegene** and the order of the samples should be the same.
#' In case you need to run a single model you can provide the response
#' variable as an atomic vector, anyway the response variables should always
#' be numeric values.
#' @param covariates Matrix or data.frame containing the covariates of the
#' models. Rows represent samples, while columns represent the different
#' the covariates. All the covariates contained in columns will be used for
#' each model, so use this argument only if all the models have the same
#' covariates. In case you need a single covariate you can provide the
#' covariate as an atomic vector. This argument is ignored if
#' **design_mat_singlegene** is provided.
#' @param design_mat A design matrix (if the model to run is one) or
#' a list of design matrices for each edgeR single gene model. The number of
#' design matrices should be equal to the number of rows of **response_var**.
#' @param offset_singlegene A vector (if the model to run is one) or
#' a list of vector containing the offsets for each single gene edgeR model.
#' The number of vectors should be equal to the number of rows
#' of **response_var**.
#' @param y_all An all gene edgeR model needed to extract the common and
#' tagwise dispersion
#' @param threads Number of threads to use for parallelization (Default is 1)
#' @import parallel edgeR

singlegene_edgeR_model <- function(response_var,
    covariates = NULL, design_mat = NULL,
    offset_singlegene = NULL, y_all, threads = 1) {

    if (is.null(covariates) & is.null(design_mat))
        stop("Design matrix should be supplied if covariates are not specified")

    if (is.null(design_mat)) {
        cov <- colnames(covariates)
        cov <- gsub("-", "_", cov)
        cov <- formula(paste0("~", paste0(cov,
            collapse = "+")))
        tmp <- covariates
        colnames(tmp) <- gsub("-", "_", colnames(tmp))
        design_mat <- model.matrix(cov, data = tmp)
        colnames(design_mat)[2:ncol(design_mat)] <- colnames(covariates)
    }

    fit_list <- mclapply(1:nrow(response_var),
        function(x) {
            if(ncol(response_var)!=nrow(design_mat[[x]]))
                stop(paste("Number of samples differ between respone_var",
                           "and design_mat"))
            y_gene <- DGEList(counts = t(response_var[x,
                ]))
            if (!is.null(offset_singlegene)) {
                if (is.list(offset_singlegene)) {
                  y_gene$offset <- offset_singlegene[[x]]
                } else {
                  y_gene$offset <- offset_singlegene
                }
            } else {
                y_gene$samples$norm.factors <- y_all$samples$norm.factors
            }
            y_gene$common.dispersion <- y_all$common.dispersion
            y_gene$tagwise.dispersion <- y_all$tagwise.dispersion[x]
            if (is.list(design_mat)) {
                design <- design_mat[[x]]
                fit <- glmFit(y_gene, design)
            } else {
                fit <- glmFit(y_gene, design_mat)
            }
            return(fit)
        }, mc.cores = threads)
    names(fit_list) <- rownames(response_var)
    return(fit_list)
}
