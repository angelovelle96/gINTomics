singlegene_edgeR_model <- function(response_var,
    covariates, design_mat = NULL, offset_singlegene = NULL,
    y_all, threads = 1) {

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
            y_gene <- DGEList(counts = t(response_var[x,
                ]))
            if (!is.null(offset_singlegene)) {
                if (is.list(offset_singlegene)) {
                  y_gene$offset <- offset_singlegene[[x]]
                } else {
                  y_gene$offset <- offset_singlegene
                }
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
