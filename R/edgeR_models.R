# Running edgeR integration
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom plyr rbind.fill
#' @importFrom stats residuals
#' @importFrom BiocParallel bpparam SerialParam
.run_edgeR_integration <- function(response_var,
    covariates, interactions = "auto", design_mat_allgene = NULL,
    offset_allgene = NULL, offset_singlegene = NULL,
    normalize = TRUE, norm_method = "TMM", steady_covariates = NULL,
    reference = NULL, BPPARAM = SerialParam()) {
    tmp <- .data_check(response_var = response_var,
        covariates = covariates, interactions = interactions)
    if (is.null(tmp))
        return(NULL)
    response_var <- tmp$response_var
    covariates <- tmp$covariates
    interactions <- tmp$interactions
    tmp <- unlist(lapply(interactions, length))
    single_cov <- FALSE
    if (sum(tmp == 1) == length(tmp)) {
        single_cov <- TRUE
    }
    if (normalize == FALSE)
        offset_allgene <- matrix(1, ncol(response_var),
            nrow(response_var))
    fit_all <- .allgene_edgeR_model(response_var = response_var,
        design_mat_allgene = design_mat_allgene,
        offset_allgene = offset_allgene, norm_method = norm_method)
    tmp <- .covariates_check(response_var = response_var,
        covariates = covariates, interactions = interactions,
        steady_covariates = steady_covariates, reference = reference)
    if (is.null(tmp))
        return(NULL)
    covariates <- tmp$covariates
    response_var <- tmp$response_var
    interactions <- tmp$interactions
    original_id <- tmp$original_id
    fformula <- .generate_formula(interactions = interactions)
    fit_gene <- .singlegene_edgeR_model(response_var = response_var,
        covariates = covariates, fit_all = fit_all,
        fformula = fformula, offset_singlegene = offset_singlegene,
        BPPARAM = BPPARAM)
    model_res <- .edger_coef_test(fit_gene, BPPARAM = BPPARAM)
    coef_pval_mat <- .building_result_matrices(model_results = model_res,
        type = "edgeR", single_cov = single_cov)
    tmp <- bplapply(fit_gene, function(x) as.data.frame(residuals(x)),
        BPPARAM = BPPARAM)
    rresiduals <- rbind.fill(tmp)
    colnames(rresiduals) <- rownames(response_var)
    rownames(rresiduals) <- names(tmp)
    results <- list(coef_data = coef_pval_mat$coef,
        pval_data = coef_pval_mat$pval, fdr_data = fdr(coef_pval_mat$pval),
        residuals = rresiduals, data = list(response_var = response_var,
            covariates = covariates, formula = fformula))
    if (nrow(original_id) > 0) {
        results <- .id_conversion(dictionary = original_id,
            results = results)
    }
    return(results)
}
# All edgeR gene models
#' @importFrom stats model.matrix
#' @importFrom edgeR glmFit DGEList calcNormFactors estimateGLMCommonDisp
#'  estimateGLMTagwiseDisp
.allgene_edgeR_model <- function(response_var, design_mat_allgene = NULL,
                                 offset_allgene = NULL,
    norm_method = "TMM") {
    response_var <- t(response_var + 1)
    if (!norm_method %in% c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
        stop(str_wrap("norm_method should be one of TMM, TMMwsp, RLE,
                        upperquartile, none"))
    }
    if (is.null(design_mat_allgene)) {
        design_mat_allgene <- model.matrix(~1,
                                           data = as.data.frame(t(
                                             response_var)))
    }
    if (ncol(response_var) != nrow(design_mat_allgene)) {
        stop(str_wrap("Number of samples differ between respone_var and
                        design_mat_allgene"))
    }
    y_all <- DGEList(counts = response_var)
    y_all <- calcNormFactors(y_all, method = norm_method)
    if (!is.null(offset_allgene))
        y_all$offset <- offset_allgene
    y_all <- estimateGLMCommonDisp(y_all, design_mat_allgene)
    y_all <- estimateGLMTagwiseDisp(y_all, design_mat_allgene)
    fit_all <- glmFit(y_all, design_mat_allgene)
    fit_all$common.dispersion <- y_all$common.dispersion
    fit_all$tagwise.dispersion <- y_all$tagwise.dispersion
    return(fit_all)
}
# Single gene edgeR model
#' @importFrom BiocParallel bplapply
.singlegene_edgeR_model <- function(response_var, covariates,
                                    offset_singlegene = NULL,
    fformula, fit_all, BPPARAM) {
    response_var <- t(response_var + 1)
    fformula <- lapply(seq_along(fformula), function(x) {
        list(formula = fformula[[x]], var_name = names(fformula)[x])
    })
    names(fformula) <- vapply(fformula, function(x) x$var_name, FUN.VALUE = "A")
    if (!is.null(offset_singlegene)) {
        if (is.list(offset_singlegene)) {
            if (length(offset_singlegene) != nrow(response_var))
                stop(str_wrap("If provided as list, the length of
                              offset_singlegene should be equal to the number
                              of genes in response_var"))
        }
    }
    fit_list <- bplapply(fformula, .def_edger, response_var = response_var,
                         covariates = covariates,
        fit_all = fit_all, offset_singlegene = offset_singlegene,
        BPPARAM = BPPARAM)
    names(fit_list) <- rownames(response_var)
    return(fit_list)
}
# edger
#' @importFrom stats model.matrix
#' @importFrom edgeR glmFit DGEList
#' @importFrom stats as.formula terms
.def_edger <- function(formula, response_var, covariates,
    offset_singlegene = NULL, fit_all) {
    if (length(attr(terms(formula$formula), "term.labels")) >
        as.integer(ncol(response_var) * 0.4)) {
        bad_str <- paste0(c("-", ";", ":", "\\*", "%in%",
            "\\^"), collapse = "|")
        tmp <- t(response_var)
        colnames(tmp) <- gsub(bad_str, "_", colnames(tmp))
        tmp <- cbind(tmp, covariates)
        tmp2 <- gsub(bad_str, "_", formula$var_name)
        tmp2 <- as.formula(paste(as.character(tmp2), "~",
            as.character(formula$formula[2])))
        tmp <- .rf_selection(formula = tmp2, data = tmp)
        formula$formula <- tmp
    }
    design <- model.matrix(formula$formula, data = covariates)
    dge <- DGEList(counts = t(response_var[formula$var_name,
        ]))
    if (!is.null(offset_singlegene)) {
        if (is.list(offset_singlegene)) {
            dge$offset <- offset_singlegene[[formula$var_name]]
        } else {
            dge$offset <- offset_singlegene
        }
    } else {
        dge$offset <- fit_all$offset[1, ]
    }
    dge$samples$norm.factors <- fit_all$samples$norm.factors
    dge$common.dispersion <- fit_all$common.dispersion
    tmp <- which(rownames(fit_all$coefficients) == formula$var_name)
    dge$tagwise.dispersion <- fit_all$tagwise.dispersion[tmp]
    fit <- glmFit(dge, design)
    return(fit)
}
# edgeR coefficient test
#' @importFrom BiocParallel bplapply
.edger_coef_test <- function(fit_list, BPPARAM) {
    top_list <- bplapply(fit_list, .def_coef_test, BPPARAM = BPPARAM)
    return(top_list)
}
# coefficient test
#' @importFrom edgeR glmLRT
#' @importFrom plyr rbind.fill
.def_coef_test <- function(fit) {
    coef <- colnames(fit$coefficients)
    lrt <- lapply(coef, function(z) {
        glmLRT(fit, coef = z)
    })
    names(lrt) <- coef
    top <- lapply(names(lrt), function(z) topTags(lrt[[z]])$table)
    top <- rbind.fill(top)
    rownames(top) <- names(lrt)
    top <- as.matrix(t(top))
    coef <- vapply(names(lrt), function(z) {
        lrt[[z]]$coefficients[, z]
    }, FUN.VALUE = numeric(1))
    top <- rbind(top, coef)
    top <- as.data.frame(top)
    return(top)
}
