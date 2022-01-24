offset <- matrix(1, 1, nrow(mirna_exp_model))
a <- lapply(1:ncol(mirna_exp_model), function(y) {
    tf <- expanded_tf_mirna_couples[[y]][,
        "tf"]
    tmp <- gsub("-", "_", tf)
    tmp <- formula(paste0("~", paste0(tmp,
        collapse = "+")))
    tmp2 <- tf_expression_model
    colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
    design <- model.matrix(tmp, data = tmp2)
    colnames(design)[2:ncol(design)] <- tf
    return(design)
})
b <- run_singlegene_edgeR_integration(response_var = mirna_exp_model,
    covariates = tf_expression_model, design_mat_singlegene = a,
    offset_singlegene = offset, threads = 16)
