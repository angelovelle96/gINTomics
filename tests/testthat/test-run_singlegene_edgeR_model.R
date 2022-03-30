test_that("single gene complete integration gives correct results - Test1", {
    counts <- run_edgeR_test1_input$counts
    data <- run_edgeR_test1_input$data
    expectedres <- run_edgeR_test1_output
    design_all <- model.matrix(~W_1+W_2+W_3, data = data)
    y_all_offset <- matrix(1, ncol(counts), nrow(counts))
    y_gene_offset <- lapply(1:nrow(counts), function(x) {
        matrix(1, 1, ncol(counts))
    })
    design_single <- list()
    for (i in rownames(counts)) {
        design_single[[i]] <- model.matrix(~data[, i]+W_1+W_2+W_3, data = data)
        colnames(design_single[[i]])[2] <- "cnv"
    }
    myres <- run_edgeR_integration( response_var = t(counts),
                                    design_mat_allgene = design_all,
                                    design_mat_singlegene = design_single,
                                    offset_allgene = t(y_all_offset),
                                    offset_singlegene = y_gene_offset,
                                    threads = 16,
                                    norm_method = "TMM")
    expect_identical(myres$coef_data$cnv, expectedres$coef_ace)
    expect_identical(myres$coef_data$`(Intercept)`, expectedres$coef_intercept)
    expect_identical(myres$pval_data$cnv, expectedres$pval_ace)
    expect_identical(myres$pval_data$`(Intercept)`, expectedres$pval_intercept)
})
############

test_that("single gene complete TF integration gives correct results - Test2", {
    mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
    tf_expression_model <- run_edgeR_test2_input$tf_expression_model
    expanded_tf_mirna_couples <- run_edgeR_test2_input$expanded_tf_mirna_couples
    expectedres <- run_edgeR_test2_output
    design_single <- lapply(1:ncol(mirna_exp_model), function(y) {
        tf <- expanded_tf_mirna_couples[[y]][, "tf"]
        tmp <- gsub("-", "_", tf)
        tmp <- formula(paste0("~", paste0(tmp, collapse = "+")))
        tmp2 <- tf_expression_model
        colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
        design <- model.matrix(tmp, data = tmp2)
        colnames(design)[2:ncol(design)] <- tf
        return(design)
    })

    y_all_offset <- matrix(1, ncol(mirna_exp_model), nrow(mirna_exp_model))
    y_gene_offset <- lapply(1:nrow(mirna_exp_model), function(x) {
        matrix(1, 1, ncol(mirna_exp_model))
    })
    design_all <- model.matrix(~1, data = tf_expression_model)
    myres <- run_edgeR_integration( response_var = mirna_exp_model,
                                    design_mat_allgene = design_all,
                                    design_mat_singlegene = design_single,
                                    threads = 16,
                                    norm_method = "TMM")

    colnames(myres$coef_data) <- paste0(colnames(myres$coef_data), "_coef")
    colnames(myres$pval_data) <- paste0(colnames(myres$pval_data), "_pvalue")
    tmp <- cbind(myres$coef_data, myres$pval_data)
    tmp <- tmp[rownames(expectedres), colnames(expectedres)]
    expect_identical(tmp, expectedres)
})
