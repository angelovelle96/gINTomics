test_that("single gene complete integration gives correct results - Test1", {
    counts <- run_edgeR_test1_input$counts
    data <- run_edgeR_test1_input$data
    expectedres <- run_edgeR_test1_output
    design_all <- model.matrix(~W_1+W_2+W_3, data = data)
    y_all_offset <- matrix(1, ncol(counts), nrow(counts))
    y_gene_offset <- lapply(1:nrow(counts), function(x) {
        matrix(1, 1, ncol(counts))
    })

    myres <- run_edgeR_integration( response_var = t(counts),
                                    covariates = data,
                                    steady_covariates = c("W_1", "W_2", "W_3"),
                                    design_mat_allgene = design_all,
                                    offset_allgene = t(y_all_offset),
                                    offset_singlegene = y_gene_offset,
                                    threads = 16,
                                    norm_method = "TMM",
                                    cnv_mode = T)
    expect_identical(myres$coef_data$cnv, expectedres$coef_ace)
    expect_identical(myres$coef_data$`(Intercept)`, expectedres$coef_intercept)
    expect_identical(myres$pval_data$cnv, expectedres$pval_ace)
    expect_identical(myres$pval_data$`(Intercept)`, expectedres$pval_intercept)

    design_single <- list()
    for (i in rownames(counts)) {
        design_single[[i]] <- model.matrix(~data[, i]+W_1+W_2+W_3, data = data)
        colnames(design_single[[i]])[2] <- "cnv"
    }

    myres2 <- run_edgeR_integration( response_var = t(counts),
                                    design_mat_singlegene = design_single,
                                    design_mat_allgene = design_all,
                                    offset_allgene = t(y_all_offset),
                                    offset_singlegene = y_gene_offset,
                                    threads = 16,
                                    norm_method = "TMM",
                                    cnv_mode = T)
    expect_identical(myres2$coef_data$cnv, expectedres$coef_ace)
    expect_identical(myres2$coef_data$`(Intercept)`, expectedres$coef_intercept)
    expect_identical(myres2$pval_data$cnv, expectedres$pval_ace)
    expect_identical(myres2$pval_data$`(Intercept)`, expectedres$pval_intercept)
})
############

test_that("single gene complete TF integration gives correct results - Test2", {
    mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
    tf_expression_model <- run_edgeR_test2_input$tf_expression_model
    interactions <- run_edgeR_test2_input$expanded_tf_mirna_couples
    interactions <- lapply(interactions, function(x) x[,1])
    expectedres <- run_edgeR_test2_output
    design_single <- lapply(1:ncol(mirna_exp_model), function(y) {
        tf <- interactions[[y]]
        tmp <- gsub("-", "_", tf)
        tmp <- formula(paste0("~", paste0(tmp, collapse = "+")))
        tmp2 <- tf_expression_model
        colnames(tmp2) <- gsub("-", "_", colnames(tmp2))
        design <- model.matrix(tmp, data = tmp2)
        colnames(design)[2:ncol(design)] <- tf
        return(design)
    })


    myres <- run_edgeR_integration( response_var = mirna_exp_model,
                                    covariates = tf_expression_model,
                                    design_mat_singlegene = design_single,
                                    threads = 16,
                                    norm_method = "TMM")
    myres2 <- run_edgeR_integration(response_var = mirna_exp_model,
                                    interactions = interactions,
                                    covariates = tf_expression_model,
                                    threads = 16,
                                    norm_method = "TMM")

    colnames(myres$coef_data) <- paste0(colnames(myres$coef_data), "_coef")
    colnames(myres$pval_data) <- paste0(colnames(myres$pval_data), "_pvalue")
    tmp <- cbind(myres$coef_data, myres$pval_data)
    tmp <- tmp[rownames(expectedres), colnames(expectedres)]
    colnames(myres2$coef_data) <- paste0(colnames(myres2$coef_data), "_coef")
    colnames(myres2$pval_data) <- paste0(colnames(myres2$pval_data), "_pvalue")
    tmp2 <- cbind(myres2$coef_data, myres2$pval_data)
    tmp2 <- tmp2[rownames(expectedres), colnames(expectedres)]
    expect_identical(tmp, expectedres)
    expect_identical(tmp2, expectedres)
    expect_identical(myres, myres2)
})
