test_that("single gene complete integration gives correct results - Test1", {
    data("test_results")
    counts <- run_edgeR_test1_input$counts
    data <- run_edgeR_test1_input$data
    expectedres <- run_edgeR_test1_output
    design_all <- model.matrix(~W_1+W_2+W_3, data = data)
    y_all_offset <- matrix(1, ncol(counts), nrow(counts))
    y_gene_offset <- lapply(1:nrow(counts), function(x) {
        matrix(1, 1, ncol(counts))
    })
    names(y_gene_offset) <- rownames(counts)

    myres <- run_edgeR_integration( response_var = t(counts),
                                    covariates = data,
                                    steady_covariates = c("W_1", "W_2", "W_3"),
                                    design_mat_allgene = design_all,
                                    offset_allgene = t(y_all_offset),
                                    offset_singlegene = y_gene_offset,
                                    normalize = F)
    expect_identical(myres$coef_data$cov, expectedres$coef_ace)
    expect_identical(myres$coef_data$`(Intercept)`, expectedres$coef_intercept)
    expect_equal(myres$pval_data$cov, expectedres$pval_ace, tolerance = 10^-12)
    expect_equal(myres$pval_data$`(Intercept)`, expectedres$pval_intercept,
                 tolerance = 10^-12)

})
