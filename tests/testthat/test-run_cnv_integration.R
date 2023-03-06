test_that("run_cnv_integration check", {

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

  myres2 <- run_cnv_integration(expression = t(counts),
                                cnv_data = data,
                                sequencing_data = T,
                                steady_covariates = c("W_1", "W_2", "W_3"),
                                design_mat_allgene = design_all,
                                offset_allgene = t(y_all_offset),
                                offset_singlegene = y_gene_offset,
                                threads = 16,
                                norm_method = "TMM")

    expect_identical(myres, myres2)

})
