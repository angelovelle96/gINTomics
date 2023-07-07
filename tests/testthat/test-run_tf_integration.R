test_that("run_tf_integration sequencing data check", {

  mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
  tf_expression_model <- run_edgeR_test2_input$tf_expression_model
  interactions <- run_edgeR_test2_input$expanded_tf_mirna_couples
  interactions <- lapply(interactions, function(x) x[,1])
  expectedres <- run_edgeR_test2_output


  myres <- run_edgeR_integration(response_var = mirna_exp_model,
                                  interactions = interactions,
                                  covariates = tf_expression_model,
                                  norm_method = "TMM")

  myres2 <- run_tf_integration(expression = mirna_exp_model,
                                tf_expression = tf_expression_model,
                                interactions = interactions,
                                sequencing_data = T,
                                norm_method = "TMM", normalize_cov = F)

  expect_identical(myres$coef_data, myres2$coef_data)
})


test_that("run_tf_integration microarray data check", {

  mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
  tf_expression_model <- run_edgeR_test2_input$tf_expression_model
  interactions <- run_edgeR_test2_input$expanded_tf_mirna_couples
  interactions <- lapply(interactions, function(x) x[,1])
  expectedres <- run_edgeR_test2_output


  myres <- suppressWarnings(run_lm_integration(response_var = mirna_exp_model,
                                 interactions = interactions,
                                 covariates = tf_expression_model,
                                 step = T))

  myres2 <- suppressWarnings(run_tf_integration(expression = mirna_exp_model,
                               tf_expression = tf_expression_model,
                               interactions = interactions,
                               sequencing_data = F,
                               step=T, normalize_cov = F))

  tmp <- lapply(myres, function(x)  x$coefficients)
  tmp2 <- lapply(myres2, function(x) x$coefficients)
  expect_identical(tmp, tmp2)

})
