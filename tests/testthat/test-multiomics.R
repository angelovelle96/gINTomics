library(BiocParallel)

test_that("run_cnv_integration works", {

  data("test_results")
  counts <- run_edgeR_test1_input$counts
  data <- run_edgeR_test1_input$data
  expectedres <- run_edgeR_test1_output
  design_all <- model.matrix(~W_1+W_2+W_3, data = data)
  y_all_offset <- matrix(1, ncol(counts), nrow(counts))
  y_gene_offset <- lapply(1:nrow(counts), function(x) {
    matrix(1, 1, ncol(counts))
  })

  myres <- .run_edgeR_integration( response_var = t(counts),
                                  covariates = data,
                                  steady_covariates = c("W_1", "W_2", "W_3"),
                                  design_mat_allgene = design_all,
                                  offset_allgene = t(y_all_offset),
                                  offset_singlegene = y_gene_offset,
                                  norm_method = "TMM")

  myres2 <- run_cnv_integration(expression = t(counts),
                                cnv_data = data,
                                sequencing_data = T,
                                steady_covariates = c("W_1", "W_2", "W_3"),
                                design_mat_allgene = design_all,
                                offset_allgene = t(y_all_offset),
                                offset_singlegene = y_gene_offset,
                                norm_method = "TMM")

    expect_identical(myres$coef_data, myres2$coef_data)

})

# test_that("run_tf_integration sequencing data check", {
#
#   data("test_results")
#   mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
#   tf_expression_model <- run_edgeR_test2_input$tf_expression_model
#   interactions <- run_edgeR_test2_input$expanded_tf_mirna_couples
#   interactions <- lapply(interactions, function(x) x[,1])
#
#
#   myres <- suppressWarnings(.run_edgeR_integration(
#                                   response_var = mirna_exp_model,
#                                   interactions = interactions,
#                                   covariates = tf_expression_model,
#                                   norm_method = "TMM"))
#
#   myres2 <- suppressWarnings(run_tf_integration(expression = mirna_exp_model,
#                                 tf_expression = tf_expression_model,
#                                 interactions = interactions,
#                                 sequencing_data = T,
#                                 norm_method = "TMM", normalize_cov = F))
#
#   expect_identical(myres$model_results$`hsa-miR-29c-3p;iso_3p:a`,
#                    myres2$model_results$`hsa-miR-29c-3p;iso_3p:a`)
#   expect_identical(names(myres$model_results), names(myres2$model_results))
# })
#
#
# test_that("run_tf_integration microarray data check", {
#
#   data("test_results")
#   mirna_exp_model <- run_edgeR_test2_input$mirna_exp_model
#   tf_expression_model <- run_edgeR_test2_input$tf_expression_model
#   interactions <- run_edgeR_test2_input$expanded_tf_mirna_couples
#   interactions <- lapply(interactions, function(x) x[,1])
#
#
#   myres <- suppressWarnings(.run_lm_integration(response_var = mirna_exp_model,
#                                  interactions = interactions,
#                                  covariates = tf_expression_model))
#
#   myres2 <- suppressWarnings(run_tf_integration(expression = mirna_exp_model,
#                                tf_expression = tf_expression_model,
#                                interactions = interactions,
#                                sequencing_data = F, normalize_cov = F))
#
#
#   expect_identical(myres$model_results[[5]]$coefficients,
#                    myres2$model_results[[5]]$coefficients)
#
# })


test_that("run_multiomics works", {
  data("ov_test_tcga_omics")
 multiomics_integration <-run_multiomics(data = mmultiassay_ov)
  expect_true(is(multiomics_integration, "MultiOmics"),
              "The function should return a MultiOmics object")
})


test_that(".def_cnv_integration works", {
   data("ov_test_tcga_omics")
  gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
  gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
  expect_no_error({
    cnv_res <- .def_cnv_integration(expression = gene_exp_matrix,
                                    cnv_data = gene_cnv_matrix,
                                    sequencing_data = TRUE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    BPPARAM = SerialParam())
    expect_type(cnv_res, "list")
  })
  expect_no_error({
    cnv_res <- .def_cnv_integration(expression = gene_exp_matrix,
                                    cnv_data = gene_cnv_matrix,
                                    sequencing_data = FALSE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    BPPARAM =  SerialParam())
    expect_type(cnv_res, "list")
  })
})


test_that("run_cnv_integration function works", {
  data("ov_test_tcga_omics")
  gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
  gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
  class <- rep(c('A', 'B'),
             each = 10)
names(class) <- colnames(mmultiassay_ov[[1]])
  expect_no_error({
    cnv_res <- run_cnv_integration(expression = gene_exp_matrix,
                                    cnv_data = gene_cnv_matrix,
                                    sequencing_data = TRUE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    class = NULL,
                                    run_deg = TRUE,
                                    BPPARAM = SerialParam())
    expect_type(cnv_res, "list")
  })
  expect_no_error({
    cnv_res <- run_cnv_integration(expression = gene_exp_matrix,
                                    cnv_data = gene_cnv_matrix,
                                    sequencing_data = FALSE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    run_deg = TRUE,
                                    BPPARAM = SerialParam())
    expect_type(cnv_res, "list")
  })
})

test_that(".def_met_integration works", {
   data("ov_test_tcga_omics")
  gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
  expect_no_error({
    met_res <- .def_met_integration(expression = gene_exp_matrix,
                                    met_data = gene_met_matrix,
                                    sequencing_data = TRUE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    BPPARAM = SerialParam())
    expect_type(met_res, "list")
  })
  expect_no_error({
    met_res <- .def_met_integration(expression = gene_exp_matrix,
                                    met_data = gene_met_matrix,
                                    sequencing_data = FALSE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    BPPARAM =  SerialParam())
    expect_type(met_res, "list")
  })
})

test_that("run_met_integration function works", {
  data("ov_test_tcga_omics")
  gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
  class <- rep(c('A', 'B'),
             each = 10)
names(class) <- colnames(mmultiassay_ov[[1]])
  expect_no_error({
    met_res <- run_met_integration(expression = gene_exp_matrix,
                                    met_data = gene_met_matrix,
                                    sequencing_data = TRUE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    class = NULL,
                                    run_deg = TRUE,
                                    BPPARAM = SerialParam())
    expect_type(met_res, "list")
  })
  expect_no_error({
    met_res <- run_cnv_integration(expression = gene_exp_matrix,
                                    met_data = gene_met_matrix,
                                    sequencing_data = FALSE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    run_deg = TRUE,
                                    BPPARAM = SerialParam())
    expect_type(met_res, "list")
  })
})

test_that(".def_genomic_integration works", {
  data("ov_test_tcga_omics")
  gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
  gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
  class <- rep(c('A', 'B'),
             each = 10)
names(class) <- colnames(mmultiassay_ov[[1]])
  gen_res <- .def_genomic_integration(expression = gene_exp_matrix,
                                      cnv_data = gene_cnv_matrix,
                                      methylation = gene_met_matrix,
                                      sequencing_data = TRUE,
                                      normalize = TRUE,
                                      norm_method = "TMM",
                                      interactions = NULL,
                                      scale = TRUE,
                                      BPPARAM = SerialParam())
  expect_type(gen_res, "list")
})

test_that("run_genomic_integration works", {
  gene_exp_matrix <- as.matrix(assay(mmultiassay_ov[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(mmultiassay_ov[["methylation"]]))
  gene_cnv_matrix <- as.matrix(assay(mmultiassay_ov[["cnv_data"]]))
  gen_res <- run_genomic_integration(expression = gene_exp_matrix,
                                     cnv_data = gene_cnv_matrix,
                                     methylation = gene_met_matrix,
                                     sequencing_data = TRUE,
                                     normalize = TRUE,
                                     norm_method = "TMM",
                                     interactions = NULL,
                                     scale = TRUE,
                                     run_deg = TRUE,
                                     BPPARAM = SerialParam())
  expect_is(gen_res, "list",
            info = "The function should return a list")
})

# Run the tests
test_check("run_genomic_integration")


