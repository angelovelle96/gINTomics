library(BiocParallel)

test_that("run_multiomics works", {
  data <- data_shiny_tests$multiassay
  multiomics_integration <-run_multiomics(data = data)
  expect_true(is(multiomics_integration, "MultiOmics"),
              "The function should return a MultiOmics object")
})


test_that(".def_cnv_integration works", {
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(assay(data[["gene_exp"]]))
  gene_cnv_matrix <- as.matrix(assay(data[["cnv_data"]]))
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
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(assay(data[["gene_exp"]]))
  gene_cnv_matrix <- as.matrix(assay(data[["cnv_data"]]))
  class <- rep(c('A', 'B'), each = 10)
  names(class) <- colnames(data[[1]])
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
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(assay(data[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(data[["methylation"]]))
  expect_no_error({
    met_res <- .def_met_integration(expression = gene_exp_matrix,
                                    methylation = gene_met_matrix,
                                    sequencing_data = TRUE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    BPPARAM = SerialParam())
    expect_type(met_res, "list")
  })
  expect_no_error({
    met_res <- .def_met_integration(expression = gene_exp_matrix,
                                    methylation = gene_met_matrix,
                                    sequencing_data = FALSE,
                                    normalize = TRUE,
                                    norm_method = "TMM",
                                    BPPARAM =  SerialParam())
    expect_type(met_res, "list")
  })
})

test_that("run_met_integration function works", {
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(assay(data[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(data[["methylation"]]))
  class <- rep(c('A', 'B'), each = 10)
names(class) <- colnames(data[[1]])
  expect_no_error({
    met_res <- run_met_integration(expression = gene_exp_matrix,
                                   methylation = gene_met_matrix,
                                   sequencing_data = TRUE,
                                   normalize = TRUE,
                                   norm_method = "TMM",
                                   class = NULL,
                                   run_deg = TRUE,
                                   BPPARAM = SerialParam())
    expect_type(met_res, "list")
  })
  expect_no_error({
    met_res <- run_met_integration(expression = gene_exp_matrix,
                                   methylation = gene_met_matrix,
                                   sequencing_data = FALSE,
                                   normalize = TRUE,
                                   norm_method = "TMM",
                                   run_deg = TRUE,
                                   BPPARAM = SerialParam())
    expect_type(met_res, "list")
  })
})

test_that(".def_genomic_integration works", {
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(assay(data[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(data[["methylation"]]))
  gene_cnv_matrix <- as.matrix(assay(data[["cnv_data"]]))
  class <- rep(c('A', 'B'), each = 10)
  names(class) <- colnames(data[[1]])
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
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(assay(data[["gene_exp"]]))
  gene_met_matrix <- as.matrix(assay(data[["methylation"]]))
  gene_cnv_matrix <- as.matrix(assay(data[["cnv_data"]]))
  tested <- run_genomic_integration(expression = gene_exp_matrix,
                                     cnv_data = gene_cnv_matrix,
                                     methylation = gene_met_matrix,
                                     sequencing_data = TRUE,
                                     normalize = TRUE,
                                     norm_method = "TMM",
                                     interactions = NULL,
                                     scale = TRUE,
                                     run_deg = TRUE,
                                     BPPARAM = SerialParam())
  expect_type(tested, "list")
})


test_that(".def_tf_integration works", {
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(t(assay(data[["gene_exp"]])))
  tested <- .def_tf_integration(expression=gene_exp_matrix,
                             tf_expression=gene_exp_matrix,
                             interactions=NULL,
                             type="tf",
                             sequencing_data=TRUE,
                             species="hsa",
                             normalize=TRUE,
                             norm_method="TMM",
                             normalize_cov=TRUE,
                             norm_method_cov="TMM",
                             BPPARAM=SerialParam())
  expect_type(tested, "list")
  tested <- .def_tf_integration(expression=gene_exp_matrix,
                             tf_expression=gene_exp_matrix,
                             interactions=NULL,
                             type="tf",
                             sequencing_data=FALSE,
                             species="hsa",
                             normalize=TRUE,
                             norm_method="TMM",
                             normalize_cov=TRUE,
                             norm_method_cov="TMM",
                             BPPARAM=SerialParam())
  expect_type(tested, "list")
  expect_named(tested$coef_data, c("(Intercept)", "cov"))
})


test_that("run_tf_integration works", {
  data <- data_shiny_tests$multiassay
  gene_exp_matrix <- as.matrix(t(assay(data[["gene_exp"]])))
  class <- rep(c('A', 'B'), each = 10)
  names(class) <- colnames(data[[1]])
  tested <- run_tf_integration(expression = gene_exp_matrix,
                            tf_expression = gene_exp_matrix,
                            interactions = NULL,
                            type = "tf",
                            sequencing_data = TRUE,
                            species = "hsa",
                            normalize = TRUE,
                            norm_method = "TMM",
                            normalize_cov = TRUE,
                            norm_method_cov = "TMM",
                            class = NULL,
                            run_deg = TRUE,
                            BPPARAM = SerialParam())
 expect_named(tested$coef_data, c("(Intercept)", "cov"))
  tested <- run_tf_integration(expression = gene_exp_matrix,
                            tf_expression = gene_exp_matrix,
                            interactions = NULL,
                            type = "tf",
                            sequencing_data = TRUE,
                            species = "hsa",
                            normalize = TRUE,
                            norm_method = "TMM",
                            normalize_cov = TRUE,
                            norm_method_cov = "TMM",
                            class = class,
                            run_deg = TRUE,
                            BPPARAM = SerialParam())
  expect_type(tested, "list")
  expect_true("deg" %in% names(tested))
})
