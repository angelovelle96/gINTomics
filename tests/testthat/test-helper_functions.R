library(BiocParallel)

test_that(".data_check works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tested <- .data_check(response_var=response_var,
                        covariates=covariates,
                        interactions=interactions)

  expect_type(tested$covariates, "list")
  expect_named(tested[1], "response_var")
})

test_that(".covariates_check works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tested <- .covariates_check(response_var=response_var,
                              covariates=covariates,
                              interactions=interactions,
                              steady_covariates=NULL,
                              linear=FALSE,
                              reference=NULL)
  expect_type(tested, "list")
  expect_named(tested[1], "covariates")
})

test_that(".generate_formula works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tested <- .generate_formula(interactions=interactions,
                              linear=FALSE)
  expect_type(tested, "list")

})

test_that(".building_result_matrices works", { ##########
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tmp <- unlist(lapply(interactions, length))
  response_var <- .data_norm(response_var,
                             method = norm_method,
                             RNAseq = FALSE)
  tmp <- .covariates_check(response_var = response_var,
                           covariates = covariates,
                           interactions = interactions,
                           steady_covariates = NULL,
                           linear = TRUE)
  tested <- .building_result_matrices(model_results=model_results,
                                      type="lm",
                                      single_cov=FALSE)
  fformula <- .generate_formula(interactions = interactions,
                                linear=TRUE)
  data <- cbind(response_var, covariates)
  lm_results <-  bplapply(fformula, .def_lm,
                          data=data,
                          BPPARAM = SerialParam())
  names(lm_results) <- names(interactions)

  tested <- .building_result_matrices(model_results = lm_results,
                                             type = "lm",
                                             single_cov = TRUE)
  expect_type(tested, "list")
})

test_that("create_multiassay works", {
  methylation <-as.matrix(assay(data_shiny_tests$multiassay[["methylation"]]))
  gene_exp_matrix <- as.matrix(assay(data_shiny_tests$multiassay[["gene_exp"]]))
  miRNA_exp_matrix <- as.matrix(assay(data_shiny_tests$multiassay[["miRNA_exp"]]))
  gene_cnv_matrix <- as.matrix(assay(data_shiny_tests$multiassay[["cnv_data"]]))
  miRNA_cnv_matrix <- as.matrix(assay(data_shiny_tests$multiassay[["miRNA_cnv_data"]]))
  tested <- create_multiassay(methylation=methylation,
                              cnv_data=gene_cnv_matrix,
                              gene_exp=gene_exp_matrix,
                              miRNA_exp=miRNA_exp_matrix,
                              miRNA_cnv_data=miRNA_cnv_matrix)
  expect_type(tested, "S4")
})

test_that(".data_norm works", {  ###########
  data <- data_shiny_tests$
  tested <- .data_norm(data,
                       method="TMM",
                       RNAseq=TRUE)
})

test_that(".id_conversion works", { #############
  tested <- .id_conversion(dictionary,
                           results)
})

test_that(".rf_selection works", { ##############
  data <- data_shiny_tests$data
  tested <- .rf_selection(data,
                          formula)
})

test_that("search_gene works", { ##############
  tested <- search_gene(genes,
                        model_res=NULL,
                        data_frame=NULL)
})

test_that(".find_deg works", {   #############
  tested <- .find_deg(eexpression,
                      class,
                      RNAseq=TRUE,
                      norm_method="TMM",
                      normalize=TRUE)
})

test_that(".shiny_preprocess works", {
  data <- data_shiny_tests$data
  tested <- .shiny_preprocess(data)
  expect_type(tested, "list")
})

test_that(".change_int_names works", {
  nnames <- c("gene_genomic_res",
              "gene_cnv_res",
              "gene_met_res",
              "mirna_cnv_res",
              "tf_res",
              "tf_mirna_res",
              "mirna_target_res")
  tested <- .change_int_names(nnames)
  expect_type(tested, "character")
})

test_that("residuals.DGEGLM works" , {  ###########
  tested <- residuals.DGEGLM(object, type=c("deviance", "pearson"))
})
