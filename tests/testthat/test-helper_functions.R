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

test_that(".building_result_matrices works", {
  steady_covariates=NULL
  reference=NULL
  normalize=TRUE
  norm_method="TMM"
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tmp <- unlist(lapply(interactions, length))
  single_cov <- FALSE
  if(sum(tmp==1)==length(tmp)){
    single_cov <- TRUE
  }
  tmp <- .data_check(response_var = response_var,
                     interactions = interactions,
                     covariates = covariates)
  if(is.null(tmp)) return(NULL)
  response_var <- tmp$response_var
  interactions <- tmp$interactions
  covariates <- tmp$covariates

  if(normalize==TRUE) response_var <- .data_norm(response_var,
                                                 method = norm_method,
                                                 RNAseq = FALSE)
  tmp <- .covariates_check(response_var = response_var,
                           covariates = covariates,
                           interactions = interactions,
                           steady_covariates = steady_covariates,
                           linear = TRUE)
  if(is.null(tmp)) return(NULL)
  covariates <- tmp$covariates
  response_var <- tmp$response_var
  interactions <- tmp$interactions
  steady_covariates <- tmp$steady_covariates
  original_id <- tmp$original_id
  fformula <- .generate_formula(interactions = interactions,
                                linear=TRUE)
  data <- cbind(response_var, covariates)
  lm_results <-  bplapply(fformula, .def_lm,
                          data=data,
                          BPPARAM = bpparam())
  names(lm_results) <- names(interactions)
  tested <- .building_result_matrices(model_results=lm_results,
                                      type="lm",
                                      single_cov=FALSE)
  expect_type(tested, "list")
  expect_identical(colnames(tested$coef)[1], "(Intercept)")
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

test_that(".data_norm works", {
  data <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  tested <- .data_norm(data,
                       method="TMM",
                       RNAseq=TRUE)
  expect_type(tested, "double")
  expect_identical(nrow(data), nrow(tested))
})

test_that(".id_conversion works", {
  results <- list()
  results$residuals <- data.frame(row.names = "a", a=1)
  results$coef_data <- data.frame(row.names = "a", a=1)
  results$pval_data <- data.frame(row.names = "a", a=1)
  results$fdr_data <- data.frame(row.names = "a", a=1)
  results$data <- list()
  results$data$response_var <- data.frame(row.names = "a", a=1)
  results$data$covariates <- data.frame(row.names = "a", a=1)
  dictionary <- list()
  dictionary$tranformed <- "a"
  dictionary$original <- "b"
  tested <- .id_conversion(dictionary = dictionary,
                           results = results)
  expect_identical(rownames(tested$residuals), "b")
  expect_identical(colnames(tested$coef_data), "b")
  expect_identical(rownames(tested$fdr_data), "b")
  })

test_that(".rf_selection works", {
  data <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  fformula <- as.formula("TSPAN6~FGR+FUCA2+NIPAL3+CFTR+CYP51A1+BAD")
  tested <- .rf_selection(data = data,
                          formula = fformula)
  expect_s3_class(tested, "formula")
})

test_that(".search_gene works", {
  data <- data_shiny_tests$data_table
  genes <- "FGR"
  tested <- .search_gene(genes,
                        model_res=NULL,
                        data_frame=data)
  expect_s3_class(tested, "data.frame")
})

test_that(".find_deg works", {
  data <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  class <- c(rep("A", 5), rep("B", 5))
  names(class) <- rownames(data)
  tested <- .find_deg(eexpression = t(data),
                      class=class,
                      RNAseq=TRUE,
                      norm_method="TMM",
                      normalize=TRUE)
  expect_s3_class(tested, "data.frame")
  expect_identical(colnames(tested)[1], "logFC")
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
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions = "auto"
  design_mat_allgene = NULL
  offset_allgene = NULL
  offset_singlegene = NULL
  normalize=TRUE
  norm_method = "TMM"
  steady_covariates = NULL
  reference = NULL
  BPPARAM = SerialParam()
  tmp <- .data_check(response_var = response_var,
                     covariates = covariates,
                     interactions = interactions)
  if(is.null(tmp)) return(NULL)
  response_var <- tmp$response_var
  covariates <- tmp$covariates
  interactions <- tmp$interactions
  tmp <- unlist(lapply(interactions, length))
  single_cov <- FALSE
  if(sum(tmp==1)==length(tmp)){
    single_cov <- TRUE
  }

  if(normalize==FALSE) offset_allgene <- matrix(1, ncol(response_var),
                                                nrow(response_var))
  fit_all <- .allgene_edgeR_model(
    response_var = response_var,
    design_mat_allgene = design_mat_allgene,
    offset_allgene = offset_allgene,
    norm_method = norm_method
  )
  tested <- residuals(fit_all)
  expect_type(tested, "double")
  expect_identical(rownames(tested), colnames(response_var))
})
