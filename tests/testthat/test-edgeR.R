library(BiocParallel)

test_that(".run_edgeR_integration works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tested <- .run_edgeR_integration(response_var,
                                   covariates,
                                   interactions = "auto",
                                   design_mat_allgene = NULL,
                                   offset_allgene = NULL,
                                   offset_singlegene = NULL,
                                   normalize=TRUE,
                                   norm_method = "TMM",
                                   steady_covariates = NULL,
                                   reference = NULL,
                                   BPPARAM = SerialParam())
  expect_type(tested, "list")
  expect_named(tested[1], "coef_data")
  expect_named(tested[2], "pval_data")
})

test_that(".allgene_edgeR_model works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  tested <- .allgene_edgeR_model(response_var,
                                 design_mat_allgene = NULL,
                                 offset_allgene = NULL,
                                 norm_method = "TMM")
  expect_type(tested, "list")
  expect_true(class(tested) == "DGEGLM")
  expect_false(is.null(tested$samples))
})

test_that(".singlegene_edgeR_model works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  fit_all <- .allgene_edgeR_model(response_var = response_var,
                                  design_mat_allgene = NULL,
                                  offset_allgene = NULL,
                                  norm_method = "TMM")
  design_mat_allgene <- NULL
  offset_allgene <- NULL
  reference <- NULL
  steady_covariates <- NULL
  normalize <- TRUE
  norm_method <- "TMM"

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

  tmp <- .covariates_check(response_var=response_var,
                           covariates=covariates,
                           interactions = interactions,
                           steady_covariates=steady_covariates,
                           reference=reference)
  if(is.null(tmp)) return(NULL)
  covariates <- tmp$covariates
  response_var <- tmp$response_var
  interactions <- tmp$interactions
  original_id <- tmp$original_id
  fformula <- .generate_formula(interactions = interactions)
  tested <- .singlegene_edgeR_model(response_var=response_var,
                                    covariates=covariates,
                                    offset_singlegene = NULL,
                                    fformula = fformula,
                                    fit_all = fit_all,
                                    SerialParam())
  expect_type(tested, "list")
  expect_true(length(tested) > 0)
})

test_that(".def_edger works", {

  offset_singlegene  <-  NULL
  design_mat_allgene <- NULL
  offset_allgene<- NULL
  offset_singlegene<- NULL
  normalize <- TRUE
  norm_method <-  "TMM"
  steady_covariates<- NULL
  reference<- NULL

  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
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
  tmp <- .covariates_check(response_var=response_var,
                           covariates=covariates,
                           interactions = interactions,
                           steady_covariates=steady_covariates,
                           reference=reference)
  if(is.null(tmp)) return(NULL)
  covariates <- tmp$covariates
  response_var <- tmp$response_var
  interactions <- tmp$interactions
  original_id <- tmp$original_id
  fformula <- .generate_formula(interactions = interactions)
  response_var <- t(response_var+1)
  fformula <- lapply(seq_along(fformula), function(x){
    list(formula=fformula[[x]], var_name=names(fformula)[x])
  })
  names(fformula) <- vapply(fformula, function(x) x$var_name, FUN.VALUE = "A")
  if(!is.null(offset_singlegene)) {
    if(is.list(offset_singlegene)) {
      if(length(offset_singlegene)!=nrow(response_var)) stop(str_wrap(
        "If provided as list, the length of offset_singlegene
            should be equal to the number of genes in response_var"))
    }
  }
  tested <-  bplapply(fformula, .def_edger,
                       response_var=response_var,
                       covariates=covariates,
                       fit_all=fit_all,
                       offset_singlegene=offset_singlegene,
                       BPPARAM = SerialParam())
  expect_type(tested, "list")
  expect_true(length(tested) > 0)
})


test_that(".edger_coef_test works", {

  interactions <- "auto"
  design_mat_allgene <- NULL
  offset_allgene <- NULL
  offset_singlegene <- NULL
  normalize <- TRUE
  norm_method <- "TMM"
  steady_covariates <-  NULL
  reference <-  NULL

  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)

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

  tmp <- .covariates_check(response_var=response_var,
                           covariates=covariates,
                           interactions = interactions,
                           steady_covariates=steady_covariates,
                           reference=reference)
  if(is.null(tmp)) return(NULL)
  covariates <- tmp$covariates
  response_var <- tmp$response_var
  interactions <- tmp$interactions
  original_id <- tmp$original_id
  fformula <- .generate_formula(interactions = interactions)
  fit_gene <- .singlegene_edgeR_model(
    response_var = response_var,
    covariates=covariates,
    fit_all = fit_all,
    fformula = fformula,
    offset_singlegene = offset_singlegene,
    BPPARAM = SerialParam()
  )
  tested <- .edger_coef_test(fit_list=fit_gene,
                             SerialParam())
  expect_type(tested, "list")
  expect_true(length(tested) > 0)
})

test_that(".def_coef_test works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)

  interactions <- "auto"
  design_mat_allgene <- NULL
  offset_allgene <- NULL
  offset_singlegene <- NULL
  normalize <- TRUE
  norm_method <- "TMM"
  steady_covariates <-  NULL
  reference <-  NULL

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

  tmp <- .covariates_check(response_var=response_var,
                           covariates=covariates,
                           interactions = interactions,
                           steady_covariates=steady_covariates,
                           reference=reference)
  if(is.null(tmp)) return(NULL)
  covariates <- tmp$covariates
  response_var <- tmp$response_var
  interactions <- tmp$interactions
  original_id <- tmp$original_id
  fformula <- .generate_formula(interactions = interactions)

  fit_gene <- .singlegene_edgeR_model(
    response_var = response_var,
    covariates=covariates,
    fit_all = fit_all,
    fformula = fformula,
    offset_singlegene = offset_singlegene,
    BPPARAM = SerialParam()
  )

  tested <- bplapply(fit_gene, .def_coef_test,
                       BPPARAM = SerialParam())
  expect_type(tested, "list")
  expect_true(length(tested) > 0)
})
