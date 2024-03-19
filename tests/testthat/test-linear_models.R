library(BiocParallel)

test_that(".run_lm_integration works", {
  response_var <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$response_var
  covariates <- data_shiny_tests$multiomics_integration$gene_genomic_res$A$data$covariates
  tmp <- grep("_met$", colnames(covariates))
  covariates <- covariates[, -tmp]
  colnames(covariates) <- gsub("_cnv", "", colnames(covariates))
  interactions <- as.list(colnames(covariates))
  names(interactions) <- colnames(covariates)
  tested <- .run_lm_integration(response_var,
                                covariates,
                                interactions="auto",
                                steady_covariates=NULL,
                                reference=NULL,
                                normalize=TRUE,
                                norm_method="TMM",
                                BPPARAM=SerialParam())
  expect_type(tested, "list")
  expect_true(length(tested) > 0)
})


test_that(".def_lm works", {
  interactions <- "auto"
  steady_covariates <- NULL
  norm_method <- "TMM"
  normalize <- TRUE
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
  tested <- bplapply(fformula, .def_lm,
                                    data=data,
                                    BPPARAM = SerialParam())
  expect_type(tested, "list")
  expect_true(length(tested) > 0)
})
