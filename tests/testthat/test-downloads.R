test_that(".download_mirna_target works", {

   data <- data_shiny_tests$multiassay
  tf_expression <- as.matrix(t(assay(data[["gene_exp"]])))
  species <- "hsa"
  tested <- .download_mirna_target(miRNAs = colnames(tf_expression),
                                             species = species)
  expect_type(tested, "list")
  expect_true(all(vapply(tested, function(x) length(unique(x)), FUN.VALUE = list()) == 1))
})

test_that(".download_tf_mirna works", {
  data <- data_shiny_tests$multiassay
  expression <- as.matrix(t(assay(data[["gene_exp"]])))
  species <- "hsa"
  tested <- .download_tf_mirna(miRNAs = colnames(expression),
                                         species = species)
  expect_type(tested, "list")
  expect_true(all(vapply(tested, function(x) length(unique(x)), FUN.VALUE = list()) == 1))
})

test_that(".download_tf works", {
  data <- data_shiny_tests$multiassay
  expression <- as.matrix(t(assay(data[["gene_exp"]])))
  species <- "hsa"
  tested <- .download_tf(genes = colnames(expression),
                                   species = species)
  expect_type(tested, "list")
})


test_that(".download_gene_info works", {
  tested1 <- .download_gene_info(genes="TP53", biomaRt = FALSE)
  expect_type(tested1, "list")
  expected_columns <- c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id",
                        "chromosome_name", "start_position", "end_position", "band")
  expect_true(all(expected_columns %in% colnames(tested1)))
  expect_gt(nrow(tested1), 0)
})
