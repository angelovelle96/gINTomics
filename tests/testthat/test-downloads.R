test_that("Test .download_mirna_target function", {
  miRNAs <- c("mir-21", "mir-22")
  species <- "hsa"
  result <- .download_mirna_target(miRNAs, species)
  expect_type(result, "list")
  expect_true(all(sapply(result, function(x) length(unique(x))) == 1))
})

test_that("Test .download_tf_mirna function", {
  miRNAs <- c("mir-21", "mir-22")
  species <- "hsa"
  result <- .download_tf_mirna(miRNAs, species)
  expect_type(result, "list")
  expect_true(all(sapply(result, function(x) length(unique(x))) == 1))
})

test_that("Test .download_tf function", {
  genes <- c("egfr", "myc2", "mcga")
  species <- "hsa"
  result <- .download_tf(genes, species)
  expect_type(result, "list")
  expect_true(all(sapply(result, function(x) length(unique(x))) == 1))
})

test_that("Test .download_gene_info_biomart function", {
  genes <- c("ENSG00000157764", "ENSG00000169174", "ENSG00000115209")
  species <- "hsa"
  filters <- c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id")
  result <- .download_gene_info_biomart(genes, species, filters)
  expect_type(result, "list")
  expected_columns <- c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id",
                        "chromosome_name", "start_position", "end_position", "band")
  expect_true(all(expected_columns %in% colnames(result)))
  expect_gt(nrow(result), 0)
})

test_that("Test .download_gene_info_org function", {
  genes <- c("ENSG00000157764", "ENSG00000169174", "ENSG00000115209")
  species <- "hsa"
  result <- .download_gene_info_org(genes, species)
  expect_type(result, "list")
  expected_columns <- c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id",
                        "chromosome_name", "start_position", "end_position", "band")
  expect_true(all(expected_columns %in% colnames(result)))
  expect_gt(nrow(result), 0)
})

test_that("Test .download_gene_info function", {
  genes <- c("ENSG00000157764", "ENSG00000169174", "ENSG00000115209")
  species <- "hsa"
  result1 <- .download_gene_info(genes, species, biomaRt = FALSE)
  result2 <- .download_gene_info(genes, species, biomaRt = TRUE)
  expect_type(result1, "list")
  expect_type(result2, "list")
  expect_false(identical(result1, result2))
  expected_columns <- c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id",
                        "chromosome_name", "start_position", "end_position", "band")
  expect_true(all(expected_columns %in% colnames(result1)))
  expect_true(all(expected_columns %in% colnames(result2)))
  expect_gt(nrow(result1), 0)
  expect_gt(nrow(result2), 0)
})
