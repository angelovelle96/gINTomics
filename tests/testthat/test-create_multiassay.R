test_that("multiassay generation works", {
  mmultiassay_ov2 <- create_multiassay(gene_exp = mmultiassay_ov@ExperimentList$gene_exp,
                                       methylation = mmultiassay_ov@ExperimentList$methylation,
                                       cnv_data = mmultiassay_ov@ExperimentList$cnv_data,
                                       miRNA_exp = mmultiassay_ov@ExperimentList$miRNA_exp,
                                       miRNA_cnv_data = mmultiassay_ov@ExperimentList$miRNA_cnv_data)

  expect_identical(mmultiassay_ov, mmultiassay_ov2)


})
