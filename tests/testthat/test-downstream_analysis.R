


test_that("run_tf_enrich works", {
  data <- data_shiny_tests$data
  data <- data[data$omics=="tf_res",]
  tested <- run_tf_enrich(extracted_data=data,
                          species="hsa",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1,
                          pAdjustMethod="BH",
                          ont = "BP",
                          run_kegg=FALSE)
  expect_type(tested, "list")
})

test_that("run_genomic_enrich works", {
  expect_error(run_genomic_enrich(extracted_data=NULL,
                              species="hsa",
                              pvalueCutoff = 0.1,
                              qvalueCutoff = 0.1,
                              pAdjustMethod="none",
                              ont = "BP",
                              run_kegg=FALSE))
})


