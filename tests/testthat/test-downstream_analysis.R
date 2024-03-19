library(BiocParallel)

test_that(".def_enrich works", {
  data <- data_shiny_tests$data
  tested <- .def_enrich(data=data,
                        species="hsa",
                        pvalueCutoff=0.05,
                        pAdjustMethod="BH",
                        qvalueCutoff=0.05,
                        ont="BP",
                        run_go=TRUE,
                        run_kegg=TRUE,
                        run_reactome=FALSE)
  expect_type(tested, "list")
})



test_that("run_tf_enrich works", {
  model_results <- data_shiny_tests$multiomics_integration
  tested <- run_tf_enrich(model_results=model_results,
                          species="hsa",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1,
                          pAdjustMethod="BH",
                          ont = "all",
                          BPPARAM = SerialParam(),
                          extracted_data=NULL)
  expect_type(tested, "list")
})
