library(shiny)
test_that(".render_reactive_network works", {

  input <- reactiveValues(layout=FALSE,
                          physics=FALSE,
                          numInteractions=200,
                          SignificativityCriteria="pval",
                          PvalRange=0.05,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  nnet <- .prepare_network(data_shiny_tests$data_table)
  reactive_network <- .select_network(data_table = data_shiny_tests$data_table,
                                      input = input,
                                      output = output,
                                      network_data = nnet,
                                      deg = FALSE)

  tested <- .render_reactive_network(input = input,
                                     output = output,
                                     reactive_network = reactive_network)
  tested <- tested$run()
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".select_network works", {

  input <- reactiveValues(layout=FALSE,
                          physics=FALSE,
                          numInteractions=200,
                          SignificativityCriteria="pval",
                          PvalRange=1,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  nnet <- .prepare_network(data_shiny_tests$data_table)
  tested <- .select_network(data_table = data_shiny_tests$data_table,
                            input = input,
                            output = output,
                            network_data = nnet,
                            deg = FALSE)
  ans <- isolate(tested())
  expect_identical(names(ans), c("nodes", "edges",
                                 "legend_edges", "legend_nodes"))
})


test_that(".prepare_reactive_venn works", {

  input <- reactiveValues(SignificativityCriteria="pval",
                          PvalRange=c(0,0.1),
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  tested <- .prepare_reactive_venn(data_table = data_shiny_tests$data_table,
                                  input = input,
                                  output = output,
                                  deg = FALSE)
  ans <- isolate(tested())
  expect_type(ans, "list")
  expect_identical(names(ans), c("cnv_sign_genes", "met_sign_genes"))
})


test_that(".prepare_reactive_volcano works", {

  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=0.1,
                          FdrRange=0.05)
  output <- reactiveValues()
  tested <- .prepare_reactive_volcano(data_shiny_tests$data_table,
                                      input = input,
                                      output = output,
                                      deg = FALSE)
  ans <- isolate(tested())
  expect_s3_class(ans, "data.frame")
})

test_that(".prepare_reactive_heatmap works", {
  input <- reactiveValues('test-IntegrationSelect'="gene_genomic_res",
                          'test-numTopGenesHeatmapCNV'=10,
                          'test-numTopGenesHeatmapMET'=10,
                          'test-numTopGenesHeatmapCNVonly'=10,
                          'test-numTopGenesHeatmapMETonly'=10,
                          'test-numTopGenesHeatmapmirna_cnv'=10,
                          'test-scaleHeatmap'="row",
                          'test-SignificativityCriteria'="pval",
                          'test-PvalRange'=0.05,
                          'test-FdrRange'=0.05,
                          'test-ClassSelect'="A")
  output <- reactiveValues()
  tested <- .prepare_reactive_heatmap(data_table=data_shiny_tests$data_table,
              multiomics_integration = data_shiny_tests$multiomics_integration,
              input=input,
              output=output,
              session = NULL,
              deg = FALSE,
              ns = NS("test"))
  tested1 <- tested$.label
  expect_identical(grep("significativityCriteria", tested1), as.integer(1))
})


test_that(".prepare_reactive_ridge works", {

  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  tested <- .prepare_reactive_ridge(data_shiny_tests$data_table,
                                    input = input,
                                    output = output,
                                    deg = FALSE)
  ans <- isolate(tested())
  expect_s3_class(ans$df, "data.frame")
  expect_identical(colnames(ans$df)[1], "response")
})


test_that(".prepare_reactive_histo works", {

  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  tested <- .prepare_reactive_histo(data_shiny_tests$data_table,
                                    input = input,
                                    output = output,
                                    deg = FALSE)
  ans <- isolate(tested())
  expect_s3_class(ans, "data.frame")
  expect_identical(colnames(ans)[1], "response")
})


test_that(".prepare_reactive_histo_tf works", {

  input <- reactiveValues(IntegrationSelect="tf_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  tested <- .prepare_reactive_histo_tf(data_shiny_tests$data_table,
                                      input = input,
                                      output = output,
                                      deg = FALSE)
  ans <- isolate(tested())
  expect_s3_class(ans, "data.frame")
  expect_identical(colnames(ans)[1], "TF")
})


test_that(".prepare_reactive_table works", {

  input <- reactiveValues(IntegrationSelect="tf_res",
                          SignificativityCriteria="pval",
                          ChrSelect="All",
                          degSelect="All",
                          PvalRange=c(0,1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  tested <- .prepare_reactive_table(data_shiny_tests$data_table,
                                   input = input,
                                   output = output)
  ans <- isolate(tested())
  expect_s3_class(ans, "data.frame")
  expect_identical(colnames(ans)[1], "response")
})


test_that(".check_reactive_bg_enrich works", {

  input <- reactiveValues()
  output <- reactiveValues()
  bg_enr <- .run_bg(FFUN = print,
                    input = input,
                    output = output,
                    args = list("a"))
  tested <- .check_reactive_bg_enrich(bg_enrich = bg_enr,
                                     input = input,
                                     output = output,
                                     session = NULL)
  ans <- isolate(tested())
  expect_true(ans%in%c("Enrichment completed",
                       paste0("Enrichment running in background, ",
                              "this may take several minutes")))
})


test_that(".reactive_gen_enrich works", {

  input <- reactiveValues(ClassSelect="A",
                          DBSelectEnrich="go",
                          genomicTypeSelect="cnv")
  output <- reactiveValues()
  bg_enr <- .run_bg(FFUN = as.list,
                    input = input,
                    output = output,
                    args = list("a"))
  tested <- .reactive_gen_enrich(bg_enrich = bg_enr,
                                 input = input,
                                 output = output,
                                 session = NULL)
  ans <- isolate(tested())
  expect_null(ans)
})


test_that(".reactive_tf_enrich works", {

  input <- reactiveValues(ClassSelect="A",
                          DBSelectEnrich="go",
                          genomicTypeSelect="cnv")
  output <- reactiveValues()
  bg_enr <- .run_bg(FFUN = as.list,
                    input = input,
                    output = output,
                    args = list("a"))
  tested <- .reactive_tf_enrich(bg_enrich = bg_enr,
                               input = input,
                               output = output,
                               session = NULL)
  ans <- isolate(tested())
  expect_null(ans)
})


test_that(".prepare_reactive_circos works", {
  library(shiny.gosling)
  use_gosling()
  input <- reactiveValues(ClassSelect="A",
                          layout="circular",
                          circosType="Gene",
                          ChrSelect="All")
  output <- reactiveValues()
  tested <- .prepare_reactive_circos(data = data_shiny_tests$data,
                                    input = input,
                                    output = output)
  ans <- isolate(tested())
  expect_type(ans, "list")
  expect_null(ans$title)
})
