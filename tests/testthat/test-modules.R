library(shiny)

test_that(".server_histo works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_histo <- .server_histo(data_table=data_table,
                                input=input,
                                output=output,
                                session=session)
  expect_s3_class(server_histo, "shiny.render.function")
})


test_that(".server_histo2 works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_histo <- .server_histo2(data_table=data_table,
                                input=input,
                                output=output,
                                session=session)
  expect_s3_class(server_histo, "shiny.render.function")
})

test_that(".server_network works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(layout=FALSE,
                          physics=FALSE,
                          numInteractions=200,
                          SignificativityCriteria="pval",
                          PvalRange=0.05,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_network <- .server_network(data_table=data_table,
                                input=input,
                                output=output,
                                session=session)
  tested <- server_network$run()
  expect_s3_class(tested, "shiny.render.function")
})


test_that(".server_venn works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(SignificativityCriteria="pval",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_venn <- .server_venn(data_table=data_table,
                                 input=input,
                                 output=output,
                                 session=session)
  expect_s3_class(server_venn, "shiny.render.function")
})

test_that(".server_volcano works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_volcano <- .server_volcano(data_table=data_table,
                                 input=input,
                                 output=output,
                                 session=session)
  expect_s3_class(server_volcano, "shiny.render.function")
})

test_that(".server_ridge works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_ridge <- .server_ridge(data_table=data_table,
                                 input=input,
                                 output=output,
                                 session=session)
  expect_s3_class(server_ridge, "shiny.render.function")
})

test_that(".server_table works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          ChrSelect="All",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_table <- .server_table(data_table=data_table,
                                 input=input,
                                 output=output,
                                 session=session)
  expect_s3_class(server_table, "shiny.render.function")
})

test_that(".server_enrich_bg works", {
  data_table <- data_shiny_tests$data_table
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          SignificativityCriteria="pval",
                          genomicTypeSelect="cnv",
                          PvalRange=c(0,0.1),
                          FdrRange=c(0,0.05),
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  data_tf_enrich <- data_table[data_table$omics=="tf_res",]
  tested <- .server_enrich_bg(input=input,
                               output=output,
                               session=session,
                               extracted_data = NULL)
  expect_null(tested)

})

test_that(".server_circos works", {
  data <- data_shiny_tests$data
  input <- reactiveValues(circosType="gene",
                          ChrSelect="All",
                          layout="circular",
                          ClassSelect="A")
  output <- reactiveValues()
  session <- list()
  server_circos <- .server_circos(input=input,
                                  output=output,
                                  session=session,
                                  data=data)

  tested <- server_circos$.label
  expect_identical(grep("export_png", tested), as.integer(1))
})

