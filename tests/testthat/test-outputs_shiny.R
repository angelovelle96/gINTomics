library(shiny)
library(shiny.gosling)
test_that(".prepare_network works", {
  data_table <- data_shiny_tests$data_table
  network <- .prepare_network(data_table)
  expect_type(network, "list")
  expected_components <- c("nodes", "edges", "legend_edges", "legend_nodes")
  expect_true(all(expected_components %in% names(network)))
  expect_true("id" %in% colnames(network$nodes))
  expect_true("label" %in% colnames(network$nodes))
  expect_true("shape" %in% colnames(network$nodes))
  expect_true("title" %in% colnames(network$nodes))
  expect_true("shadow" %in% colnames(network$nodes))
  expect_true("color" %in% colnames(network$nodes))
  expect_true("width" %in% colnames(network$nodes))
  expect_true("from" %in% colnames(network$edges))
  expect_true("to" %in% colnames(network$edges))
  expect_true("width" %in% colnames(network$edges))
  expect_true("color" %in% colnames(network$edges))
  expect_true("length" %in% colnames(network$edges))
  expect_true("title" %in% colnames(network$edges))
  expect_true("label" %in% colnames(network$legend_nodes))
  expect_true("shape" %in% colnames(network$legend_nodes))
  expect_true("color" %in% colnames(network$legend_nodes))
  expect_true("color" %in% colnames(network$legend_edges))
  expect_true("label" %in% colnames(network$legend_edges))
  expect_true("arrows" %in% colnames(network$legend_edges))
})

test_that(".build_network works", {
  input <- reactiveValues(layout=FALSE,
                          physics=FALSE,
                          numInteractions=200,
                          SignificativityCriteria="pval",
                          PvalRange=0.5,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  nnet <- .prepare_network(data_shiny_tests$data_table)
  reactive_network <- .select_network(data_table = data_shiny_tests$data_table,
                                      input = input,
                                      output = output,
                                      network_data = nnet,
                                      deg = FALSE)
  tmp <- isolate(reactive_network())
  network <- .build_network(nodes = tmp$nodes,
                            edges = tmp$edges,
                            legend_nodes = tmp$legend_nodes,
                            legend_edges = tmp$legend_edges)
  tmp$edges <- tmp$edges[!is.na(tmp$edges$from),]
  expect_type(network, "list")
  expect_equal(length(network$x$nodes$id), nrow(tmp$nodes))
  expect_equal(length(network$x$edges$from), nrow(tmp$edges))
})

test_that(".build_venn works", {
  input <- reactiveValues(SignificativityCriteria="pval",
                          PvalRange=0.5,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  data_table <- data_shiny_tests$data_table
  venn_data <- .prepare_reactive_venn(data_table,
                                      input = input,
                                      output = output)
  venn_data <- isolate(venn_data())
  venn <- .build_venn(venn_data)
  expect_type(venn, "list")
  expect_s3_class(venn, "gg")
  expect_equal(names(venn)[1], "data")
})

test_that(".render_venn works", {
  input <- reactiveValues(SignificativityCriteria="pval",
                          PvalRange=0.5,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  data_table <- data_shiny_tests$data_table
  venn_data <- .prepare_reactive_venn(data_table,
                                      input = input,
                                      output = output)
  venn_data <- isolate(venn_data())
  tested <- .render_venn(venn_data)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".render_venn_table works", {
  input <- reactiveValues(SignificativityCriteria="pval",
                          PvalRange=0.5,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  data_table <- data_shiny_tests$data_table
  venn_data <- .prepare_reactive_venn(data_table,
                                      input = input,
                                      output = output)
  venn_data <- isolate(venn_data())
  tested <- .render_venn_table(venn_data)
  expect_s3_class(tested, "shiny.render.function")
})


test_that(".build_volcano works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=0.5,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  reactive_volcano <- .prepare_reactive_volcano(data_shiny_tests$data_table,
                                                input = input,
                                                output = output,
                                                deg = FALSE)
  tmp <- isolate(reactive_volcano())
  tested <- .build_volcano(tmp)
  expect_type(tested, "list")
  expect_s3_class(tested, "plotly")
  expect_equal(names(tested)[1], "x")
})

test_that(".render_volcano works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=0.5,
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  reactive_volcano <- .prepare_reactive_volcano(data_shiny_tests$data_table,
                                                input = input,
                                                output = output,
                                                deg = FALSE)
  tmp <- isolate(reactive_volcano())
  tested <- .render_volcano(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".build_ridge works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,0.5),
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  reactive_ridge <- .prepare_reactive_ridge(data_shiny_tests$data_table,
                                            input = input,
                                            output = output,
                                            deg = FALSE)
  tmp <- isolate(reactive_ridge())
  tested <- .build_ridge(ridge_data = tmp$df, quantiles = tmp$quantiles)
  expect_type(tested, "list")
  expect_s3_class(tested, "ggplot")
  expect_equal(names(tested)[1], "data")
})

test_that(".render_ridge works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,0.5),
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  reactive_ridge <- .prepare_reactive_ridge(data_shiny_tests$data_table,
                                            input = input,
                                            output = output,
                                            deg = FALSE)
  tmp <- isolate(reactive_ridge())
  tested <- .render_ridge(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".build_histo works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,0.5),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All")
  output <- reactiveValues()
  reactive_histo <- .prepare_reactive_histo(data_shiny_tests$data_table,
                                            input = input,
                                            output = output,
                                            deg = FALSE)
  tmp <- isolate(reactive_histo())
  tested <- .build_histo(tmp)
  expect_type(tested, "list")
  expect_s3_class(tested, "ggplot")
  expect_equal(names(tested)[1], "data")
})

test_that(".render_histo works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,0.5),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All")
  output <- reactiveValues()
  reactive_histo <- .prepare_reactive_histo(data_shiny_tests$data_table,
                                            input = input,
                                            output = output,
                                            deg = FALSE)
  tmp <- isolate(reactive_histo())
  tested <- .render_histo(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".render_histo_table works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,0.5),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All")
  output <- reactiveValues()
  reactive_histo <- .prepare_reactive_histo(data_shiny_tests$data_table,
                                            input = input,
                                            output = output,
                                            deg = FALSE,
                                            table = TRUE)
  tmp <- isolate(reactive_histo())
  tested <- .render_histo_table(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".build_histo_TFbyChr works", {
  input <- reactiveValues(IntegrationSelect="tf_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,1),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All")
  output <- reactiveValues
  reactive_histo_tf_transcript <- .prepare_reactive_histo_tf(data_shiny_tests$data_table,
                                                             input = input,
                                                             output = output,
                                                             deg = FALSE)
  tmp <- isolate(reactive_histo_tf_transcript())
  tested <- .build_histo_TFbyChr(tmp)
  expect_type(tested, "list")
  expect_s3_class(tested, "plotly")
  expect_equal(names(tested)[1], "x")
})

test_that(".render_histo_TF works", {
  input <- reactiveValues(IntegrationSelect="tf_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,1),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All")
  output <- reactiveValues
  reactive_histo_tf_transcript <- .prepare_reactive_histo_tf(data_shiny_tests$data_table,
                                                             input = input,
                                                             output = output,
                                                             deg = FALSE)
  tmp <- isolate(reactive_histo_tf_transcript())
  tested <- .render_histo_TF(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".render_histo_tf_table works", {
  input <- reactiveValues(IntegrationSelect="tf_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,1),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All")
  output <- reactiveValues
  reactive_histo_tf_transcript <- .prepare_reactive_histo_tf(data_shiny_tests$data_table,
                                                             input = input,
                                                             output = output,
                                                             deg = FALSE)
  tmp <- isolate(reactive_histo_tf_transcript())
  tested <- .render_histo_tf_table(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".render_ridge_table works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,0.5),
                          FdrRange=0.05,
                          ClassSelect="A")
  output <- reactiveValues()
  reactive_ridge <- .prepare_reactive_ridge(data_shiny_tests$data_table,
                                            input = input,
                                            output = output,
                                            deg = FALSE,
                                            table = TRUE)
  tmp <- isolate(reactive_ridge())
  tested <- .render_ridge_table(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".build_table works", {
  tested <- .build_table(data_shiny_tests$data_table)
  expect_s3_class(tested, "datatables")
  expect_s3_class(tested$x$data, "data.frame")
})

test_that(".render_table works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,1),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All",
                          degSelect="All")
  output <- reactiveValues()
  reactive_table <- .prepare_reactive_table(data_shiny_tests$data_table,
                                            input=input,
                                            output=output)
  tmp <- isolate(reactive_table())
  tested <- .render_table(tmp)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".download_csv works", {
  input <- reactiveValues(IntegrationSelect="gene_genomic_res",
                          genomicTypeSelect="cnv",
                          SignificativityCriteria="pval",
                          PvalRange=c(0,1),
                          FdrRange=0.05,
                          ClassSelect="A",
                          ChrSelect="All",
                          degSelect="All")
  output <- reactiveValues()
  tested <- .download_csv(plotType = "table",
                         input=input,
                         output=output,
                         data_table=data_shiny_tests$data_table)
  expect_s3_class(tested, "shiny.render.function")
  tested <- .download_csv(plotType = "histo",
                          input=input,
                          output=output,
                          data_table=data_shiny_tests$data_table)
  expect_s3_class(tested, "shiny.render.function")
  tested <- .download_csv(plotType = "venn",
                          input=input,
                          output=output,
                          data_table=data_shiny_tests$data_table)
  expect_s3_class(tested, "shiny.render.function")
})

test_that(".circos_preprocess works", {
  tested <- .circos_preprocess(data_shiny_tests$data)
  expect_s4_class(tested[[1]], "GenomicRanges")
  expect_s4_class(tested[[2]], "GenomicRanges")
  expect_s4_class(tested[[2]]@seqnames, "Rle")
})

test_that(".circos_track_cols works", {
  tested <- .circos_track_cols(seq_len(100))
  expect_length(tested, 100)
})

test_that(".create_single_track works", {
  use_gosling()
  tmp <- .circos_preprocess(data_shiny_tests$data)
  gr <- tmp$df_cnv
  gr$cov_value2 <- as.character(gr$cov_value_original)
  ccol <- .circos_track_cols(vvalues = gr$cov_value2)
  tested <- .create_single_track(data=gr,
                                    dataValue='cov_value',
                                    x_axis="none",
                                    xe_axis="none",
                                    y_axis="none",
                                    colorField="cov_value2",
                                    colorDomain=unique(gr$cov_value2),
                                    colorRange=ccol,
                                    tooltipField1="cov_value",
                                    tooltipTitle="cnv",
                                    tooltipAlt1="CNV Value:",
                                    tooltipField2="gene",
                                    tooltipAlt2="Gene Name:",
                                    tooltipField3="class",
                                    tooltipAlt3="Class:",
                                    tooltipField4="cnv_met",
                                    tooltipAlt4="Integration Type:",
                                    tooltipField5="cov_value_original",
                                    tooltipAlt5="Original CNV Value:",
                                    legend=FALSE,
                                    colorType="nominal",
                                    title="CNV")
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "data")
})

test_that(".create_tracks works", {
  use_gosling()
  gr <- .circos_preprocess(data_shiny_tests$data)
  tested <- .create_tracks(data = data_shiny_tests$data, gr = gr)
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "track_cnv")
  expect_length(tested, 9)
})

test_that(".create_cnv_track works", {
  use_gosling()
  gr <- .circos_preprocess(data_shiny_tests$data)
  tested <- .create_cnv_track(gr = gr$df_cnv)
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "data")
})

test_that(".create_met_track works", {
  use_gosling()
  gr <- .circos_preprocess(data_shiny_tests$data)
  tested <- .create_met_track(gr = gr$df_met)
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "data")
})

test_that(".create_expr_track works", {
  use_gosling()
  gr <- .circos_preprocess(data_shiny_tests$data)
  tested <- .create_expr_track(gr = gr$df_cnv)
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "data")
})

test_that(".create_coef_track works", {
  use_gosling()
  gr <- .circos_preprocess(data_shiny_tests$data)
  tested <- .create_coef_track(gr = gr$df_met)
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "data")
})

test_that(".create_cyto_track works", {
  use_gosling()
  tested <- .create_cyto_track()
  expect_type(tested, "list")
  expect_identical(names(tested)[1], "id")
})

test_that(".create_composed_view works", {
  use_gosling()
  gr <- .circos_preprocess(data_shiny_tests$data)
  tracks <- .create_tracks(data = data_shiny_tests$data, gr = gr)
  tested <- .create_composed_view(tracks = tracks, width = 100, height = 100)
  expect_type(tested, "list")
  expect_type(tested[[1]], "list")
  expect_identical(names(tested)[1], "circos_genomic")
})

test_that(".prepare_gen_heatmap works", {
  ht_opt$message <- FALSE
  integrationSelect <- "gene_genomic_res"
  numTopCNV <- 20
  numTopMET <- 20
  numTopCNVonly <- 20
  numTopMETonly <- 20
  numTopMiCNV <- 20
  classSelect <- "A"
  significativityCriteria <- "pval"
  pvalRange <- 0.1
  fdrRange <- 0.1
  scale <- "row"
  numSamples <- 50
  df_heatmap <- data_shiny_tests$multiomics_integration[[integrationSelect]][[
    classSelect]]$data$response_var
  data_table <- data_shiny_tests$data_table[data_shiny_tests$data_table$omics == integrationSelect,]
  data_table <- data_table[data_table$class == classSelect,]
  df_heatmap_t <- t(as.matrix(df_heatmap))
  tested <- .prepare_gen_heatmap(data_table = data_table,
                              df_heatmap = df_heatmap,
                              df_heatmap_t = df_heatmap_t,
                              significativityCriteria=significativityCriteria,
                              pvalRange = pvalRange,
                              fdrRange = fdrRange,
                              numTopCNV = numTopCNV,
                              numTopMET = numTopMET,
                              scale = scale,
                              numSamples = numSamples)
  expect_s4_class(tested, "Heatmap")
})

test_that(".prepare_cnv_heatmap works", {
  ht_opt$message <- FALSE
  integrationSelect <- "gene_genomic_res"
  numTopCNV <- 20
  numTopMET <- 20
  numTopCNVonly <- 20
  numTopMETonly <- 20
  numTopMiCNV <- 20
  classSelect <- "A"
  significativityCriteria <- "pval"
  pvalRange <- 0.1
  fdrRange <- 0.1
  scale <- "row"
  numSamples <- 50
  df_heatmap <- data_shiny_tests$multiomics_integration[[integrationSelect]][[
    classSelect]]$data$response_var
  data_table <- data_shiny_tests$data_table[data_shiny_tests$data_table$omics == integrationSelect,]
  data_table <- data_table[data_table$class == classSelect,]
  df_heatmap_t <- t(as.matrix(df_heatmap))
  data_table <- data_table[data_table$cnv_met=="cnv",]
  data_table$omics <- "gene_cnv_res"
  tested <- .prepare_cnv_heatmap(data_table = data_table,
                                 df_heatmap = df_heatmap,
                                 df_heatmap_t = df_heatmap_t,
                                 significativityCriteria=significativityCriteria,
                                 pvalRange = pvalRange,
                                 fdrRange = fdrRange,
                                 numTopCNVonly = numTopCNVonly,
                                 scale = scale,
                                 numSamples = numSamples)
  expect_s4_class(tested, "Heatmap")
})

test_that(".prepare_met_heatmap works", {
  ht_opt$message <- FALSE
  integrationSelect <- "gene_genomic_res"
  numTopCNV <- 20
  numTopMET <- 20
  numTopCNVonly <- 20
  numTopMETonly <- 20
  numTopMiCNV <- 20
  classSelect <- "A"
  significativityCriteria <- "pval"
  pvalRange <- 0.1
  fdrRange <- 0.1
  scale <- "row"
  numSamples <- 50
  df_heatmap <- data_shiny_tests$multiomics_integration[[integrationSelect]][[
    classSelect]]$data$response_var
  data_table <- data_shiny_tests$data_table[data_shiny_tests$data_table$omics == integrationSelect,]
  data_table <- data_table[data_table$class == classSelect,]
  df_heatmap_t <- t(as.matrix(df_heatmap))
  data_table <- data_table[data_table$cnv_met=="met",]
  data_table$omics <- "gene_met_res"
  tested <- .prepare_met_heatmap(data_table = data_table,
                                 df_heatmap = df_heatmap,
                                 df_heatmap_t = df_heatmap_t,
                                 significativityCriteria=significativityCriteria,
                                 pvalRange = pvalRange,
                                 fdrRange = fdrRange,
                                 numTopMETonly = numTopMETonly,
                                 scale = scale,
                                 numSamples = numSamples)
  expect_s4_class(tested, "Heatmap")
})

test_that(".prepare_gen_heatmap works", {
  ht_opt$message <- FALSE
  integrationSelect <- "mirna_cnv_res"
  numTopCNV <- 20
  numTopMET <- 20
  numTopCNVonly <- 20
  numTopMETonly <- 20
  numTopMiCNV <- 20
  classSelect <- "A"
  significativityCriteria <- "pval"
  pvalRange <- 0.1
  fdrRange <- 0.1
  scale <- "row"
  numSamples <- 50
  df_heatmap <- data_shiny_tests$multiomics_integration[[integrationSelect]][[
    classSelect]]$data$response_var
  data_table <- data_shiny_tests$data_table[data_shiny_tests$data_table$omics == integrationSelect,]
  data_table <- data_table[data_table$class == classSelect,]
  df_heatmap_t <- t(as.matrix(df_heatmap))
  tested <- .prepare_mirna_heatmap(data_table = data_table,
                                   df_heatmap = df_heatmap,
                                   df_heatmap_t = df_heatmap_t,
                                   significativityCriteria=significativityCriteria,
                                   pvalRange = pvalRange,
                                   fdrRange = fdrRange,
                                   numTopMiCNV = numTopMiCNV,
                                   scale = scale,
                                   numSamples = numSamples)
  expect_s4_class(tested, "Heatmap")
})

test_that(".run_bg works", {
  tested <- .run_bg(FFUN = paste0,
                    args = list("test"),
                    input = list(),
                    output = list())
  expect_true(tested$is_alive())
  i=0
  while(tested$is_alive()){
    i=i+1
  }
  expect_identical(tested$get_result(), "test")

})

test_that("run_shiny works", {
  tested <- run_shiny(data_shiny_tests$multiomics_integration)
  expect_s3_class(tested, "shiny.appobj")
  expect_identical(names(tested)[1], "httpHandler")
})

#file.remove("tests/testthat/Rplots.pdf")

