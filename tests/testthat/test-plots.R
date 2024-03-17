test_that("plot_network generates a visNetwork object", {
  data_table <- data_shiny_tests$data_table
  tested <- plot_network(data_table)
  expect_true(inherits(tested, "visNetwork"),
              "The result should be a visNetwork object")
})

library(shiny)
test_that("plot_venn generates a plotly object", {
   data_table <- data_shiny_tests$data_table
  tested <- plot_venn(data_table, class="A")
  expect_true(inherits(tested, "plotly"),
              "The result should be a plotly object")
})


test_that("plot_volcano generates a plotly object", {
   data_table <- data_shiny_tests$data_table
  plot <- plot_volcano(data_table, class = "A", omics = "gene_genomic_res", cnv_met = "cnv")
  expect_true(inherits(plot, "plotly"),
              "The result should be a ggplot object")
})


test_that("plot_ridge generates a ggplot object", {
  data_table <- data_shiny_tests$data_table
  plot <- plot_ridge(data_table, class = "A", omics = "omics1", cnv_met = "cnv")
  expect_true(inherits(plot, "ggplot"),
              "The result should be a ggplot object")
})


test_that("plot_heatmap generates a HeatmapList object", {
  data_table <- data_shiny_tests$data_table
  multiomics_integration <- data_shiny_tests$multiomics_integration
  plot <- plot_heatmap(multiomics_integration, data_table, omics = "gene_genomic_res", class="A")
   expect_true(inherits(plot, "HeatmapList"),
              "The result should be a HeatmapList object")
})


test_that("plot_chr_distribution generates a ggplot object", {
  data_table <- data_shiny_tests$data_table
  plot <- plot_chr_distribution(data_table, class="A", omics="gene_genomic_res")
  expect_true(inherits(plot, "ggplot"),
              "The result should be a ggplot object")
})


test_that("plot_tf_distribution generates a ggplot object", {
  data_table <- data_shiny_tests$data_table
  plot <- plot_tf_distribution(data_table, class="A")
  expect_true(inherits(plot, "ggplot"),
              "The result should be a ggplot object")
})


test_that("dot_plotly generates a plotly object", {
  enrich_result <- data.frame(
    Description = paste("Gene", 1:10),
    GeneRatio = runif(10),
    Count = sample(1:100, 10, replace = TRUE),
    pvalue = runif(10),
    qvalue = runif(10)
  )
  plot <- dot_plotly(enrich_result)
  expect_true(inherits(plot, "list"),
              "The result should be a list object")
})
