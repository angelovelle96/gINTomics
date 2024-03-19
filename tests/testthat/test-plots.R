library(shiny)

test_that("plot_network works", {
  data_table <- data_shiny_tests$data_table
  tested <- plot_network(data_table)
  expect_true(inherits(tested, "visNetwork"),
              "The result should be a visNetwork object")
})

test_that("plot_venn works", {
   data_table <- data_shiny_tests$data_table
  tested <- plot_venn(data_table, class="A")
  expect_true(inherits(tested, "plotly"),
              "The result should be a plotly object")
})


test_that("plot_volcano works", {
   data_table <- data_shiny_tests$data_table
  plot <- plot_volcano(data_table, class = "A", omics = "gene_genomic_res", cnv_met = "cnv")
  expect_true(inherits(plot, "plotly"),
              "The result should be a ggplot object")
})


test_that("plot_ridge works", {
  data_table <- data_shiny_tests$data_table
  plot <- plot_ridge(data_table, class = "A", omics = "omics1", cnv_met = "cnv")
  expect_true(inherits(plot, "ggplot"),
              "The result should be a ggplot object")
})


test_that("plot_heatmap works", {
  data_table <- data_shiny_tests$data_table
  multiomics_integration <- data_shiny_tests$multiomics_integration
  plot <- plot_heatmap(multiomics_integration, data_table, omics = "gene_genomic_res", class="A")
   expect_true(inherits(plot, "Heatmap"),
              "The result should be a HeatmapList object")
})


test_that("plot_chr_distribution works", {
  data_table <- data_shiny_tests$data_table
  plot <- plot_chr_distribution(data_table, class="A", omics="gene_genomic_res")
  expect_true(inherits(plot, "ggplot"),
              "The result should be a ggplot object")
})


test_that("plot_tf_distribution works", { ######3
  data_table <- data_shiny_tests$data_table
   class <- "A"
  plot <- plot_tf_distribution(data_table=data_table, class=class, pval=0.05)
  expect_true(inherits(plot, "ggplot"),
              "The result should be a ggplot object")
})


test_that("dot_plotly works", {  ########33


  plot <- dot_plotly(enrich_result=enrich_result,
                     title=NULL,
                     showCategory=10,
                     width=800,
                     height=700)
  expect_true(inherits(plot, "list"),
              "The result should be a list object")
})
