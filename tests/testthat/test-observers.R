
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




