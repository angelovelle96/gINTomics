

test_that(".gint_dashboardsidebar works", {
  out <- .gint_dashboardsidebar(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_tabitem_home works", {
  out <- .gint_tabitem_home(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_coefDistribGenomic works", {
  out <- .gint_subItem_coefDistribGenomic(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_HeatmapGenomic works", {
  out <- .gint_subItem_HeatmapGenomic(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_chrDistribGenomic works", {
  out <- .gint_subItem_chrDistribGenomic(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_tabitem_enr works", {
  out <- .gint_tabitem_enr(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_enrichTranscript works", {
  out <- .gint_subItem_enrichTranscript(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_coefDistribTranscript works", {
  out <- .gint_subItem_coefDistribTranscript(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_chrDistribTranscript works", {
  out <- .gint_subItem_chrDistribTranscript(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_networkTranscript works", {
  out <- .gint_subItem_networkTranscript(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_coefDistribDEGs works", {
  out <- .gint_subItem_coefDistribDEGs(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_HeatmapDEGs works", {
  out <- .gint_subItem_HeatmapDEGs(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_chrDistribDEGs works", {
  out <- .gint_subItem_chrDistribDEGs(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_networkDEGs works", {
  out <- .gint_subItem_networkDEGs(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_circosCompleteInt works", {
  out <- .gint_subItem_circosCompleteInt(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".gint_subItem_tableCompleteInt works", {
  out <- .gint_subItem_tableCompleteInt(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})

test_that(".create_ui works", {
  out <- .create_ui(data_shiny_tests$data_table)
  expect_s3_class(out, "shiny.tag")
})
