test_that("Test .prepare_network function", {
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

test_that("Test .build_network function", {
  nodes <- data.frame(id = c("TF1", "TF2", "Target1", "miRNA1"),
                      label = c("TF1", "TF2", "Target1", "miRNA1"),
                      shape = c("diamond", "diamond", "circle", "triangle"),
                      title = c("TF1", "TF2", "Target1", "miRNA1"),
                      shadow = c(TRUE, TRUE, TRUE, TRUE),
                      color = c("#FFBA01", "#FFBA01", "#CCE7C9", "#9969C7"),
                      width = c(10, 10, 10, 10),
                      stringsAsFactors = FALSE)
  edges <- data.frame(from = c("TF1", "TF2", "miRNA1"),
                      to = c("Target1", "Target1", "miRNA1"),
                      width = c(4, 6, 8),
                      color = c("#4169E1", "#4169E1", "#ED5564"),
                      length = c(500, 500, 500),
                      title = c("pval: 0.01, coef: 0.5, FDR: 0.02",
                                "pval: 0.05, coef: -0.7, FDR: 0.07",
                                "pval: 0.02, coef: 0.3, FDR: 0.03"),
                      stringsAsFactors = FALSE)
  legend_nodes <- data.frame(label = c("TF", "Target", "miRNA"),
                             shape = c("diamond", "circle", "triangle"),
                             color = c("#FFBA01", "#CCE7C9", "#9969C7"),
                             stringsAsFactors = FALSE)
  legend_edges <- data.frame(color = c("#4169E1", "#ED5564", "black"),
                             label = c("UPregulate", "DOWNregulate", "no-effect"),
                             arrows = c("to", "to", "to"),
                             stringsAsFactors = FALSE)
  network <- .build_network(nodes, edges, legend_nodes, legend_edges)
  expect_type(network, "list")
  expect_equal(length(network$x$data$nodes$id), nrow(nodes))
  expect_equal(length(network$x$data$edges$from), nrow(edges))
  expect_equal(length(network$x$data$legend$groups$groupname), nrow(legend_nodes))
  expect_equal(length(network$x$data$legend$edges$color), nrow(legend_edges))
})

test_that("Test .build_venn function", {
  data_table <- data_shiny_tests$data_table
  venn_data <- .prepare_reactive_venn(data_table)
  venn <- .build_venn(venn_data)
  expect_type(venn, "ggplot")
  expect_equal(length(unique(venn_data$cnv_sign_genes$cov)), length(venn_data$cnv_sign_genes$cov))
  expect_equal(length(unique(venn_data$met_sign_genes$cov)), length(venn_data$met_sign_genes$cov))
})

test_that(".render_venn function renders plotly plot", {
  # Dati di esempio
  venn_data <- list(
    cnv_sign_genes = data.frame(cov = c("Gene1", "Gene2", "Gene3", "Gene4")),
    met_sign_genes = data.frame(cov = c("Gene2", "Gene3", "Gene5", "Gene6"))
  )
  reactive_venn <- function() {
    return(venn_data)
  }

  # Chiamata alla funzione
  plot <- .render_venn(reactive_venn)

  # Verifica che l'output sia di tipo plotly
  expect_is(plot, "plotly::plotly")

  # Verifica che il diagramma sia stato costruito correttamente
  expect_true(!is.null(plot$x$data))
})

test_that(".render_venn function handles null data correctly", {
  # Funzione reattiva che restituisce dati nulli
  reactive_venn_null <- function() {
    return(NULL)
  }

  # Chiamata alla funzione con dati nulli
  plot_null <- .render_venn(reactive_venn_null)

  # Verifica che il plot restituito sia vuoto
  expect_true(is.null(plot_null$x$data))
})

test_that(".render_venn_table function renders data table", {
  # Dati di esempio
  venn_data <- list(
    cnv_sign_genes = data.frame(cov = c("Gene1", "Gene2", "Gene3", "Gene4")),
    met_sign_genes = data.frame(cov = c("Gene2", "Gene3", "Gene5", "Gene6"))
  )
  reactive_venn <- function() {
    return(venn_data)
  }

  # Chiamata alla funzione
  datatable <- .render_venn_table(reactive_venn)

  # Verifica che l'output sia di tipo DataTable
  expect_is(datatable, "data.frame")

  # Verifica che la tabella sia stata costruita correttamente
  expect_equal(nrow(datatable), 2) # Numero di righe atteso
  expect_equal(ncol(datatable), 1) # Numero di colonne atteso
  expect_equal(names(datatable), "Genes") # Nome della colonna atteso
  expect_equal(datatable$Genes, c("Gene2", "Gene3")) # Contenuto atteso
})

test_that(".render_venn_table function handles null data correctly", {
  # Funzione reattiva che restituisce dati nulli
  reactive_venn_null <- function() {
    return(NULL)
  }

  # Chiamata alla funzione con dati nulli
  datatable_null <- .render_venn_table(reactive_venn_null)

  # Verifica che la tabella restituita sia vuota
  expect_equal(nrow(datatable_null), 0)
  expect_equal(ncol(datatable_null), 0)
})

test_that(".build_volcano function generates volcano plot", {
  # Dati di esempio
  volcano_data <- data.frame(
    cov = c("Gene1", "Gene2", "Gene3", "Gene4"),
    coef = c(1.5, -2.3, 0.8, -1.9),
    pval_fdr = c(0.001, 0.0001, 0.01, 0.05),
    group = c("Group1", "Group2", "Group1", "Group2"),
    class = c("Class1", "Class2", "Class1", "Class2")
  )

  # Chiamata alla funzione
  plot <- .build_volcano(volcano_data)

  # Verifica che l'output sia di tipo Plotly
  expect_is(plot, "list")
  expect_is(plot$x$data[[1]], "list")

  # Verifica che il titolo del grafico sia corretto
  expect_equal(plot$x$layout$title$text, "Volcano Plot")

  # Verifica che i dati siano stati correttamente trasferiti al plot
  expect_equal(plot$x$data[[1]]$x, c(1.5, -2.3, 0.8, -1.9))
  expect_equal(plot$x$data[[1]]$y, c(-log10(0.001), -log10(0.0001), -log10(0.01), -log10(0.05)))
  expect_equal(plot$x$data[[1]]$marker$color, c("Group1", "Group2", "Group1", "Group2"))
  expect_equal(plot$x$data[[1]]$marker$symbol, c("Class1", "Class2", "Class1", "Class2"))
  expect_equal(plot$x$data[[1]]$text, c("Group: Group1 <br> Class: Class1 <br> Name: Gene1 <br> Pval/FDR(-log10): 3.0 <br> coef 1.5",
                                         "Group: Group2 <br> Class: Class2 <br> Name: Gene2 <br> Pval/FDR(-log10): 4.0 <br> coef -2.3",
                                         "Group: Group1 <br> Class: Class1 <br> Name: Gene3 <br> Pval/FDR(-log10): 2.0 <br> coef 0.8",
                                         "Group: Group2 <br> Class: Class2 <br> Name: Gene4 <br> Pval/FDR(-log10): 1.301 <br> coef -1.9"))
})

test_that(".render_volcano function renders volcano plot", {
  # Dati di esempio
  volcano_data <- data.frame(
    cov = c("Gene1", "Gene2", "Gene3", "Gene4"),
    coef = c(1.5, -2.3, 0.8, -1.9),
    pval_fdr = c(0.001, 0.0001, 0.01, 0.05),
    group = c("Group1", "Group2", "Group1", "Group2"),
    class = c("Class1", "Class2", "Class1", "Class2")
  )

  # Funzione reattiva di esempio
  reactive_volcano <- function() {
    return(volcano_data)
  }

  # Chiamata alla funzione
  plot <- .render_volcano(reactive_volcano)

  # Verifica che l'output sia di tipo Plotly
  expect_is(plot, "htmlwidget")
  expect_true(any(grepl("plotly", class(plot))))

  # Verifica che il grafico sia stato correttamente generato
  plotly_output <- plotly:::plotly_build(plot)
  expect_equal(plotly_output$data[[1]]$x, c(1.5, -2.3, 0.8, -1.9))
  expect_equal(plotly_output$data[[1]]$y, c(-log10(0.001), -log10(0.0001), -log10(0.01), -log10(0.05)))
  expect_equal(plotly_output$data[[1]]$marker$color, c("Group1", "Group2", "Group1", "Group2"))
  expect_equal(plotly_output$data[[1]]$marker$symbol, c("Class1", "Class2", "Class1", "Class2"))
  expect_equal(plotly_output$data[[1]]$text, c("Group: Group1 <br> Class: Class1 <br> Name: Gene1 <br> Pval/FDR(-log10): 3.0 <br> coef 1.5",
                                                "Group: Group2 <br> Class: Class2 <br> Name: Gene2 <br> Pval/FDR(-log10): 4.0 <br> coef -2.3",
                                                "Group: Group1 <br> Class: Class1 <br> Name: Gene3 <br> Pval/FDR(-log10): 2.0 <br> coef 0.8",
                                                "Group: Group2 <br> Class: Class2 <br> Name: Gene4 <br> Pval/FDR(-log10): 1.301 <br> coef -1.9"))
})

test_that(".build_ridge function generates ridgeline plot", {
  # Dati di esempio
  ridge_data <- data.frame(
    coef = rnorm(100),
    significance = runif(100)
  )
  quantiles <- c(-2, 2)  # Limiti dell'asse x

  # Chiamata alla funzione
  plot <- .build_ridge(ridge_data, quantiles)

  # Verifica che l'output sia di tipo ggplot
  expect_is(plot, "ggplot")

  # Verifica che il grafico abbia i titoli corretti
  expect_equal(ggplot2::ggtitle(plot), "Ridgeline Plot")
  expect_equal(names(ggplot2::xlab(plot)), "label")
  expect_equal(names(ggplot2::ylab(plot)), "label")
  expect_equal(ggplot2::xlab(plot)$label, "Value")
  expect_equal(ggplot2::ylab(plot)$label, "Significativity")

  # Verifica che i limiti dell'asse x siano stati impostati correttamente
  expect_equal(attr(ggplot2::layer_scales(plot)$x, "limits"), quantiles)

  # Verifica che il grafico includa una densità delle righe
  expect_true(any(grepl("geom_density_ridges", ggplot2::ggplot_build(plot)$layout$name)))

  # Verifica che il grafico abbia la giusta scala di colori
  expect_true(any(grepl("scale_fill_manual", ggplot2::ggplot_build(plot)$layout$name)))
})

test_that(".render_ridge function renders ridgeline plot", {
  # Dati di esempio
  ridge_data <- data.frame(
    coef = rnorm(100),
    significance = runif(100)
  )
  quantiles <- c(-2, 2)  # Limiti dell'asse x
  reactive_ridge <- function() {
    list(df = ridge_data, quantiles = quantiles)
  }

  # Chiamata alla funzione
  plot <- .render_ridge(reactive_ridge)

  # Verifica che l'output sia di tipo ggplot
  expect_is(plot, "ggplot")

  # Verifica che il grafico abbia i titoli corretti
  expect_equal(ggplot2::ggtitle(plot), "Ridgeline Plot")
  expect_equal(names(ggplot2::xlab(plot)), "label")
  expect_equal(names(ggplot2::ylab(plot)), "label")
  expect_equal(ggplot2::xlab(plot)$label, "Value")
  expect_equal(ggplot2::ylab(plot)$label, "Significativity")

  # Verifica che i limiti dell'asse x siano stati impostati correttamente
  expect_equal(attr(ggplot2::layer_scales(plot)$x, "limits"), quantiles)

  # Verifica che il grafico includa una densità delle righe
  expect_true(any(grepl("geom_density_ridges", ggplot2::ggplot_build(plot)$layout$name)))
})

test_that(".render_histo function renders histogram plot", {
  # Dati di esempio
  histo_data <- data.frame(
    chr_cov = sample(1:22, 100, replace = TRUE),
    significance = sample(c("significant", "not significant"), 100, replace = TRUE)
  )

  reactive_histo <- function() {
    histo_data
  }

  # Chiamata alla funzione
  plot <- .render_histo(reactive_histo)

  # Verifica che l'output sia di tipo ggplot
  expect_is(plot, "ggplot")

  # Verifica che il grafico abbia i titoli corretti
  expect_equal(ggplot2::ggtitle(plot), "Number of Genes with Significant Coefficients by Chromosome")
  expect_equal(names(ggplot2::xlab(plot)), "label")
  expect_equal(names(ggplot2::ylab(plot)), "label")
  expect_equal(ggplot2::xlab(plot)$label, "Chromosome")
  expect_equal(ggplot2::ylab(plot)$label, "Count")

  # Verifica che il grafico includa un layer di barre
  expect_true(any(grepl("geom_bar", ggplot2::ggplot_build(plot)$layout$name)))
})

test_that(".render_histo function renders histogram plot", {
  # Dati di esempio
  histo_data <- data.frame(
    chr_cov = sample(1:22, 100, replace = TRUE),
    significance = sample(c("significant", "not significant"), 100, replace = TRUE)
  )

  reactive_histo <- function() {
    histo_data
  }

  # Chiamata alla funzione
  plot <- .render_histo(reactive_histo)

  # Verifica che l'output sia di tipo plotly
  expect_is(plot, "plotly")

  # Verifica che il grafico abbia il titolo corretto
  expect_equal(plot$x$title$text, "Chromosome")
  expect_equal(plot$y$title$text, "Count")

  # Verifica che il grafico includa un layer di barre
  expect_true("type" %in% plot$x$data[[1]] && plot$x$data[[1]]$type == "bar")
})

test_that(".render_histo_table function renders data table", {
  # Dati di esempio
  table_data <- data.frame(
    chr_cov = sample(1:22, 100, replace = TRUE),
    significance = sample(c("significant", "not significant"), 100, replace = TRUE)
  )

  reactive_histo_table <- function() {
    table_data
  }

  # Chiamata alla funzione
  table <- .render_histo_table(reactive_histo_table)

  # Verifica che l'output sia di tipo DataTable
  expect_is(table, "dataTable")

  # Verifica che la tabella contenga i dati corretti
  expect_equal(nrow(table), nrow(table_data))
  expect_equal(ncol(table), ncol(table_data))
  expect_equal(colnames(table), colnames(table_data))
  expect_equal(rownames(table), rownames(table_data))
})

test_that(".build_histo_TFbyChr function creates correct plot", {
  # Dati di esempio
  histo_data <- data.frame(
    TF = c("TF1", "TF2", "TF3", "TF1", "TF2", "TF3"),
    Chromosome = c("chr1", "chr2", "chr1", "chr2", "chr1", "chr2"),
    Count = c(10, 15, 20, 25, 30, 35)
  )

  # Chiamata alla funzione
  plot <- .build_histo_TFbyChr(histo_data)

  # Verifica che l'output sia un oggetto plotly
  expect_is(plot, "list")
  expect_is(plot$x$data[[1]]$x, "formula")
  expect_is(plot$x$data[[1]]$y, "formula")
  expect_is(plot$x$data[[1]]$type, "bar")

  # Verifica che il layout sia stato impostato correttamente
  expect_equal(plot$x$layout$title$text, "Number of genes targeted by TFs/miRNAs")
  expect_equal(plot$x$layout$xaxis$title$text, "TF/miRNAs")
  expect_equal(plot$x$layout$yaxis$title$text, "Number of targets")
  expect_equal(plot$x$layout$barmode, "group")
})

