



.gint_dashboardsidebar <- function(data_table){
  myImgResources <- "imgResources/logo_gINTomics.png"
  addResourcePath(prefix = "imgResources", directoryPath = "inst/www/")
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Data Table", tabName = "page_table", icon = icon("table")),
      menuItem("Histo Plot", tabName = "page_histo"),
      menuItem("Ridgeline Plot", tabName = "page_ridge"),
      menuItem("Venn Diagram", tabName = "page_venn"),
      menuItem("Heatmap", tabName = "page_heatmap"),
      menuItem("Volcano Plot", tabName = "page_volcano"),
      menuItem("Circos Plot", tabName = "page_circos"),
      menuItem("Network", tabName = "page_network"),
      menuItem("Enrichment", tabName = "page_enrichment",
               menuSubItem("Gene Ontology",
                           tabName = "GO"),
               menuSubItem("Gene set",
                           tabName = "gene_set"),
               menuSubItem("Pathway",
                           tabName = "path")),
      menuItem("Class Comparison", tabName = "page_compare_class")
    )
    ,img(src = myImgResources, height = 100, width = 100, align="center"),
    tags$style(".left-side, .main-sidebar {padding-top: 80px}"))
}

#############################################################
##############################################################

.gint_tabitem_home <- function(data_table){
  myImgResources <- "imgResources/logo_gINTomics.png"
  addResourcePath(prefix = "imgResources", directoryPath = "inst/www/")
  tabItem(tabName = "home",
          fluidRow(box(title = "Welcome to gINTomics Interactive Visualizer",
            "This is the home page content."),
            column(width = 1, img(src = myImgResources,
                                  height = 100, width = 100))))
}

###################################################################
####################################################################

.gint_tabitem_page_table <- function(data_table){
  tabItem(tabName = "page_table",
          fluidRow(
            box(title = "Page containing data tables.",
                "1. ."),
            mainPanel(

              sidebarLayout(
                sidebarPanel(
                  # Input
                  selectInput(inputId = 'integrationSelectTable',
                              label = 'Integration Type:',
                              choices = unique(data_table$omics)),
                  selectInput(inputId = 'classSelectTable',
                              label = 'Class:',
                              choices = unique(data_table$class)),
                  selectInput(inputId = 'chrSelectTable',
                              label = 'Chr:',
                              choices = unique(data_table$chr_cov)),
                  selectInput(inputId = 'degSelectTable',
                              label = 'DEGs:',
                              choices = c('All', 'Only DEGs')),
                  selectInput(inputId = 'significativityCriteriaTable',
                              label = 'Significativity criteria:',
                              choices = c('pval', 'FDR')),
                  conditionalPanel(
                    condition = "input.significativityCriteriaTable == 'pval'",
                    sliderInput("pvalRangeTable",
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005)),
                  conditionalPanel(
                    condition = "input.significativityCriteriaTable == 'FDR'",
                    sliderInput("FDRRangeTable",
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005))
                ),
                fluidRow(
                  column(
                    width = 12,
                    offset = 1,
                    tags$div(
                      style = 'overflow-x: auto;',
                      DT::dataTableOutput('res_table')
                    )
                  )
                )


              )
            )
          )
  )
}

####################################################
#####################################################


.gint_tabitem_page_histo <- function(data_table){
  tabItem(tabName = "page_histo",
      fluidRow(
        box(title = "Page containing data table",
            "This is the content of Page 1."),
        mainPanel(
          tabsetPanel(type = 'tabs',
            tabPanel('XXX',
             sidebarLayout(
               sidebarPanel(
                 #input
                 selectInput(inputId = 'integrationSelectHisto',
                             label = 'Integration Type:',
                             choices = unique(data_table$omics)),
                 selectInput(inputId = 'classSelectHisto',
                             label = 'Class:',
                             choices = unique(data_table$class)),
                 selectInput(inputId = 'chrSelectHisto',
                             label = 'Chr:',
                             choices = c('All', unique(data_table$chr_cov))),
                 selectInput(inputId = 'degSelectHisto',
                             label = 'DEGs:',
                             choices = c('All', 'Only DEGs')),
                 selectInput(inputId = 'significativityCriteriaHisto',
                             label = 'Significativity criteria:',
                             choices = c('pval', 'FDR')),
                 conditionalPanel(
                   condition = "input.significativityCriteriaHisto == 'pval'",
                   sliderInput("pvalRangeHisto",
                               "P-Value Range:",
                               min = 0,
                               max = 1,
                               value = c(0, 0.05),
                               step = 0.005)),
                 conditionalPanel(
                   condition = "input.significativityCriteriaHisto == 'FDR'",
                   sliderInput("FDRRangeHisto",
                               "FDR-Value Range:",
                               min = 0,
                               max = 1,
                               value = c(0, 0.05),
                               step = 0.005))
               ),
               mainPanel(
                 plotlyOutput('histogramPlot')
               )
             )
    ),
    tabPanel('TFs',
             sidebarLayout(
               sidebarPanel(
                 #input

                 selectInput(inputId = 'classSelectHistoTFs',
                             label = 'Class:',
                             choices = unique(data_table$class)),
                 selectInput(inputId = 'degSelectHistoTFs',
                             label = 'DEGs:',
                             choices = c('All', 'Only DEGs')),
                 selectInput(inputId = 'significativityCriteriaHistoTFs',
                             label = 'Significativity criteria:',
                             choices = c('pval', 'FDR'),),
                 conditionalPanel(
                   condition = "input.significativityCriteriaHistoTFs == 'pval'",
                   sliderInput("pvalRangeHistoTFs",
                               "P-Value Range:",
                               min = 0,
                               max = 1,
                               value = c(0, 0.10),
                               step = 0.005)),
                 conditionalPanel(
                   condition = "input.significativityCriteriaHistoTFs == 'FDR'",
                   sliderInput("FDRRangeHistoTFs",
                               "FDR-Value Range:",
                               min = 0,
                               max = 1,
                               value = c(0, 0.10),
                               step = 0.005))
               ),
               mainPanel(
                 plotlyOutput('histogramPlotTFs'),
                 plotlyOutput('histogramPlotTFsByChromosome')
               )
             ))
            )
           )
          )
         )
}

##########################################################
###########################################################

.gint_tabitem_page_ridge <- function(data_table){
  tabItem(tabName = "page_ridge",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
              sidebarLayout(
                sidebarPanel(
                  selectInput(inputId = 'integrationSelectRidge',
                              label = 'Integration Type:',
                              choices = unique(data_table$omics)),
                  selectInput("classSelectRidge",
                              "Select Class:",
                              choices = unique(data_table$class)),
                  selectInput(inputId = 'degSelectRidge',
                              label = 'DEGs:',
                              choices = c('All','Only DEGs')),
                  selectInput(inputId = 'significativityCriteriaRidge',
                              label = 'Significativity criteria:',
                              choices = c('pval', 'FDR')),
                  conditionalPanel(
                    condition = "input.significativityCriteriaRidge == 'pval'",
                    sliderInput("pvalRangeRidge",
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005)),
                  conditionalPanel(
                    condition = "input.significativityCriteriaRidge == 'FDR'",
                    sliderInput("FDRRangeRidge",
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005))
                ),
                mainPanel(
                  # OUTPUT
                  plotOutput("ridgelinePlot")
                )
              )
            )
          )
  )
}
#############################################################
##############################################################

.gint_tabitem_page_venn <- function(data_table){
  tabItem(tabName = "page_venn",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(

              sidebarLayout(
                sidebarPanel(
                  #input
                  selectInput("integrationSelectVenn",
                              "Select Integration:",
                              choices = unique(data_table$omics)),
                  selectInput("classSelectVenn",
                              "Select Class:",
                              choices = unique(data_table$class)),
                  selectInput(inputId = 'degSelectVenn',
                              label = 'DEGs:',
                              choices = c('All', 'Only DEGs')),
                  selectInput(inputId = 'significativityCriteriaVenn',
                              label = 'Significativity criteria:',
                              choices = c('pval', 'FDR')),
                  conditionalPanel(
                    condition = "input.significativityCriteriaVenn == 'pval'",
                    sliderInput("pvalRangeVenn",
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005)),
                  conditionalPanel(
                    condition = "input.significativityCriteriaVenn == 'FDR'",
                    sliderInput("FDRRangeVenn",
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005))

                ),
                mainPanel(
                  # OUTPUT
                  plotlyOutput("venn_plot"),
                  dataTableOutput("common_genes_table")
                )

              )

            )
          )
  )
}
##################################################################
###################################################################

.gint_tabitem_page_heatmap <- function(data_table){
  tabItem(tabName = "page_heatmap",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(

              sidebarLayout(
                sidebarPanel(
                  #input
                  selectInput(inputId = 'integrationSelectHeatmap',
                              label = 'Integration Type:',
                              choices = unique(data_table$omics)),
                  numericInput("numTopGenesHeatmap",
                               "Number of top genes:",
                               value = 50,
                               min = 1,
                               max = 200),

                  selectInput("selectClassHeatmap",
                              "Select the class:",
                              choices = unique(data_table$class),
                              multiple = FALSE),
                  selectInput(inputId = 'degSelectHeatmap',
                              label = 'DEGs:',
                              choices = c('All','Only DEGs')),
                  selectInput(inputId = 'significativityCriteriaHeatmap',
                              label = 'Significativity criteria:',
                              choices = c('pval', 'FDR')),
                  conditionalPanel(
                    condition = "input.significativityCriteriaHeatmap == 'pval'",
                    sliderInput("pvalRangeHeatmap",
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = 0.05,
                                step = 0.005)),
                  conditionalPanel(
                    condition = "input.significativityCriteriaHeatmap == 'FDR'",
                    sliderInput("FDRRangeHeatmap",
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = 0.05,
                                step = 0.005))

                ),
                mainPanel(
                  # OUTPUT
                  InteractiveComplexHeatmapOutput('heatmap')
                )

              )
            )
          )
  )
}
################################################################
#################################################################

.gint_tabitem_page_volcano <- function(data_table){
  tabItem(tabName = "page_volcano",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
              tabsetPanel(
                type = 'tabs',
                tabPanel('XXXX',
              sidebarLayout(
                sidebarPanel(
                  #input
                  selectInput(inputId = 'integrationSelectVolcano',
                              label = 'Integration Type:',
                              choices = unique(data_table$omics)),
                  selectInput(inputId = 'genomicTypeSelectVolcano',
                              label = 'Integration Type:',
                              choices = unique(data_table$cnv_met)),
                  numericInput("num_top_genes",
                               "Number of top genes to visualize:",
                               value = 50),
                  sliderInput("pvalRangeVolcano",
                              "P-Value Range:",
                              min = 0,
                              max = 1,
                              value = c(0, 0.05),
                              step = 0.005)
                ),
                mainPanel(
                  # OUTPUT
                  plotlyOutput('volcanoPlot')
                )
              )
            )
          , tabPanel('DEGs'))
        )
      )
  )
}
########################################################################
#########################################################################

.gint_tabitem_page_network <- function(data_table){
  tabItem(tabName = "page_network",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
              tabsetPanel(
                type = 'tabs',
                tabPanel('XXXX',
              # OUTPUT
              checkboxInput("layoutNetwork",
                           label = "Layout:",
                           value = FALSE),
              checkboxInput("deg",
                            label = "deg:",
                            value = FALSE),
              visNetworkOutput("networkPlot",
                               height = 800,
                               width = 1600)
            )
          , tabPanel('DEGs'))
        )
      )
  )
}
##################################################################
###################################################################

.gint_tabitem_page_circos <- function(data_table){
  tabItem(tabName = "page_circos",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
                do.call("tabsetPanel", lapply(as.list(unique(data_table$omics)), function(x) tabPanel(x)))
        )
      )
  )
}
################################################################
#################################################################

.create_ui <- function(data_table){
  myImgResources <- "imgResources/logo_gINTomics.png"
  addResourcePath(prefix = "imgResources", directoryPath = "inst/www/")
  dashboardPage(
    dashboardHeader(title = span("gINTomics",
                                 span("Visualizer 1.0",
                                      style = "color: gray; font-size: 16px")),
                    tags$li(a(href="https://github.com/angelovelle96/gINTomics",
                              img(src = myImgResources,
                                  height = 80,
                                  width = 80),
                              style = "padding-top:1px; padding-bottom:1px;"),
                            class = "dropdown"),
                    tags$li(class = "dropdown",
                            tags$style(".main-header {max-height: 80px}"),
                            tags$style(".main-header .logo {height: 80px;}"),
                            tags$style(".sidebar-toggle {height: 80px;}"))
                            #tags$style(".navbar {min-height:80px !important}")
                     ),
    .gint_dashboardsidebar(),
    dashboardBody(
      tabItems(
        .gint_tabitem_home(data_table),
        .gint_tabitem_page_table(data_table),
        .gint_tabitem_page_histo(data_table),
        .gint_tabitem_page_ridge(data_table),
        .gint_tabitem_page_venn(data_table),
        .gint_tabitem_page_heatmap(data_table),
        .gint_tabitem_page_volcano(data_table),
        .gint_tabitem_page_network(data_table),
        .gint_tabitem_page_circos(data_table)
      )
    )
  ,skin = "purple")
}





