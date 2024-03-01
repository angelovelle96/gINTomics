
.gint_dashboardsidebar <- function(data_table){
  myImgResources <- "imgResources/logo_gINTomics.png"
  addResourcePath(prefix = "imgResources", directoryPath = "inst/www/")
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Genomic Integration", tabName = "genomicIntegrationPage", startExpanded = TRUE,
               menuSubItem("Coefficients Distribution", tabName = "coefDistribGenomic"),
               menuSubItem("Heatmap", tabName = "heatmapGenomic"),
               menuSubItem("Chromosome Distribution", tabName = "histoGenomic"),
               menuSubItem("Enrichment", tabName = "enrichGenomic")),
      menuItem("Transcription Integration", tabName = "transcriptionIntegrationPage", startExpanded = TRUE,
               menuSubItem("Coefficients Distribution", tabName = "coefDistribTranscript"),
               menuSubItem("Chromosome Distribution", tabName = "histoTranscript"),
               menuSubItem("Network", tabName = "networkTranscript"),
               menuSubItem("Enrichment", tabName = "enrichTrascript")),
      menuItem("DEGs", tabName = "degsPage", startExpanded = TRUE,
               menuSubItem("Coefficients Distribution", tabName = "coefDistribDEGs"),
               menuSubItem("Heatmap", tabName = "heatmapDEGs"),
               menuSubItem("Chromosome Distribution", tabName = "histoDEGs"),
               menuSubItem("Network", tabName = "networkDEGs")),
      menuItem("Complete Integration", tabName = "completeIntegrationPage", startExpanded = TRUE,
               menuSubItem("Circos Plots", tabName = "circosIntegration"),
               menuSubItem("Data Table", tabName = "fullTable")),
      menuItem("Class Comparison", tabName = "compareClassPage")
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

.gint_subItem_coefDistribGenomic <- function(data_table){
  tabItem(tabName = "coefDistribGenomic",
          fluidRow(
            box(title = "Page containing coefs distribution plots (Venn, Volcano and RidgeLine plots).",
                "1. ricordare di eliminare bottone deg (le analisi deg vanno nella sezione apposita"),
            mainPanel(
              tabsetPanel(type = "tabs",
                          tabPanel("Venn Diagram",

              sidebarLayout(
                sidebarPanel(
                  selectInput("classSelectVenn",
                              "Select Class:",
                              choices = unique(data_table$class)),
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
                    sliderInput("fdrRangeVenn",
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005))
                ),
               mainPanel(
                 plotlyOutput("venn_plot"),
                 dataTableOutput("common_genes_table")
               )
              )
            ),
            tabPanel("Volcano Plot",
                     sidebarLayout(
                       sidebarPanel(
                         selectInput("genomicIntegrationSelectVolcano",
                                     label = "Gene/miRNA:",
                                     choices = intersect(unique(data_table$omics), c("gene_genomic_res", "cnv_gene_res", "met_gene_res", "mirna_cnv_res"))),
                         conditionalPanel(
                           condition = "input.genomicIntegrationSelectVolcano=='gene_genomic_res'",
                           selectInput(inputId = 'genomicTypeSelectVolcano',
                                       label = 'Integration Type:',
                                       choices = intersect(unique(data_table$cnv_met), c("cnv", "met")))),
                         selectInput(inputId = 'genomicSignificativityCriteriaVolcano',
                                     label = 'Significativity criteria:',
                                     choices = c('pval', 'FDR')),
                         conditionalPanel(
                           condition = "input.genomicSignificativityCriteriaVolcano == 'pval'",
                           sliderInput("genomicPvalRangeVolcano",
                                       "P-Value Range:",
                                       min = 0,
                                       max = 1,
                                       value = c(0.05),
                                       step = 0.005)),
                         conditionalPanel(
                           condition = "input.genomicSignificativityCriteriaVolcano == 'FDR'",
                           sliderInput("genomicFdrRangeVolcano",
                                       "FDR-Value Range:",
                                       min = 0,
                                       max = 1,
                                       value = c(0.05),
                                       step = 0.005))
                       ),
                       mainPanel(
                         plotlyOutput('volcanoPlot')
                       )
                     )),
            tabPanel("RidgeLine Plot",
                     sidebarLayout(
                       sidebarPanel(
                         selectInput("genomicClassSelectRidge",
                                     "Select Class:",
                                     choices = unique(data_table$class)),
                         selectInput("genomicIntegrationSelectRidge",
                                     label = "Gene/miRNA:",
                                     choices = intersect(unique(data_table$omics), c("gene_genomic_res", "cnv_gene_res", "met_gene_res", "mirna_cnv_res"))),
                         conditionalPanel(
                           condition = "input.genomicIntegrationSelectRidge=='gene_genomic_res'",

                           selectInput(inputId = 'genomicTypeSelectRidge',
                                       label = 'Genomic Type:',
                                       choices = intersect(unique(data_table$cnv_met), c("met", "cnv")))),
                         selectInput(inputId = 'genomicSignificativityCriteriaRidge',
                                     label = 'Significativity criteria:',
                                     choices = c('pval', 'FDR')),
                         conditionalPanel(
                           condition = "input.genomicSignificativityCriteriaRidge == 'pval'",
                           sliderInput("genomicPvalRangeRidge",
                                       "P-Value Range:",
                                       min = 0,
                                       max = 1,
                                       value = c(0, 0.05),
                                       step = 0.005)),
                         conditionalPanel(
                           condition = "input.genomicSignificativityCriteriaRidge == 'FDR'",
                           sliderInput("genomicFdrRangeRidge",
                                       "FDR-Value Range:",
                                       min = 0,
                                       max = 1,
                                       value = c(0, 0.05),
                                       step = 0.005))
                       ),
                       mainPanel(
                         plotOutput("ridgelinePlot"),
                         DT::dataTableOutput('ridgelineTable')
                       )
                     ))
          )
  )
          )
  )
}

####################################################
#####################################################

.gint_subItem_HeatmapGenomic  <- function(data_table){

  tabItem(tabName = "heatmapGenomic",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(

              sidebarLayout(
                sidebarPanel(
                  selectInput(inputId = 'integrationSelectHeatmap',
                              label = 'Integration Type:',
                              choices = intersect(c("gene_genomic_res", "gene_met_res",
                                                    "gene_cnv_res", "mirna_cnv_res"),
                                                  unique(data_table$omics))),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmap == 'gene_genomic_res'",
                    sliderInput("numTopGenesHeatmapCNV",
                                "Number of top genes (CNV):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10),
                    sliderInput("numTopGenesHeatmapMET",
                                "Number of top genes (MET):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmap == 'cnv_gene_res'",
                    sliderInput("numTopGenesHeatmapCNVonly",
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmap == 'met_gene_res'",
                    sliderInput("numTopGenesHeatmapMETonly",
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmap == 'mirna_cnv_res'",
                    sliderInput("numTopGenesHeatmapmirna_cnv",
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  selectInput("selectClassHeatmap",
                              "Select the class:",
                              choices = unique(data_table$class),
                              multiple = FALSE),
                  selectInput(inputId = 'significativityCriteriaHeatmap',
                              label = 'Significativity criteria:',
                              choices = c('pval', 'FDR')),
                  conditionalPanel(
                    condition = "input.significativityCriteriaHeatmap == 'pval'",
                    sliderInput("pvalRangeHeatmap",
                                "P-Value Range:",
                                min = 0.01,
                                max = 0.5,
                                value = 0.05,
                                step = 0.005)),
                  conditionalPanel(
                    condition = "input.significativityCriteriaHeatmap == 'FDR'",
                    sliderInput("FDRRangeHeatmap",
                                "FDR-Value Range:",
                                min = 0.01,
                                max = 0.5,
                                value = 0.05,
                                step = 0.005))

                ),
                mainPanel(
                  InteractiveComplexHeatmapOutput('heatmap')
                )

              )
            )
          )
  )
}

##########################################################
###########################################################

.gint_subItem_chrDistribGenomic <- function(data_table){
  tabItem(tabName = "histoGenomic",
          fluidRow(
            box(title = "Page containing data table",
                "This is the content of Page 1."),
            mainPanel(
              tabsetPanel(type = 'tabs',
                          tabPanel('XXX',
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectInput(inputId = 'genomicIntegrationSelectHisto',
                                                   label = 'Integration Type:',
                                                   choices = intersect(c("gene_genomic_res", "gene_met_res",
                                                                         "gene_cnv_res", "mirna_cnv_res"),
                                                                       unique(data_table$omics))),
                                       conditionalPanel(
                                         condition = "input.genomicIntegrationSelectHisto=='gene_genomic_res'",
                                         selectInput(inputId = "genomicTypeSelect",
                                                     label = "Type Selection",
                                                     choices = intersect(unique(data_table$cnv_met), c("met", "cnv")))),
                                       selectInput(inputId = 'genomicClassSelectHisto',
                                                   label = 'Class:',
                                                   choices = unique(data_table$class)),
                                       selectInput(inputId = 'genomicChrSelectHisto',
                                                   label = 'Chr:',
                                                   choices = c('All', unique(data_table$chr_cov))),
                                       selectInput(inputId = 'genomicSignificativityCriteriaHisto',
                                                   label = 'Significativity criteria:',
                                                   choices = c('pval', 'FDR')),
                                       conditionalPanel(
                                         condition = "input.genomicSignificativityCriteriaHisto == 'pval'",
                                         sliderInput("genomicPvalRangeHisto",
                                                     "P-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.genomicSignificativityCriteriaHisto == 'FDR'",
                                         sliderInput("genomicFdrRangeHisto",
                                                     "FDR-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005))
                                     ),
                                     mainPanel(
                                       plotlyOutput('histogramPlot'),
                                       DT::dataTableOutput('histogramTable')
                                     )
                                   )
                          ),
                          tabPanel('TFs',
                                   sidebarLayout(
                                     sidebarPanel(
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

#############################################################
##############################################################

.gint_tabitem_enr <- function(data_table){
  ns <- NS("prova")
  tabItem(tabName = "enrichGenomic",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(

                         sidebarLayout(
                           sidebarPanel(
                             #input
                             selectInput(inputId = ns('genomicTypeSelectEnrich'),
                                         label = 'Integration Type:',
                                         choices = unique(data_table$cnv_met)[!is.na(unique(data_table$cnv_met))]),
                             selectInput(inputId = 'ClassSelectEnrich',
                                         label = 'Class:',
                                         choices = unique(data_table$class))

                           ),
                           mainPanel(textOutput(ns("gen_enrichment")),
                                     plotlyOutput(ns("gen_dotplot"))))
                )


          ))
}

##################################################################
###################################################################



.gint_subItem_coefDistribTranscript <- function(data_table){
  tabItem(tabName = "coefDistribTranscript",
          fluidRow(
            box(title = "Page containing coefs distribution plots (Venn, Volcano and RidgeLine plots).",
                "1. ricordare di eliminare bottone deg (le analisi deg vanno nella sezione apposita"),
            mainPanel(
              tabsetPanel(type = "tabs",
                          tabPanel("Volcano Plot",
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectInput("transcriptIntegrationSelectVolcano",
                                                   label = "Gene/miRNA:",
                                                   choices = intersect(unique(data_table$omics), c("tf_res", "tf_mirna_res", "mirna_target_res"))),
                                       selectInput(inputId = 'transcriptSignificativityCriteriaVolcano',
                                                   label = 'Significativity criteria:',
                                                   choices = c('pval', 'FDR')),
                                       conditionalPanel(
                                         condition = "input.transcriptSignificativityCriteriaVolcano == 'pval'",
                                         sliderInput("transcriptPvalRangeVolcano",
                                                     "P-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.transcriptSignificativityCriteriaVolcano == 'FDR'",
                                         sliderInput("transcriptFdrRangeVolcano",
                                                     "FDR-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0.05),
                                                     step = 0.005))
                                     ),
                                     mainPanel(
                                       plotlyOutput('transcriptVolcanoPlot')
                                     )
                                   )),
                          tabPanel("RidgeLine Plot",
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectInput(inputId = 'transcriptIntegrationSelectRidge',
                                                   label = 'Integration Type:',
                                                   choices = intersect(unique(data_table$omics), c("tf_res", "tf_mirna_res", "mirna_target_res"))),
                                       selectInput("transcriptClassSelectRidge",
                                                   "Select Class:",
                                                   choices = unique(data_table$class)),

                                       selectInput(inputId = 'transcriptSignificativityCriteriaRidge',
                                                   label = 'Significativity criteria:',
                                                   choices = c('pval', 'FDR')),
                                       conditionalPanel(
                                         condition = "input.transcriptSignificativityCriteriaRidge == 'pval'",
                                         sliderInput("transcriptPvalRangeRidge",
                                                     "P-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.transcriptSignificativityCriteriaRidge == 'FDR'",
                                         sliderInput("transcriptFdrRangeRidge",
                                                     "FDR-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005))
                                     ),
                                     mainPanel(
                                       plotOutput("ridgelinePlotTranscript"),
                                       DT::dataTableOutput('ridgelineTableTranscript')
                                     )
                                   ))
              )
            )
          )
  )
}
################################################################
#################################################################

.gint_subItem_chrDistribTranscript <- function(data_table){
  tabItem(tabName = "histoTranscript",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
              tabsetPanel(
                type = 'tabs',
                tabPanel('XXX',
                         sidebarLayout(
                           sidebarPanel(
                             selectInput(inputId = 'transcriptIntegrationSelectHisto',
                                         label = 'Integration Type:',
                                         choices = intersect(c("tf_res", "tf_mirna_res", "mirna_target_res"),
                                                             unique(data_table$omics))),
                             selectInput(inputId = 'transcriptClassSelectHisto',
                                         label = 'Class:',
                                         choices = unique(data_table$class)),
                             selectInput(inputId = 'transcriptChrSelectHisto',
                                         label = 'Chr:',
                                         choices = c('All', unique(data_table$chr_cov))),
                             selectInput(inputId = 'transcriptSignificativityCriteriaHisto',
                                         label = 'Significativity criteria:',
                                         choices = c('pval', 'FDR')),
                             conditionalPanel(
                               condition = "input.transcriptSignificativityCriteriaHisto == 'pval'",
                               sliderInput("transcriptPvalRangeHisto",
                                           "P-Value Range:",
                                           min = 0,
                                           max = 1,
                                           value = c(0, 0.05),
                                           step = 0.005)),
                             conditionalPanel(
                               condition = "input.transcriptSignificativityCriteriaHisto == 'FDR'",
                               sliderInput("transcriptFdrRangeHisto",
                                           "FDR-Value Range:",
                                           min = 0,
                                           max = 1,
                                           value = c(0, 0.05),
                                           step = 0.005))
                           ),
                           mainPanel(
                             plotlyOutput('histogramPlotTranscript'),
                             DT::dataTableOutput('histogramTableTranscript')
                           )
                         )
                )
              )
        )
      )
  )
}
########################################################################
#########################################################################

.gint_subItem_networkTranscript <- function(data_table){
  tabItem(tabName = "networkTranscript",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
              sliderInput("numNodes",
                          label = "Number of Nodes",
                          min = 10,
                          max = nrow(data_table),
                          value = 300),
              selectInput(inputId = 'significativityCriteriaNetwork',
                          label = 'Significativity criteria:',
                          choices = c('pval', 'FDR')),
              conditionalPanel(
                condition = "input.significativityCriteriaNetwork == 'pval'",
                sliderInput("pvalNetwork",
                            "P-Value:",
                            min = 0,
                            max = 1,
                            value = c(0.05),
                            step = 0.005)),
              conditionalPanel(
                condition = "input.significativityCriteriaNetwork == 'FDR'",
                sliderInput("fdrNetwork",
                            "FDR:",
                            min = 0,
                            max = 1,
                            value = c(0.05),
                            step = 0.005)),
              checkboxInput("layoutNetwork",
                           label = "Switch to tree Layout:",
                           value = FALSE),
              checkboxInput("physics",
                            label = "Physics",
                            value = FALSE),
              visNetworkOutput("networkPlot",
                               height = 800,
                               width = 1600)
            )
      )
  )
}

# ##################################################################
# ###################################################################
#
# .gint_subItem_enrichTranscript <- function(data_table){
#   tabItem(tabName = "enrichTrascript",
#           fluidRow(
#             box(title = "Subsection 1",
#                 "This is the content of Subsection 1."),
#             mainPanel(
#               tabsetPanel(
#                 type = 'tabs',
#                 tabPanel('XXXX',
#                          sliderInput("numNodes2",
#                                      label = "Number of Nodes",
#                                      min = 10,
#                                      max = nrow(data_table),
#                                      value = 300),
#                          selectInput(inputId = 'significativityCriteriaNetwork2',
#                                      label = 'Significativity criteria:',
#                                      choices = c('pval', 'FDR')),
#                          conditionalPanel(
#                            condition = "input.significativityCriteriaNetwork == 'pval'",
#                            sliderInput("pvalNetwork2",
#                                        "P-Value:",
#                                        min = 0,
#                                        max = 1,
#                                        value = c(0.05),
#                                        step = 0.005)),
#                          conditionalPanel(
#                            condition = "input.significativityCriteriaNetwork == 'FDR'",
#                            sliderInput("fdrNetwork2",
#                                        "FDR:",
#                                        min = 0,
#                                        max = 1,
#                                        value = c(0.05),
#                                        step = 0.005)),
#                          checkboxInput("layoutNetwork2",
#                                        label = "Switch to tree Layout:",
#                                        value = FALSE),
#                          checkboxInput("degNetwork2",
#                                        label = "Visualize only DEGs:",
#                                        value = FALSE),
#                          checkboxInput("physics2",
#                                        label = "Physics",
#                                        value = FALSE),
#                          visNetworkOutput("networkPlot",
#                                           height = 800,
#                                           width = 1600)
#                 )
#                 , tabPanel('DEGs'))
#             )
#           )
#   )
# }
##################################################################
###################################################################

.gint_subItem_coefDistribDEGs <- function(data_table){
  tabItem(tabName = "coefDistribDEGs",
          fluidRow(
            box(title = "Page containing coefs distribution plots (Venn, Volcano and RidgeLine plots).",
                "1. ricordare di eliminare bottone deg (le analisi deg vanno nella sezione apposita"),
            mainPanel(
              tabsetPanel(type = "tabs",
                          tabPanel("Venn Diagram",
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectInput("classSelectVennDEG",
                                                   "Select Class:",
                                                   choices = unique(data_table$class)),
                                       selectInput(inputId = 'significativityCriteriaVennDEG',
                                                   label = 'Significativity criteria:',
                                                   choices = c('pval', 'FDR')),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaVennDEG == 'pval'",
                                         sliderInput("pvalRangeVennDEG",
                                                     "P-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaVennDEG == 'FDR'",
                                         sliderInput("fdrRangeVennDEG",
                                                     "FDR-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005))
                                     ),
                                     mainPanel(
                                       plotlyOutput("venn_plotDEG"),
                                       dataTableOutput("venn_tableDEG")
                                     )
                                   )
                          ),
                          tabPanel("Volcano Plot",
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectInput("integrationSelectVolcanoDEG",
                                                   label = "Gene/miRNA:",
                                                   choices = unique(data_table$omics)),
                                       conditionalPanel(
                                         condition = "input.integrationSelectVolcanoDEG=='gene_genomic_res'",
                                         selectInput(inputId = 'typeSelectVolcanoDEG',
                                                     label = 'Integration Type:',
                                                     choices = intersect(unique(data_table$cnv_met), c("cnv", "met")))),
                                       selectInput(inputId = 'significativityCriteriaVolcanoDEG',
                                                   label = 'Significativity criteria:',
                                                   choices = c('pval', 'FDR')),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaVolcanoDEG == 'pval'",
                                         sliderInput("pvalRangeVolcanoDEG",
                                                     "P-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaVolcanoDEG == 'FDR'",
                                         sliderInput("fdrRangeVolcanoDEG",
                                                     "FDR-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0.05),
                                                     step = 0.005))
                                     ),
                                     mainPanel(
                                       plotlyOutput('volcanoPlotDEG')
                                     )
                                   )),
                          tabPanel("RidgeLine Plot",
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectInput("integrationSelectRidgeDEG",
                                                   label = "Gene/miRNA:",
                                                   choices = unique(data_table$omics)),
                                       conditionalPanel(
                                         condition = "input.integrationSelectRidgeDEG=='gene_genomic_res'",

                                         selectInput(inputId = 'typeSelectRidgeDEG',
                                                     label = 'Genomic Type:',
                                                     choices = intersect(unique(data_table$cnv_met), c("met", "cnv")))),
                                       selectInput(inputId = 'significativityCriteriaRidgeDEG',
                                                   label = 'Significativity criteria:',
                                                   choices = c('pval', 'FDR')),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaRidgeDEG == 'pval'",
                                         sliderInput("pvalRangeRidgeDEG",
                                                     "P-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaRidgeDEG == 'FDR'",
                                         sliderInput("fdrRangeRidgeDEG",
                                                     "FDR-Value Range:",
                                                     min = 0,
                                                     max = 1,
                                                     value = c(0, 0.05),
                                                     step = 0.005))
                                     ),
                                     mainPanel(
                                       plotOutput("ridgelinePlotDEG"),
                                       DT::dataTableOutput('ridgelineTableDEG')
                                     )
                                   ))
              )
            )
          )
  )
}

##################################################################
###################################################################

.gint_subItem_HeatmapDEGs <- function(data_table){
  tabItem(tabName = "heatmapDEGs",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(

              sidebarLayout(
                sidebarPanel(
                  selectInput(inputId = 'integrationSelectHeatmapDEG',
                              label = 'Integration Type:',
                              choices = unique(data_table$omics)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmapDEG == 'gene_genomic_res'",
                    sliderInput("numTopGenesHeatmapCNVDEG",
                                "Number of top genes (CNV):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10),
                    sliderInput("numTopGenesHeatmapMETDEG",
                                "Number of top genes (MET):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmapDEG == 'cnv_gene_res'",
                    sliderInput("numTopGenesHeatmapCNVonlyDEG",
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmapDEG == 'met_gene_res'",
                    sliderInput("numTopGenesHeatmapMETonlyDEG",
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  conditionalPanel(
                    condition = "input.integrationSelectHeatmapDEG == 'mirna_cnv_res'",
                    sliderInput("numTopGenesHeatmapmirna_cnvDEG",
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10)),
                  selectInput("selectClassHeatmapDEG",
                              "Select the class:",
                              choices = unique(data_table$class),
                              multiple = FALSE),
                  selectInput(inputId = 'significativityCriteriaHeatmapDEG',
                              label = 'Significativity criteria:',
                              choices = c('pval', 'FDR')),
                  conditionalPanel(
                    condition = "input.significativityCriteriaHeatmapDEG == 'pval'",
                    sliderInput("pvalRangeHeatmapDEG",
                                "P-Value Range:",
                                min = 0.01,
                                max = 0.5,
                                value = 0.05,
                                step = 0.005)),
                  conditionalPanel(
                    condition = "input.significativityCriteriaHeatmapDEG == 'FDR'",
                    sliderInput("FDRRangeHeatmapDEG",
                                "FDR-Value Range:",
                                min = 0.01,
                                max = 0.5,
                                value = 0.05,
                                step = 0.005))

                ),
                mainPanel(
                  InteractiveComplexHeatmapOutput('heatmapDEG')
                )

              )
            )
          )
  )
}

#############################################################
##############################################################

.gint_subItem_chrDistribDEGs <- function(data_table){
  tabItem(tabName = "histoDEGs",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
              tabsetPanel(
                type = 'tabs',
                tabPanel('XXX',
                         sidebarLayout(
                           sidebarPanel(
                             selectInput(inputId = 'integrationSelectHistoDEG',
                                         label = 'Integration Type:',
                                         choices = unique(data_table$omics)),
                             conditionalPanel(
                               condition = "input.integrationSelectHistoDEG=='gene_genomic_res'",
                               selectInput(inputId = "typeSelectDEG",
                                           label = "Type Selection",
                                           choices = intersect(unique(data_table$cnv_met), c("met", "cnv")))),
                             selectInput(inputId = 'classSelectHistoDEG',
                                         label = 'Class:',
                                         choices = unique(data_table$class)),
                             selectInput(inputId = 'chrSelectHistoDEG',
                                         label = 'Chr:',
                                         choices = c('All', unique(data_table$chr_cov))),
                             selectInput(inputId = 'significativityCriteriaHistoDEG',
                                         label = 'Significativity criteria:',
                                         choices = c('pval', 'FDR')),
                             conditionalPanel(
                               condition = "input.significativityCriteriaHistoDEG == 'pval'",
                               sliderInput("pvalRangeHistoDEG",
                                           "P-Value Range:",
                                           min = 0,
                                           max = 1,
                                           value = c(0, 0.05),
                                           step = 0.005)),
                             conditionalPanel(
                               condition = "input.significativityCriteriaHistoDEG == 'FDR'",
                               sliderInput("fdrRangeHistoDEG",
                                           "FDR-Value Range:",
                                           min = 0,
                                           max = 1,
                                           value = c(0, 0.05),
                                           step = 0.005))
                           ),
                           mainPanel(
                             plotlyOutput('histogramPlotDEG'),
                             DT::dataTableOutput('histogramTableDEG')
                           )
                         )
                )
              )
            )
          )
  )
}
#############################################################
##############################################################

.gint_subItem_networkDEGs <- function(data_table){
  tabItem(tabName = "networkDEGs",
          fluidRow(
            box(title = "Subsection 1",
                "This is the content of Subsection 1."),
            mainPanel(
                         sliderInput("numNodesDEG",
                                     label = "Number of Nodes",
                                     min = 10,
                                     max = nrow(data_table),
                                     value = 300),
                         selectInput(inputId = 'significativityCriteriaNetworkDEG',
                                     label = 'Significativity criteria:',
                                     choices = c('pval', 'FDR')),
                         conditionalPanel(
                           condition = "input.significativityCriteriaNetworkDEG == 'pval'",
                           sliderInput("pvalNetworkDEG",
                                       "P-Value:",
                                       min = 0,
                                       max = 1,
                                       value = c(0.05),
                                       step = 0.005)),
                         conditionalPanel(
                           condition = "input.significativityCriteriaNetworkDEG == 'FDR'",
                           sliderInput("fdrNetworkDEG",
                                       "FDR:",
                                       min = 0,
                                       max = 1,
                                       value = c(0.05),
                                       step = 0.005)),
                         checkboxInput("layoutNetworkDEG",
                                       label = "Switch to tree Layout:",
                                       value = FALSE),
                         checkboxInput("physicsDEG",
                                       label = "Physics",
                                       value = FALSE),
                         visNetworkOutput("networkPlotDEG",
                                          height = 800,
                                          width = 1600)
                )
          )
  )
}

#############################################################
##############################################################
.gint_subItem_circosCompleteInt <- function(data_table){
  # tabItem(tabName = "page_circos",
  #         fluidRow(
  #           box(title = "Subsection 1",
  #               "This is the content of Subsection 1."),
  #           mainPanel(
  #               do.call("tabsetPanel", lapply(as.list(unique(data_table$omics)), function(x) tabPanel(x, goslingOutput(paste0("circos_",x)))))
  #       )          ## USARE uiOutput("gosling_plot_circos_ui")
  #     )
  # )
  chr <- unique(data_table$chr_response)
  chr <- c("All", paste0("chr", gtools::mixedsort(chr[!is.na(chr)])))
  tabItem(tabName = "circosIntegration",fluidPage(use_gosling(),
            sidebarLayout(
              sidebarPanel(
                selectInput(inputId = "circosClass",choices = unique(data_table$class), label = "Class"),
                selectInput(inputId = "circosType",choices = c("Gene", "miRNA"), label = "Gene/miRNA"),
                selectInput(inputId = "circosChr",choices = chr, label = "Chromosome"),
                actionButton("circosDownload_png","PNG",icon = icon("cloud-arrow-down")),
                actionButton("circosDownload_pdf","PDF",icon = icon("cloud-arrow-down")),
                selectInput(inputId = "circosLayout",choices = c("circular", "linear"), label = "Layout"), width = 3),
              mainPanel(column(12,goslingOutput('gosling_plot_circos'),
                               div(style="height: 200px;")))),
            tags$head(tags$style(
              HTML('.container-fluid {background-color: #ffffff;}')
            )
          )))



}
# '/* body */.content-wrapper, .right-side {background-color: #7d848b;}'
#"body {background-color: #ffffff; /* white background color */}"
#############################################################
##############################################################
#
# .gint_tabitem_enr <- function(data_table){
#   ns <- NS("prova")
#   tabItem(tabName = "page_enrichment",
#           fluidRow(
#             box(title = "Subsection 1",
#                 "This is the content of Subsection 1."),
#             mainPanel(
#               tabsetPanel(
#                 type = 'tabs',
#                 tabPanel('Genomic',
#                          sidebarLayout(
#                            sidebarPanel(
#                              #input
#                              selectInput(inputId = ns('genomicTypeSelectEnrich'),
#                                          label = 'Integration Type:',
#                                          choices = unique(data_table$cnv_met)[!is.na(unique(data_table$cnv_met))]),
#                              selectInput(inputId = 'ClassSelectEnrich',
#                                          label = 'Class:',
#                                          choices = unique(data_table$class))
#
#                            ),
#                          mainPanel(textOutput(ns("gen_enrichment")),
#                                    plotlyOutput(ns("gen_dotplot"))))
#                         )
#                 , tabPanel('TF',
#                            sidebarLayout(
#                              sidebarPanel(
#                                #input
#                                selectInput(inputId = 'ClassSelectEnrich',
#                                            label = 'Class:',
#                                            choices = unique(data_table$class))
#
#                              ),
#                              mainPanel(textOutput(ns("tf_enrichment")),
#                                        plotlyOutput(ns("tf_dotplot"))))))
#             )
#           ))
# }
################################################################
#################################################################

.gint_subItem_tableCompleteInt <- function(data_table){
  tabItem(tabName = "fullTable",
          fluidRow(
            box(title = "Page containing data tables.",
                "1. ."),
            mainPanel(
                          sidebarLayout(
                                     sidebarPanel(
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
                                                     min = 0.001,
                                                     max = 1,
                                                     value = c(0.001, 0.05),
                                                     step = 0.005)),
                                       conditionalPanel(
                                         condition = "input.significativityCriteriaTable == 'FDR'",
                                         sliderInput("FDRRangeTable",
                                                     "FDR-Value Range:",
                                                     min = 0.001,
                                                     max = 1,
                                                     value = c(0.001, 0.05),
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

#############################################################
##############################################################

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
        .gint_subItem_coefDistribGenomic(data_table),
        .gint_subItem_HeatmapGenomic(data_table),
        .gint_subItem_chrDistribGenomic(data_table),
        .gint_tabitem_enr(data_table),
         .gint_subItem_coefDistribTranscript(data_table),
         .gint_subItem_chrDistribTranscript(data_table),
         .gint_subItem_networkTranscript(data_table),
        # .gint_subItem_enrichTranscript(data_table),
         .gint_subItem_coefDistribDEGs(data_table),
         .gint_subItem_HeatmapDEGs(data_table),
         .gint_subItem_chrDistribDEGs(data_table),
         .gint_subItem_networkDEGs(data_table),
         .gint_subItem_circosCompleteInt(data_table),
         .gint_subItem_tableCompleteInt(data_table)
      )
    )
  ,skin = "purple")
}

