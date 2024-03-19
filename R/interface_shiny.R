# Integration Dashboard Sidebar
#' @importFrom shiny addResourcePath tags img icon
#' @importFrom shinydashboard dashboardSidebar sidebarMenu menuItem menuSubItem
.gint_dashboardsidebar <- function(data_table) {
    myImgResources <- "imgResources/logo_gINTomics.png"
    addResourcePath(
        prefix = "imgResources",
        system.file(
            directoryPath = "www/",
            package = "gINTomics"
        )
    )
    dashboardSidebar(
        sidebarMenu(
            menuItem("Home",
                tabName = "home", icon = icon("home")
            ),
            menuItem("Genomic Integration",
                tabName = "genomicIntegrationPage",
                startExpanded = TRUE,
                menuSubItem("Coefficients Distribution",
                    tabName =
                        "coefDistribGenomic"
                ),
                menuSubItem("Heatmap",
                    tabName = "heatmapGenomic"
                ),
                menuSubItem("Chromosome Distribution",
                    tabName = "histoGenomic"
                ),
                menuSubItem("Enrichment",
                    tabName = "enrichGenomic"
                )
            ),
            menuItem("Transcription Integration",
                tabName =
                    "transcriptionIntegrationPage",
                startExpanded = TRUE,
                menuSubItem("Coefficients Distribution",
                    tabName =
                        "coefDistribTranscript"
                ),
                menuSubItem("Chromosome Distribution",
                    tabName =
                        "histoTranscript"
                ),
                menuSubItem("Network",
                    tabName =
                        "networkTranscript"
                ),
                menuSubItem("Enrichment",
                    tabName =
                        "enrichTranscript"
                )
            ),
            menuItem("Class Comparison",
                tabName = "degsPage",
                startExpanded = TRUE,
                menuSubItem("Coefficients Distribution",
                    tabName =
                        "coefDistribDEGs"
                ),
                menuSubItem("Heatmap",
                    tabName = "heatmapDEGs"
                ),
                menuSubItem("Chromosome Distribution",
                    tabName = "histoDEGs"
                ),
                menuSubItem("Network",
                    tabName = "networkDEGs"
                )
            ),
            menuItem("Complete Integration",
                tabName = "completeIntegrationPage",
                startExpanded = TRUE,
                menuSubItem("Circos Plots",
                    tabName =
                        "circosIntegration"
                ),
                menuSubItem("Data Table",
                    tabName = "fullTable"
                )
            )
        ),
        img(
            src = myImgResources,
            height = 100, width = 100,
            align = "center"
        ),
        tags$style(".left-side, .main-sidebar {padding-top: 80px}")
    )
}


# Home Tab Item for gINTomics
#' @importFrom shiny addResourcePath fluidPage tags img fluidRow column
#' @importFrom shinydashboard tabItem box
.gint_tabitem_home <- function(data_table) {
    myImgResources <- "imgResources/home.png"
    addResourcePath(prefix = "imgResources", system.file(
        directoryPath = "www/",
        package = "gINTomics"
    ))
    tabItem(
        tabName = "home", fluidRow(column(width = 12, box(
            title =
                "gINTomics
                                                            Visualizer 1.0",
            footer = HTML("
                Welcome to gINTomics Visualizer 1.0. <br> <br>
                This interactive environment, based on Shiny, allows the
                user to easily visualize output results from gINTomics
                package and to save them for downstream tasks and reports
                creation."),
            width = 12
        ))), fluidRow(column(
            width = 12,
            div(
                img(
                    src = myImgResources,
                    height = 300 * 3,
                    width = 400 * 3
                ),
                class = "text-center"
            )
        )),
        .gint_home_box()
    )
}

# Coefficients Distribution Subitem for Genomic Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#'  mainPanel sidebarPanel plotOutput div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput
.gint_subItem_coefDistribGenomic <- function(data_table) {
    ns <- NS("venn_gen")
    ns2 <- NS("volcano_gen")
    ns3 <- NS("ridge_gen")
    tabItem(
        tabName = "coefDistribGenomic",
        fluidRow(
            .gint_coefDistrib_box(),
            mainPanel(
                tabsetPanel(
                    type = "tabs",
                    tabPanel(
                        "Venn Diagram",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns("ClassSelect"),
                                    "Select Class:",
                                    choices = unique(
                                        data_table$class
                                    )
                                ),
                                selectInput(
                                    inputId = ns(
                                        "SignificativityCriteria"
                                    ),
                                    label = "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                      "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns
                                )
                            ),
                            mainPanel(
                                plotlyOutput(ns("plotly")),
                                tags$div(
                                    style = "overflow-x: auto;",
                                    dataTableOutput(ns("table")),
                                    downloadButton(
                                        ns("download_csv"),
                                        "Download CSV"
                                    )
                                )
                            )
                        )
                    ),
                    tabPanel(
                        "Volcano Plot",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns2("IntegrationSelect"),
                                    label = "Gene/miRNA:",
                                    choices =
                                        .change_int_names(
                                            intersect(
                                                unique(
                                                    data_table$omics
                                                ),
                                                c(
                                                    "gene_genomic_res",
                                                    "gene_cnv_res",
                                                    "gene_met_res",
                                                    "mirna_cnv_res"
                                                )
                                            )
                                        )
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.IntegrationSelect==
                                    'gene_genomic_res'",
                                    selectInput(
                                        inputId = ns2(
                                            "genomicTypeSelect"
                                        ),
                                        label = "Integration Type:",
                                        choices =
                                            intersect(
                                                unique(
                                                    data_table$cnv_met
                                                ),
                                                c("cnv", "met")
                                            )
                                    ),
                                    ns = ns2
                                ),
                                selectInput(
                                    inputId = ns2(
                                        "SignificativityCriteria"
                                    ),
                                    label =
                                        "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria==
                                    'pval'",
                                    sliderInput(ns2("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0.05),
                                        step = 0.005
                                    ), ns = ns2
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria ==
                                    'FDR'",
                                    sliderInput(ns2("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0.05),
                                        step = 0.005
                                    ), ns = ns2
                                )
                            ),
                            mainPanel(
                                plotlyOutput(ns2("plotly")),
                                div(style = "height: 400px;")
                            )
                        )
                    ),
                    tabPanel(
                        "RidgeLine Plot",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns3("ClassSelect"),
                                    "Select Class:",
                                    choices = unique(
                                        data_table$class
                                    )
                                ),
                                selectInput(ns3("IntegrationSelect"),
                                    label = "Gene/miRNA:",
                                    choices =
                                        .change_int_names(
                                            intersect(
                                                unique(
                                                    data_table$omics
                                                ),
                                                c(
                                                    "gene_genomic_res",
                                                    "gene_cnv_res",
                                                    "gene_met_res",
                                                    "mirna_cnv_res"
                                                )
                                            )
                                        )
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.IntegrationSelect==
                                    'gene_genomic_res'",
                                    selectInput(
                                        inputId = ns3(
                                            "genomicTypeSelect"
                                        ),
                                        label = "Genomic Type:",
                                        choices = intersect(
                                            unique(
                                                data_table$cnv_met
                                            ),
                                            c("met", "cnv")
                                        )
                                    ),
                                    ns = ns3
                                ),
                                selectInput(
                                    inputId = ns3(
                                        "SignificativityCriteria"
                                    ),
                                    label =
                                        "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns3("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns3
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns3("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns3
                                )
                            ),
                            mainPanel(
                                plotOutput(ns3("plotly"), width = 700),
                                tags$div(
                                    style = "overflow-x: auto;",
                                    dataTableOutput(ns3("table")),
                                    downloadButton(
                                        ns3("download_csv"),
                                        "Download CSV"
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}

# Heatmap Subitem for Genomic Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#' mainPanel sidebarPanel
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom plotly plotlyOutput
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
.gint_subItem_HeatmapGenomic <- function(data_table) {
    ns <- NS("heat_gen")
    tabItem(
        tabName = "heatmapGenomic",
        fluidRow(
            .gint_heatmap_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("IntegrationSelect"),
                            label = "Integration Type:",
                            choices = .change_int_names(
                                intersect(
                                    c(
                                        "gene_genomic_res", "gene_met_res",
                                        "gene_cnv_res", "mirna_cnv_res"
                                    ),
                                    unique(data_table$omics)
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_genomic_res'",
                            sliderInput(ns("numTopGenesHeatmapCNV"),
                                "Number of top genes (CNV):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ),
                            sliderInput(ns("numTopGenesHeatmapMET"),
                                "Number of top genes (MET):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_cnv_res'",
                            sliderInput(ns("numTopGenesHeatmapCNVonly"),
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_met_res'",
                            sliderInput(ns("numTopGenesHeatmapMETonly"),
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'mirna_cnv_res'",
                            sliderInput(ns("numTopGenesHeatmapmirna_cnv"),
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        sliderInput(ns("numSamples"),
                            "Number of Samples:",
                            value = 50,
                            min = 1,
                            max = 1000,
                            step = 10
                        ),
                        selectInput(ns("ClassSelect"),
                            "Select the class:",
                            choices = unique(data_table$class),
                            multiple = FALSE
                        ),
                        selectInput(
                            inputId = ns("SignificativityCriteria"),
                            label = "Significativity criteria:",
                            choices = c("pval", "FDR")
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='pval'",
                            sliderInput(ns("PvalRange"),
                                "P-Value Range:",
                                min = 0.01,
                                max = 1,
                                value = 0.05,
                                step = 0.005
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='FDR'",
                            sliderInput(ns("FdrRange"),
                                "FDR-Value Range:",
                                min = 0.01,
                                max = 1,
                                value = 0.05,
                                step = 0.005
                            ), ns = ns
                        ),
                        radioButtons(ns("scaleHeatmap"),
                            "scale by:",
                            choices = c("row", "col", "none")
                        )
                    ),
                    mainPanel(
                        InteractiveComplexHeatmapOutput(ns("heatmap")),
                        dataTableOutput(ns("aa"))
                    )
                )
            )
        )
    )
}
# Chr Distribution Subitem for Genomic Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#' mainPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom gtools mixedsort
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput
#' @importFrom stats na.omit
.gint_subItem_chrDistribGenomic <- function(data_table) {
    ns <- NS("histo_gen")
    chr <- c("All", mixedsort(unique(na.omit(data_table$chr_cov))))
    tabItem(
        tabName = "histoGenomic",
        fluidRow(
            .gint_chrDistrib_box(),
            mainPanel(sidebarLayout(
                sidebarPanel(
                    selectInput(
                        inputId = ns("IntegrationSelect"),
                        label = "Integration Type:",
                        choices = .change_int_names(
                            intersect(
                                c(
                                    "gene_genomic_res", "gene_met_res",
                                    "gene_cnv_res", "mirna_cnv_res"
                                ),
                                unique(data_table$omics)
                            )
                        )
                    ),
                    conditionalPanel(
                        condition = "input.IntegrationSelect==
                        'gene_genomic_res'",
                        selectInput(
                            inputId = ns("genomicTypeSelect"),
                            label = "Type Selection",
                            choices = intersect(
                                unique(data_table$cnv_met),
                                c("met", "cnv")
                            )
                        ), ns = ns
                    ),
                    selectInput(
                        inputId = ns("ClassSelect"),
                        label = "Class:",
                        choices = unique(data_table$class)
                    ),
                    selectInput(
                        inputId = ns("ChrSelect"),
                        label = "Chr:",
                        choices = chr
                    ),
                    selectInput(
                        inputId = ns("SignificativityCriteria"),
                        label = "Significativity criteria:",
                        choices = c("pval", "FDR")
                    ),
                    conditionalPanel(
                        condition = "input.SignificativityCriteria=='pval'",
                        sliderInput(ns("PvalRange"),
                            "P-Value Range:",
                            min = 0,
                            max = 1,
                            value = c(0, 0.05),
                            step = 0.005
                        ), ns = ns
                    ),
                    conditionalPanel(
                        condition = "input.SignificativityCriteria=='FDR'",
                        sliderInput(ns("FdrRange"),
                            "FDR-Value Range:",
                            min = 0,
                            max = 1,
                            value = c(0, 0.05),
                            step = 0.005
                        ), ns = ns
                    )
                ),
                mainPanel(
                    plotlyOutput(ns("plotly")),
                    tags$div(
                        style = "overflow-x: auto;",
                        dataTableOutput(ns("table")),
                        downloadButton(ns("download_csv"), "Download CSV")
                    )
                )
            ))
        )
    )
}

# Enrichment Tab Item for Genomic Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS HTML tags fluidRow
#' mainPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput
.gint_tabitem_enr <- function(data_table) {
    ns <- NS("enrich_gen")
    tabItem(
        tabName = "enrichGenomic",
        fluidRow(
            .gint_enrich_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("genomicTypeSelect"),
                            label = "Integration Type:",
                            choices = unique(data_table$cnv_met)[
                                !is.na(unique(data_table$cnv_met))
                            ]
                        ),
                        selectInput(
                            inputId = ns("ClassSelect"),
                            label = "Class:",
                            choices = unique(data_table$class)
                        ),
                        selectInput(
                            inputId = ns("DBSelectEnrich"),
                            label = "Database:",
                            choices = c("go", "kegg")
                        )
                    ),
                    mainPanel(
                        textOutput(ns("check")),
                        plotlyOutput(ns("dotplot")),
                        HTML(paste0(rep("<br>", 20), collapse = "")),
                        tags$div(
                            style = "overflow-x: auto;",
                            dataTableOutput(ns("table"))
                        ),
                        downloadButton(ns("download_csv"), "Download CSV")
                    )
                )
            )
        )
    )
}

# Enrichment Subitem for Transcription Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput tabPanel NS htmlOutput textOutput tags fluidRow
#'  mainPanel uiOutput
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
.gint_subItem_enrichTranscript <- function(data_table) {
    ns <- NS("enrich_tf")
    tabItem(
        tabName = "enrichTranscript",
        fluidRow(
            .gint_enrich_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("ClassSelect"),
                            label = "Class:",
                            choices = unique(data_table$class)
                        ),
                        selectInput(
                            inputId = ns("DBSelectEnrich"),
                            label = "Database:",
                            choices = c("go", "kegg")
                        )
                    ),
                    mainPanel(
                        textOutput(ns("check")),
                        uiOutput(ns("dotplot"))
                    )
                )
            )
        )
    )
}

# Coefficients Distribution Subitem for Transcription Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#'  mainPanel sidebarPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput

.gint_subItem_coefDistribTranscript <- function(data_table) {
    ns <- NS("volcano_trans")
    ns2 <- NS("ridge_trans")
    tabItem(
        tabName = "coefDistribTranscript",
        fluidRow(
            .gint_coefDistrib_box(),
            mainPanel(
                tabsetPanel(
                    type = "tabs",
                    tabPanel(
                        "Volcano Plot",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns("IntegrationSelect"),
                                    label = "Gene/miRNA:",
                                    choices =
                                        .change_int_names(
                                            intersect(
                                                unique(
                                                    data_table$omics
                                                ),
                                                c(
                                                    "tf_res",
                                                    "tf_mirna_res",
                                                    "mirna_target_res"
                                                )
                                            )
                                        )
                                ),
                                selectInput(
                                    inputId = ns(
                                        "SignificativityCriteria"
                                    ),
                                    label = "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0.05),
                                        step = 0.005
                                    ), ns = ns
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0.05),
                                        step = 0.005
                                    ), ns = ns
                                )
                            ),
                            mainPanel(
                                plotlyOutput(ns("plotly")),
                                div(style = "height: 400px;")
                            )
                        )
                    ),
                    tabPanel(
                        "RidgeLine Plot",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(
                                    inputId = ns2(
                                        "IntegrationSelect"
                                    ),
                                    label = "Integration Type:",
                                    choices =
                                        .change_int_names(
                                            intersect(
                                                unique(data_table$omics),
                                                c(
                                                    "tf_res",
                                                    "tf_mirna_res",
                                                    "mirna_target_res"
                                                )
                                            )
                                        )
                                ),
                                selectInput(
                                    inputId = ns2("ClassSelect"),
                                    "Select Class:",
                                    choices =
                                        unique(data_table$class)
                                ),
                                selectInput(
                                    inputId =
                                        ns2(
                                            "SignificativityCriteria"
                                        ),
                                    label =
                                        "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns2("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns2
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns2("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns2
                                )
                            ),
                            mainPanel(
                                plotOutput(ns2("plotly"), width = 700),
                                tags$div(
                                    style = "overflow-x: auto;",
                                    dataTableOutput(ns2("table")),
                                    downloadButton(
                                        ns2("download_csv"),
                                        "Download CSV"
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}

# Chr Distribution Subitem for Transcription Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#' mainPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom gtools mixedsort
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput
#' @importFrom stats na.omit
.gint_subItem_chrDistribTranscript <- function(data_table) {
    ns <- NS("histo_trans")
    chr <- c("All", mixedsort(unique(na.omit(data_table$chr_cov))))
    tabItem(
        tabName = "histoTranscript",
        fluidRow(
            .gint_chrDistrib_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("IntegrationSelect"),
                            label = "Integration Type:",
                            choices = .change_int_names(
                                intersect(
                                    c(
                                        "tf_res",
                                        "tf_mirna_res",
                                        "mirna_target_res"
                                    ),
                                    unique(data_table$omics)
                                )
                            )
                        ),
                        selectInput(
                            inputId = ns("ClassSelect"),
                            label = "Class:",
                            choices = unique(data_table$class)
                        ),
                        selectInput(
                            inputId = ns("ChrSelect"),
                            label = "Chr:",
                            choices = chr
                        ),
                        selectInput(
                            inputId = ns("SignificativityCriteria"),
                            label = "Significativity criteria:",
                            choices = c("pval", "FDR")
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='pval'",
                            sliderInput(ns("PvalRange"),
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='FDR'",
                            sliderInput(ns("FdrRange"),
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005
                            ), ns = ns
                        )
                    ),
                    mainPanel(
                        plotlyOutput(ns("plotly")),
                        tags$div(
                            style = "overflow-x: auto;",
                            dataTableOutput(ns("table")),
                            downloadButton(ns("download_csv"), "Download CSV")
                        ),
                        plotlyOutput(ns("plotly_tf")),
                        dataTableOutput(ns("table_tf")),
                        downloadButton(ns("download_csv_tf"), "Download CSV")
                    )
                )
            )
        )
    )
}

# Network Subitem for Transcription Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#' mainPanel sidebarPanel div inputPanel
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom visNetwork visNetworkOutput

.gint_subItem_networkTranscript <- function(data_table) {
    ns <- NS("network_trans")
    tabItem(
        tabName = "networkTranscript",
        fluidRow(
            .gint_network_box(),
            mainPanel(
                sidebarLayout(
                    div(
                        style = "width: 1600px;",
                        inputPanel(
                            selectInput(ns("ClassSelect"),
                                label = "Select the Class:",
                                choices = unique(data_table$class)
                            ),
                            selectInput(
                                inputId = ns("SignificativityCriteria"),
                                label = "Significativity criteria:",
                                choices = c("pval", "FDR"),
                            ),
                            conditionalPanel(
                                condition =
                                  "input.SignificativityCriteria=='pval'",
                                sliderInput(ns("PvalRange"),
                                    "P-Value:",
                                    min = 0,
                                    max = 1,
                                    value = c(0.05),
                                    step = 0.005
                                ), ns = ns
                            ),
                            conditionalPanel(
                                condition =
                                  "input.SignificativityCriteria=='FDR'",
                                sliderInput(ns("FdrRange"),
                                    "FDR:",
                                    min = 0,
                                    max = 1,
                                    value = c(0.05),
                                    step = 0.005
                                )
                            ),
                            sliderInput(ns("numInteractions"),
                                label = "Number of Interactions",
                                min = 10,
                                max = nrow(data_table),
                                value = 200,
                                step = 50
                            ),
                            checkboxInput(ns("layout"),
                                label = "Switch to tree Layout:",
                                value = FALSE
                            ),
                            checkboxInput(ns("physics"),
                                label = "Physics",
                                value = FALSE
                            )
                        )
                    ),
                    mainPanel(
                        visNetworkOutput(ns("networkPlot"),
                            height = 800,
                            width = 1600
                        )
                    )
                )
            )
        )
    )
}

# Coefficients Distribution Subitem for Class Comparison
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#' mainPanel sidebarPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput

.gint_subItem_coefDistribDEGs <- function(data_table) {
    ns <- NS("venn_deg")
    ns2 <- NS("volcano_deg")
    ns3 <- NS("ridge_deg")
    tabItem(
        tabName = "coefDistribDEGs",
        fluidRow(
            .gint_coefDistrib_box(),
            mainPanel(
                tabsetPanel(
                    type = "tabs",
                    tabPanel(
                        "Venn Diagram",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns("ClassSelect"),
                                    "Select Class:",
                                    choices =
                                        unique(data_table$class)
                                ),
                                selectInput(
                                    inputId =
                                        ns(
                                            "SignificativityCriteria"
                                        ),
                                    label =
                                        "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    )
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    )
                                )
                            ),
                            mainPanel(
                                plotlyOutput(ns("plotly")),
                                dataTableOutput(ns("table")),
                                downloadButton(
                                    ns("download_csv"),
                                    "Download CSV"
                                )
                            )
                        )
                    ),
                    tabPanel(
                        "Volcano Plot",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns2("IntegrationSelect"),
                                    label = "Gene/miRNA:",
                                    choices =
                                        .change_int_names(
                                            unique(data_table$omics)
                                        )
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.IntegrationSelect==
                                    'gene_genomic_res'",
                                    selectInput(
                                        inputId = ns2(
                                            "genomicTypeSelect"
                                        ),
                                        label = "Integration Type:",
                                        choices = intersect(
                                            unique(
                                                data_table$cnv_met
                                            ),
                                            c("cnv", "met")
                                        )
                                    ),
                                    ns = ns2
                                ),
                                selectInput(
                                    inputId = ns2(
                                        "SignificativityCriteria"
                                    ),
                                    label =
                                        "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns2("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0.05),
                                        step = 0.005
                                    ), ns = ns2
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns2("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0.05),
                                        step = 0.005
                                    ), ns = ns2
                                )
                            ),
                            mainPanel(
                                plotlyOutput(ns2("plotly")),
                                div(style = "height: 400px;")
                            )
                        )
                    ),
                    tabPanel(
                        "RidgeLine Plot",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput(ns3("IntegrationSelect"),
                                    label = "Gene/miRNA:",
                                    choices =
                                        .change_int_names(
                                            unique(data_table$omics)
                                        )
                                ),
                                conditionalPanel(
                                    condition =
                                        "input.IntegrationSelect==
                                    'gene_genomic_res'",
                                    selectInput(
                                        inputId = ns3("genomicTypeSelect"),
                                        label = "Genomic Type:",
                                        choices = intersect(unique(
                                          data_table$cnv_met), c("met", "cnv"))
                                    ), ns = ns3
                                ),
                                selectInput(
                                    inputId = ns3("ClassSelect"),
                                    label = "Class:",
                                    choices = unique(data_table$class)
                                ),
                                selectInput(
                                    inputId = ns3("SignificativityCriteria"),
                                    label = "Significativity criteria:",
                                    choices = c("pval", "FDR")
                                ),
                                conditionalPanel(
                                    condition =
                                      "input.SignificativityCriteria=='pval'",
                                    sliderInput(ns3("PvalRange"),
                                        "P-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns3
                                ),
                                conditionalPanel(
                                    condition =
                                      "input.SignificativityCriteria=='FDR'",
                                    sliderInput(ns3("FdrRange"),
                                        "FDR-Value Range:",
                                        min = 0,
                                        max = 1,
                                        value = c(0, 0.05),
                                        step = 0.005
                                    ), ns = ns3
                                )
                            ),
                            mainPanel(
                                plotOutput(ns3("plotly"), width = 700),
                                tags$div(
                                    style = "overflow-x: auto;",
                                    dataTableOutput(ns3("table")),
                                    downloadButton(ns3("download_csv"),
                                                   "Download CSV")
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}

# Heatmap Subitem for class comparison
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS radioButtons tags
#' fluidRow mainPanel sidebarPanel
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
.gint_subItem_HeatmapDEGs <- function(data_table) {
    ns <- NS("heat_deg")
    tabItem(
        tabName = "heatmapDEGs",
        fluidRow(
            .gint_heatmap_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("IntegrationSelect"),
                            label = "Integration Type:",
                            choices = .change_int_names(intersect(
                                c(
                                    "gene_genomic_res", "gene_met_res",
                                    "gene_cnv_res", "mirna_cnv_res"
                                ),
                                unique(data_table$omics)
                            ))
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_genomic_res'",
                            sliderInput(ns("numTopGenesHeatmapCNV"),
                                "Number of top genes (CNV):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ),
                            sliderInput(ns("numTopGenesHeatmapMET"),
                                "Number of top genes (MET):",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_cnv_res'",
                            sliderInput(ns("numTopGenesHeatmapCNVonly"),
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_met_res'",
                            sliderInput(ns("numTopGenesHeatmapMETonly"),
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'mirna_cnv_res'",
                            sliderInput(ns("numTopGenesHeatmapmirna_cnv"),
                                "Number of top genes:",
                                value = 10,
                                min = 1,
                                max = 200,
                                step = 10
                            ), ns = ns
                        ),
                        sliderInput(ns("numSamples"),
                            "Number of Samples:",
                            value = 50,
                            min = 1,
                            max = 1000,
                            step = 10
                        ),
                        selectInput(ns("ClassSelect"),
                            "Select the class:",
                            choices = unique(data_table$class),
                            multiple = FALSE
                        ),
                        selectInput(
                            inputId = ns("SignificativityCriteria"),
                            label = "Significativity criteria:",
                            choices = c("pval", "FDR")
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='pval'",
                            sliderInput(ns("PvalRange"),
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = 0.05,
                                step = 0.005
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='FDR'",
                            sliderInput(ns("FdrRange"),
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = 0.05,
                                step = 0.005
                            ), ns = ns
                        ),
                        radioButtons(ns("scaleHeatmap"),
                            "scale by:",
                            choices = c("row", "col", "none")
                        )
                    ),
                    mainPanel(
                        InteractiveComplexHeatmapOutput(ns("heatmap"))
                    )
                )
            )
        )
    )
}

# Chromosome Distribution Subitem for class comparison
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput downloadButton tabPanel NS tags fluidRow
#'  mainPanel sidebarPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom gtools mixedsort
#' @importFrom plotly plotlyOutput
#' @importFrom DT dataTableOutput
#' @importFrom stats na.omit
.gint_subItem_chrDistribDEGs <- function(data_table) {
    ns <- NS("histo_deg")
    chr <- c("All", mixedsort(unique(na.omit(data_table$chr_cov))))
    tabItem(
        tabName = "histoDEGs",
        fluidRow(
            .gint_chrDistrib_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("IntegrationSelect"),
                            label = "Integration Type:",
                            choices = .change_int_names(
                              unique(data_table$omics))
                        ),
                        conditionalPanel(
                            condition = "input.IntegrationSelect==
                            'gene_genomic_res'",
                            selectInput(
                                inputId = ns("genomicTypeSelect"),
                                label = "Type Selection",
                                choices = intersect(unique(data_table$cnv_met),
                                                    c("met", "cnv"))
                            ), ns = ns
                        ),
                        selectInput(
                            inputId = ns("ClassSelect"),
                            label = "Class:",
                            choices = unique(data_table$class)
                        ),
                        selectInput(
                            inputId = ns("ChrSelect"),
                            label = "Chr:",
                            choices = chr
                        ),
                        selectInput(
                            inputId = ns("SignificativityCriteria"),
                            label = "Significativity criteria:",
                            choices = c("pval", "FDR")
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='pval'",
                            sliderInput(ns("PvalRange"),
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='FDR'",
                            sliderInput(ns("FdrRange"),
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005
                            ), ns = ns
                        )
                    ),
                    mainPanel(
                        plotlyOutput(ns("plotly")),
                        tags$div(
                            style = "overflow-x: auto;",
                            dataTableOutput(ns("table")),
                            downloadButton(ns("download_csv"), "Download CSV")
                        ),
                        plotlyOutput(ns("plotly_tf")),
                        dataTableOutput(ns("table_tf")),
                        downloadButton(ns("download_csv_tf"), "Download CSV")
                    )
                )
            )
        )
    )
}

# Network Subitem for class comparison
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#'  conditionalPanel sliderInput downloadButton tabPanel NS checkboxInput tags
#'  fluidRow mainPanel sidebarPanel div inputPanel
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom visNetwork visNetworkOutput

.gint_subItem_networkDEGs <- function(data_table) {
    ns <- NS("network_deg")
    tabItem(
        tabName = "networkDEGs",
        fluidRow(
            .gint_network_box(),
            mainPanel(
                sidebarLayout(
                    div(
                        style = "width: 1600px;",
                        inputPanel(
                            selectInput(ns("ClassSelect"),
                                label = "Select the Class:",
                                choices = unique(data_table$class)
                            ),
                            selectInput(
                                inputId = ns("SignificativityCriteria"),
                                label = "Significativity criteria:",
                                choices = c("pval", "FDR")
                            ),
                            sliderInput(ns("numInteractions"),
                                label = "Number of Interactions",
                                min = 10,
                                max = nrow(data_table),
                                value = 200,
                                step = 50
                            ),
                            conditionalPanel(
                                condition = "input.SignificativityCriteria==
                                'pval'",
                                sliderInput(ns("PvalRange"),
                                    "P-Value:",
                                    min = 0,
                                    max = 1,
                                    value = c(0.05),
                                    step = 0.005
                                ), ns = ns
                            ),
                            conditionalPanel(
                                condition = "input.SignificativityCriteria==
                                'FDR'",
                                sliderInput(ns("FdrRange"),
                                    "FDR:",
                                    min = 0,
                                    max = 1,
                                    value = c(0.05),
                                    step = 0.005
                                ), ns = ns
                            ),
                            checkboxInput(ns("layout"),
                                label = "Switch to tree Layout:",
                                value = FALSE
                            ),
                            checkboxInput(ns("physics"),
                                label = "Physics",
                                value = FALSE
                            )
                        )
                    ),
                    mainPanel(
                        visNetworkOutput(ns("networkPlot"),
                            height = 800,
                            width = 1600
                        )
                    )
                )
            )
        )
    )
}

# Circos Subitem for Complete Integration
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput icon
#' conditionalPanel sliderInput tabPanel NS tags fluidPage column mainPanel
#' sidebarPanel div actionButton
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom gtools mixedsort
#' @importFrom shiny.gosling goslingOutput use_gosling
.gint_subItem_circosCompleteInt <- function(data_table) {
    tmp <- ifelse(sum(unique(data_table$omics) %in% c(
        "gene_genomic_res",
        "gene_cnv_res",
        "gene_met_res"
    )) > 0,
    "Gene", NA
    )
    tmp2 <- ifelse(sum(unique(data_table$omics) %in% c("mirna_cnv_res")),
        "miRNA", NA
    )
    tmp <- c(tmp, tmp2)
    tmp <- tmp[!is.na(tmp)]
    ns <- NS("circos")
    chr <- unique(data_table$chr_response)
    chr <- c("All", paste0("chr", mixedsort(chr[!is.na(chr)])))
    tabItem(
        tabName = "circosIntegration",
        fluidPage(
            use_gosling(),
            sidebarLayout(
                sidebarPanel(
                    .gint_circos_box(),
                    HTML("<br> <br> <br> <br>"),
                    selectInput(
                        inputId = ns("circosType"),
                        choices = tmp,
                        label = "Gene/miRNA"
                    ),
                    selectInput(
                        inputId = ns("ClassSelect"),
                        choices = unique(data_table$class),
                        label = "Class"
                    ),
                    selectInput(
                        inputId = ns("ChrSelect"),
                        choices = chr,
                        label = "Chromosome"
                    ),
                    selectInput(
                        inputId = ns("layout"),
                        choices = c("circular", "linear"),
                        label = "Layout"
                    ),
                    actionButton(ns("Download_png"), "PNG",
                        icon = icon("cloud-arrow-down")
                    ),
                    actionButton(ns("Download_pdf"), "PDF",
                        icon = icon("cloud-arrow-down")
                    )
                ),
                mainPanel(
                    column(
                        12, goslingOutput(ns("gosling_plot")),
                        div(style = "height: 200px;")
                    )
                )
            )
        )
    )
}

# Table Subitem with all integrations
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#' conditionalPanel sliderInput tabPanel NS downloadButton tags fluidRow
#' mainPanel sidebarPanel div
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' @importFrom DT dataTableOutput
#' @importFrom gtools mixedsort
#' @importFrom stats na.omit
.gint_subItem_tableCompleteInt <- function(data_table) {
    ns <- NS("complete_table")
    chr <- c("All", mixedsort(unique(na.omit(data_table$chr_cov))))
    tabItem(
        tabName = "fullTable",
        fluidRow(
            .gint_table_box(),
            mainPanel(
                sidebarLayout(
                    sidebarPanel(
                        selectInput(
                            inputId = ns("IntegrationSelect"),
                            label = "Integration Type:",
                            choices = .change_int_names(unique(
                              data_table$omics))
                        ),
                        selectInput(
                            inputId = ns("ClassSelect"),
                            label = "Class:",
                            choices = unique(data_table$class)
                        ),
                        selectInput(
                            inputId = ns("ChrSelect"),
                            label = "Chr:",
                            choices = chr
                        ),
                        selectInput(
                            inputId = ns("degSelect"),
                            label = "DEGs:",
                            choices = c("All", "Only DEGs")
                        ),
                        selectInput(
                            inputId = ns("SignificativityCriteria"),
                            label = "Significativity criteria:",
                            choices = c("pval", "FDR")
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='pval'",
                            sliderInput(ns("PvalRange"),
                                "P-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005
                            ), ns = ns
                        ),
                        conditionalPanel(
                            condition = "input.SignificativityCriteria=='FDR'",
                            sliderInput(ns("FdrRange"),
                                "FDR-Value Range:",
                                min = 0,
                                max = 1,
                                value = c(0, 0.05),
                                step = 0.005
                            ), ns = ns
                        )
                    ),
                    mainPanel(
                        tags$div(
                            style = "overflow-x: auto;",
                            dataTableOutput(ns("table")),
                            downloadButton(ns("download_csv"), "Download CSV")
                        )
                    )
                )
            )
        )
    )
}

# Create user interface
#' @importFrom shiny fluidPage sidebarLayout tabsetPanel selectInput
#'  conditionalPanel sliderInput tabPanel NS downloadButton addResourcePath
#'  tags HTML a img span sidebarPanel
#' @importFrom shinydashboard dashboardSidebar sidebarMenu tabItem
#' dashboardHeader dashboardBody tabItems dashboardPage
.create_ui <- function(data_table) {
    myImgResources <- "imgResources/logo_gINTomics2.png"
    addResourcePath(prefix = "imgResources", system.file(
      directoryPath = "www/", package = "gINTomics"))
    dashboardPage(
        dashboardHeader(
            title = span(
                "gINTomics",
                span("Visualizer 1.0",
                    style = "color: gray; font-size: 16px"
                )
            ),
            tags$li(
                a(
                    href = "https://github.com/angelovelle96/gINTomics",
                    img(src = myImgResources, height = 80, width = 220),
                    style = "padding-top:0px; padding-bottom:0px"
                ),
                class = "dropdown"
            ),
            tags$li(
                class = "dropdown",
                tags$style(".main-header {max-height: 80px}"),
                tags$style(".main-header .logo {height: 80px;
                                       line-height: 80px !important;
                                       padding: 0 0px;}"),
                tags$style(".sidebar-toggle {height: 80px;}")
            )
        ),
        .gint_dashboardsidebar(),
        dashboardBody(
            tags$head(tags$style(
                HTML(".content-wrapper,
                     .right-side {background-color: #ffffff;}")
            )),
            tabItems(
                .gint_tabitem_home(data_table),
                .gint_subItem_coefDistribGenomic(data_table),
                .gint_subItem_HeatmapGenomic(data_table),
                .gint_subItem_chrDistribGenomic(data_table),
                .gint_tabitem_enr(data_table),
                .gint_subItem_coefDistribTranscript(data_table),
                .gint_subItem_chrDistribTranscript(data_table),
                .gint_subItem_networkTranscript(data_table),
                .gint_subItem_enrichTranscript(data_table),
                .gint_subItem_coefDistribDEGs(data_table),
                .gint_subItem_HeatmapDEGs(data_table),
                .gint_subItem_chrDistribDEGs(data_table),
                .gint_subItem_networkDEGs(data_table),
                .gint_subItem_circosCompleteInt(data_table),
                .gint_subItem_tableCompleteInt(data_table)
            )
        ),
        skin = "purple"
    )
}

# Defining home box description
#' @importFrom shiny addResourcePath fluidPage tags img fluidRow column
#' @importFrom shinydashboard tabItem box
.gint_home_box <- function() {
    fluidRow(box(
        title = "Contacts and links",
        footer = HTML("
- Angelo Velle <br>
e-mail: angelo.velle@unipd.it <br>
github: https://github.com/angelovelle96 <br> <br>

- Francesco Patane' <br>
e-mail: francesco.patane@unipd.it <br>
github: https://github.com/francescopatane96 <br> <br>

- Chiara Romualdi <br>
e-mail: chiara.romualdi@unipd.it <br>
website: https://romualdi.bio.unipd.it/
      ")
    ))
}

# Defining coefficients distribution box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_coefDistrib_box <- function() {
    box(
        title = "Coefficients Distribution",
        collapsible = TRUE,
        collapsed = FALSE,
        footer = HTML("
                             The coefficients distribution subpanel contains
                             three plot types that are Venn Diagram, Volcano
                             Plot, and RidgeLine Plot. <br>
                             They were designed to
                             represent the integration coefficients
                             distribution; In particular, The Venn Diagram is
                             designed for the genomic integration. It can help
                             to identify genes which are significantly
                             regulated by both CNV and methylation. <br>
                             The Volcano Plot shows the distribution of
                             integration coefficients for every integration
                             type. For each integration coefficient, on the y
                             axis you have the -log10 of Pvalue/FDR and on the
                             x axis the value of the coefficient. <br>
                             Finally, the ridgeline plot is designed to
                             compare different distributions, it has been
                             integrated in the package with the aim to compare
                             the distribution of significant and non
                             significant coefficients returned by our
                             integration models. For each distribution, on the
                             y axis you have the frequencies and on the x axis
                             the values of the coefficients.
      ")
    )
}

# Defining Heatmap box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_heatmap_box <- function() {
    box(
        title = "Heatmap",
        collapsible = TRUE,
        collapsed = FALSE,
        footer = HTML("
                    Heatmap subpanel contains an interactive heatmap developed
                    using InteractiveComplexHeatmap package, and it shows the
                    value of genomic integration coefficients (only cnv or met,
                    or both of them) in relation of gene expression and along
                    the different input samples. What is more, the user can
                    scale values by rows or columns, or simply view raw values.
      ")
    )
}

# Defining Chromosome distribution box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_chrDistrib_box <- function() {
    box(
        title = "Chromosome Distribution",
        collapsible = TRUE,
        collapsed = FALSE,
        footer = HTML("
                    Chromosome distribution subpanel allows the user to
                    visualize the distribution of integration coefficients
                    along the different chromosomes and to focus on specific
                    chromosomes and genes trough an interactive and
                    downloadable data table. This kind of visualization will
                    identify chromosomes in which the type of regulation under
                    analysis is particularly active.
      ")
    )
}

# Defining Network box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_network_box <- function() {
    box(
        title = "Network Visualization",
        collapsible = TRUE,
        collapsed = TRUE,
        footer = HTML("
                    Network subpanel shows the significant interactions between
                    transcriptional regulators (TFs and miRNAs, if present)
                    and their targets genes/miRNA. Nodes and edges are selected
                    ordering them by the most high coefficient values (absolute
                    value) and by default the top 200 interactions are showed.
      ")
    )
}

# Defining Enrichment box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_enrich_box <- function() {
    box(
        title = "Enrichment Analysis",
        collapsible = TRUE,
        collapsed = FALSE,
        footer = HTML("
                    Enrichment subpanel shows the enrichment results obtained
                    with enrichGO and enrichKEGG (clusterProfiler). The genomic
                    enrichment is performed providing the list of genes
                    significantly regulated by methylation or CNV, while the
                    transcriptional one with the list of genes significantly
                    regulated by each transcription factor (we run an
                    enrichment for each TF that significantly regulates at
                    least 12 targets).
      ")
    )
}

# Defining Circos box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_circos_box <- function() {
    box(
        title = "Circos Plot",
        collapsible = TRUE,
        collapsed = TRUE,
        footer = HTML("
                    Circos plots subpanel lodges Circos plots for genomic
                    integration (it requires the presence of both cnv and met
                    data or only one of them).
                    The strenght of this type of visualization lies in the
                    possibility to visulize different omics together, moreover
                    the user of zoom in the plot focusing the attention to only
                    specific chromosome of interest or particular genomic
                    regions. What is more, the plot can be linearized and by
                    passing the pointer on genes it is possible to visualize a
                    tooltip with interesting informations like the pval,
                    FDR, gene name and symbol etc.. Once the user has choosen
                    the genomic region of interest, it is also possible to
                    dowload an instantaneous of the plot.
      "), width = "100%"
    )
}

# Defining Table box description
#' @importFrom shiny fluidPage tags
#' @importFrom shinydashboard box
.gint_table_box <- function() {
    box(
        title = "Full Integration Data Table",
        collapsible = TRUE,
        collapsed = FALSE,
        footer = HTML("
                   Data Table subpanel contains a table that shows the full and
                   initial data frame containing all data that is used by the
                   shiny app for generating plots and relate tables.
      ")
    )
}
