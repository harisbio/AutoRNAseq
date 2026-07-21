## All modules in R/modules/ are loaded automatically by R package infrastructure.## All modules in R/modules/ are loaded automatically by R package infrastructure.


#' AutoRNAseq: Sequence Expression Analyser
#'
#' AutoRNAseq enables comprehensive differential expression analysis
#' of RNA-Seq data, plus downstream functional analysis and biomarker discovery.
#' This function launches the main Shiny web application.
#'
#' @name runAnalyser
#' @return A Shiny App is launched for interactive and comprehensive
#'   differential gene expression analysis of RNA-Seq data.
#' @export
#'
#' @examples
#' \dontrun{
#' runAnalyser()
#' }
runAnalyser <- function(host = "127.0.0.1", port = NULL) {
  # ── Global Options ──────────────────────────────────────────────────────────
  max_request_size <- as.numeric(Sys.getenv("SHINY_MAX_REQUEST_SIZE", unset = as.character(12288 * 1024^2)))
  if (is.na(max_request_size) || max_request_size <= 0) {
    max_request_size <- 12288 * 1024^2
  }
  options(shiny.maxRequestSize = max_request_size) # Bytes; can be overridden by SHINY_MAX_REQUEST_SIZE
  message(sprintf("AutoRNAseq upload limit (shiny.maxRequestSize): %.2f GB", max_request_size / 1024^3))
  options(shiny.launch.browser = FALSE)

  # ── UI ──────────────────────────────────────────────────────────────────────
  ui <- shinydashboard::dashboardPage(
    skin = "black",
    title = "AutoRNAseq",

    # Header ─────────────────────────────────────────────────────────────────
    shinydashboard::dashboardHeader(
      title = tags$span(
        img(src = "AutoRNAseq/app_logo.png?v=20260320", height = "40px"),
        "AutoRNAseq"
      ),
      titleWidth = 380
    ),

    # Sidebar ─────────────────────────────────────────────────────────────────
    shinydashboard::dashboardSidebar(
      width = 350,
      tags$div(
        class = "arnx-sidebar-top",
        tags$div(
          class = "arnx-sidebar-brand",
          tags$img(src = "AutoRNAseq/app_logo.png?v=20260320", alt = "AutoRNAseq logo"),
          tags$div(
            class = "arnx-sidebar-brand-text",
            tags$div(class = "brand-title", "AutoRNAseq"),
            tags$div(class = "brand-subtitle", "AI-Driven RNA-Seq Intelligence")
          )
        ),
        tags$div(
          class = "arnx-sidebar-pills",
          tags$span("End-to-End"),
          tags$span("Clinical-Ready"),
          tags$span("Interactive")
        )
      ),
      shinydashboard::sidebarMenu(
        # ── Existing pipeline ──────────────────────────────────────────────
        shinydashboard::menuItem("Welcome !",
          tabName = "welcomePage",
          icon = shiny::icon("home")
        ),
        shinydashboard::menuItem("Guidebook",
          tabName = "guidebookPage",
          icon = shiny::icon("book")
        ),
        shinydashboard::menuItem("Automated Pipeline",
          tabName = "automated_pipeline",
          icon = shiny::icon("robot")
        ),
        shinydashboard::menuItem("Data Setup",
          tabName = "collected_data",
          icon = shiny::icon("download")
        ),
        shinydashboard::menuItem("Quality Control",
          tabName = "quality",
          icon = shiny::icon("check-circle")
        ),
        shinydashboard::menuItem("Filtering and Trimming",
          tabName = "filtering",
          icon = shiny::icon("filter")
        ),
        shinydashboard::menuItem("Alignment",
          tabName = "alignment",
          icon = shiny::icon("align-left")
        ),
        shinydashboard::menuItem("Quantification",
          tabName = "quantification",
          icon = shiny::icon("list-ol")
        ),
        shinydashboard::menuItem("Differential Expression",
          tabName = "DESeq2",
          icon = shiny::icon("balance-scale")
        ),

        # ── Divider ────────────────────────────────────────────────────────
        tags$hr(style = "border-top: 1px solid #555; margin: 10px 15px;"),

        # ── New modules ────────────────────────────────────────────────────
        shinydashboard::menuItem("Downstream Functional Analysis",
          icon = shiny::icon("search"),
          startExpanded = FALSE,
          shinydashboard::menuSubItem("GO Enrichment", tabName = "goEnrichment", icon = shiny::icon("circle")),
          shinydashboard::menuSubItem("KEGG Pathways", tabName = "keggAnalysis", icon = shiny::icon("road")),
          shinydashboard::menuSubItem("GSEA", tabName = "gseaAnalysis", icon = shiny::icon("chart-line"))
        ),
        shinydashboard::menuItem("PPI Network Analysis",
          tabName = "ppiNetwork",
          icon = shiny::icon("project-diagram")
        ),
        shinydashboard::menuItem("Biomarker Discovery (ML)",
          tabName = "mlBiomarker",
          icon = shiny::icon("microchip")
        ),
        # shinydashboard::menuItem("Reports",
        #   tabName = "reportsPage",
        #   icon = shiny::icon("file-alt")
        # ),

        # ── About ──────────────────────────────────────────────────────────
        tags$hr(style = "border-top: 1px solid #555; margin: 10px 15px;"),
        shinydashboard::menuItem("About",
          tabName = "aboutPage",
          icon = shiny::icon("university")
        )
      )
    ),

    # Body ───────────────────────────────────────────────────────────────────
    shinydashboard::dashboardBody(
      class = "autornaseq-ui",
      introjsUI(),
      tags$head(
        tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
        tags$link(rel = "icon", type = "image/png", href = "AutoRNAseq/favicon.png?v=20260419c"),
        tags$link(rel = "stylesheet", href = "AutoRNAseq/autornaseq-theme.css?v=20260419c"),
        tags$script(src = "AutoRNAseq/autornaseq-ui.js?v=20260419c"),
        tags$style(HTML("
          #myScrollBox { overflow-y: scroll; }
          .sidebar-menu .header { color: #ffffff !important; font-size: 11px; letter-spacing: 0.06em; padding: 8px 18px; text-transform: uppercase; }
        "))
      ),
      tabItems(
        # ── Existing tabs ──────────────────────────────────────────────────
        tabItem(
          tabName = "welcomePage",
          box(
            width = FALSE, solidHeader = FALSE,
            includeMarkdown(system.file("extdata", "welcome.md",
              package = "AutoRNAseq"
            ))
          )
        ),
        tabItem(
          tabName = "guidebookPage",
          box(
            width = FALSE, solidHeader = FALSE,
            includeMarkdown(system.file("extdata", "guidebook.md",
              package = "AutoRNAseq"
            ))
          )
        ),
        tabItem(tabName = "automated_pipeline", automatedPipelineUI("autoPipeline")),
        tabItem(tabName = "collected_data", importDataUI("importData")),
        tabItem(tabName = "quality", qualityControlUI("qualityControl")),
        tabItem(tabName = "filtering", Filter_TrimUI("filter_trim")),
        tabItem(tabName = "alignment", alignmentUI("alignment")),
        tabItem(tabName = "quantification", quantificationUI("quantification")),
        tabItem(tabName = "DESeq2", dif_gene_ex_analysisUI("DGEA")),

        # ── New tabs ───────────────────────────────────────────────────────
        tabItem(tabName = "goEnrichment", goEnrichmentUI("goEnrich")),
        tabItem(tabName = "keggAnalysis", keggAnalysisUI("keggAnalysis")),
        tabItem(tabName = "gseaAnalysis", gseaAnalysisUI("gseaAnalysis")),
        tabItem(tabName = "ppiNetwork", ppiNetworkUI("ppiNetwork")),
        tabItem(tabName = "mlBiomarker", mlBiomarkerUI("mlBiomarker")),

        # # ── Reports tab ────────────────────────────────────────────────────
        # tabItem(
        #   tabName = "reportsPage",
        #   box(
        #     width = NULL, title = span(icon("file-lines"), " Generate Analysis Report"),
        #     status = "primary", solidHeader = TRUE,
        #     style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
        #     sidebarLayout(
        #       sidebarPanel(
        #         style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
        #         h4("📄 Report Configuration", style = "color: #0092AC;"),
        #         helpText("Upload the results files from previous analysis steps to include in the report."),
        #         fileInput("report_deg_csv", "DEG Results CSV"),
        #         fileInput("report_go_csv", "GO Enrichment CSV"),
        #         fileInput("report_kegg_csv", "KEGG Enrichment CSV"),
        #         fileInput("report_hub_csv", "PPI Hub Genes CSV"),
        #         fileInput("report_biomarker_csv", "Biomarker Ranking CSV"),
        #         hr(),
        #         selectInput("report_format", "Output Format:",
        #           choices = c(
        #             "HTML" = "html_document",
        #             "PDF" = "pdf_document",
        #             "Word" = "word_document"
        #           ),
        #           selected = "html_document"
        #         ),
        #         hr(),
        #         div(
        #           style = "display:flex; gap:10px;",
        #           downloadButton("btn_generate_report", "Generate Report",
        #             style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:8px 16px; font-size:14px; border-radius:5px;"
        #           )
        #         )
        #       ),
        #       mainPanel(
        #         box(
        #           width = NULL, status = "success", solidHeader = TRUE,
        #           title = "Report Contents",
        #           style = "background-color:#fff; border-radius:10px; padding:15px;",
        #           tags$ul(
        #             tags$li(icon("check"), " Differential Expression Summary"),
        #             tags$li(icon("check"), " GO Enrichment Results + Plots"),
        #             tags$li(icon("check"), " KEGG Pathway Results + Plots"),
        #             tags$li(icon("check"), " GSEA Results"),
        #             tags$li(icon("check"), " PPI Hub Genes"),
        #             tags$li(icon("check"), " ML Biomarker Candidates"),
        #             tags$li(icon("check"), " Session Information")
        #           ),
        #           tags$hr(),
        #           tags$p(
        #             style = "color:#555;",
        #             "The report is generated automatically from the CSV files you upload above. ",
        #             "Formats available: HTML (interactive), PDF, Word."
        #           )
        #         )
        #       )
        #     )
        #   )
        # ),

        # ── About tab ──────────────────────────────────────────────────────
        tabItem(
          tabName = "aboutPage",
          box(
            width = FALSE, solidHeader = FALSE,
            includeMarkdown(system.file("extdata", "about.md",
              package = "AutoRNAseq"
            )),
            verbatimTextOutput("sessioninfo")
          )
        )
      ),
      footer()
    )
  )

  # ── Server ──────────────────────────────────────────────────────────────────
  server <- function(input, output, session) {
    # ── Existing modules ───────────────────────────────────────────────────
    automatedPipelineServer("autoPipeline")
    importedData <- callModule(importDataServer, "importData")
    qualityControlServer("qualityControl", importedData)
    Filter_TrimServer("filter_trim", importedData)
    alignmentServer("alignment")
    quantificationServer("quantification")
    dif_gene_ex_analysisServer("DGEA")

    # ── New modules ────────────────────────────────────────────────────────
    goEnrichmentServer("goEnrich")
    keggAnalysisServer("keggAnalysis")
    gseaAnalysisServer("gseaAnalysis")
    ppiNetworkServer("ppiNetwork")
    mlBiomarkerServer("mlBiomarker")

    # ── Session info ───────────────────────────────────────────────────────
    output$sessioninfo <- renderPrint({
      sessionInfo()
    })

    # ── Report generation ──────────────────────────────────────────────────
    output$btn_generate_report <- downloadHandler(
      filename = function() {
        ext <- switch(input$report_format,
          "html_document" = "html",
          "pdf_document"  = "pdf",
          "word_document" = "docx"
        )
        paste0("AutoRNAseq_Report_", Sys.Date(), ".", ext)
      },
      content = function(file) {
        template <- system.file("report_templates", "report.Rmd",
          package = "AutoRNAseq"
        )

        # Load CSVs if uploaded
        load_csv <- function(fileInput) {
          if (!is.null(fileInput) && !is.null(fileInput$datapath)) {
            tryCatch(read.csv(fileInput$datapath, header = TRUE), error = function(e) NULL)
          } else {
            NULL
          }
        }

        params_list <- list(
          deg_table   = load_csv(input$report_deg_csv),
          go_result   = NULL, # CSVs are pre-processed data frames
          kegg_result = NULL,
          gsea_result = NULL,
          hub_genes   = load_csv(input$report_hub_csv),
          biomarkers  = load_csv(input$report_biomarker_csv),
          ml_metrics  = NULL
        )

        withProgress(message = "Generating report...", value = 0.5, {
          tryCatch(
            {
              rmarkdown::render(
                input = template,
                output_format = input$report_format,
                output_file = file,
                params = params_list,
                envir = new.env(parent = globalenv())
              )
            },
            error = function(e) {
              showNotification(paste("Report error:", conditionMessage(e)), type = "error")
            }
          )
        })
      }
    )
  }

  # ── Launch app ──────────────────────────────────────────────────────────────
  runApp(list(ui = ui, server = server), host = host, port = port)
}

