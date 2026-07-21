footer <- function() {
  tags$div(
    class = "panel-footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      tags$div(
        tags$div("© 2026 AutoRNAseq — An AI-powered, fully automated RNA-Seq analysis pipeline for comprehensive transcriptomic insights."),
        tags$div(HTML('Developed by <a class="portfolio-link" href="https://harisbio.github.io/" target="_blank" rel="noopener noreferrer">Muhammad Haris</a> | Department of Bioinformatics & Biotechnology, Government College Faisalabad, Pakistan'))
      )
    )
  )
}

help_modal_button <- function(input_id, label = "Help") {
  shiny::actionButton(
    input_id,
    label,
    icon = shiny::icon("circle-info"),
    style = "color:#fff; background-color:#F39C12; border-color:#E67E22; padding:6px 12px; font-size:14px; border-radius:5px;"
  )
}

show_step_help <- function(title, bullets, extra = NULL) {
  items <- paste0("<li>", bullets, "</li>", collapse = "")
  body <- paste0("<div class='autornaseq-help'><ul>", items, "</ul>")
  if (!is.null(extra) && nzchar(extra)) {
    body <- paste0(body, "<p>", extra, "</p>")
  }
  body <- paste0(body, "</div>")
  shiny::showModal(shiny::modalDialog(
    title = title,
    shiny::HTML(body),
    easyClose = TRUE,
    footer = shiny::modalButton("Close"),
    size = "m"
  ))
}

downstream_organisms <- function() {
  data.frame(
    label = c(
      "Human (Homo sapiens)",
      "Mouse (Mus musculus)",
      "Rat (Rattus norvegicus)",
      "Fruit Fly (Drosophila melanogaster)",
      "Arabidopsis (Arabidopsis thaliana)",
      "Maize (Zea mays)",
      "Rice (Oryza sativa)",
      "Wheat (Triticum aestivum)",
      "Sugarcane (Saccharum officinarum)",
      "Soybean (Glycine max)",
      "Tomato (Solanum lycopersicum)",
      "Barley (Hordeum vulgare)",
      "Sorghum (Sorghum bicolor)"
    ),
    kegg = c("hsa", "mmu", "rno", "dme", "ath", "zma", "osa", "tae", "scc", "gmx", "sly", "hvu", "sbi"),
    orgdb = c(
      "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db",
      "org.At.tair.db", "org.Zmays.eg.db", "org.Osativa.eg.db",
      NA, NA, NA, NA, NA, NA
    ),
    ppi_taxid = c(9606, 10090, 10116, 7227, 3702, 4577, 4530, 4565, 4547, 3847, 4081, 4513, 4558),
    stringsAsFactors = FALSE
  )
}

downstream_kegg_choices <- function() {
  org <- downstream_organisms()
  stats::setNames(org$kegg, org$label)
}

downstream_orgdb_choices <- function(include_none = FALSE) {
  org <- downstream_organisms()
  has_orgdb <- !is.na(org$orgdb) & nzchar(org$orgdb)
  choices <- stats::setNames(org$orgdb[has_orgdb], org$label[has_orgdb])
  if (isTRUE(include_none)) {
    choices <- c(choices, "Not available / Use ENTREZ IDs" = "none")
  }
  choices
}

downstream_orgdb_from_kegg <- function(kegg_code) {
  org <- downstream_organisms()
  idx <- which(org$kegg == kegg_code)
  if (length(idx) == 0) {
    return(NA_character_)
  }
  org$orgdb[idx[1]]
}

downstream_ppi_choices <- function() {
  org <- downstream_organisms()
  has_ppi <- !is.na(org$ppi_taxid)
  labels <- paste0(org$label[has_ppi], " (", org$ppi_taxid[has_ppi], ")")
  stats::setNames(as.character(org$ppi_taxid[has_ppi]), labels)
}

select_files <- function(session, inputID, inputData, message) {
  fileChoices <- reactive({
    if (is.null(inputData())) {
      return(c())
      print(inputID)
    }
    files <- names(inputData())
    choices <- c()
    for (filename in files) {
      choices <- c(choices, filename)
    }
    return(choices)
  })
  observe({
    updateSelectInput(session, inputID, label = message, choices = fileChoices())
  })
}

blank_plot <- function(main) {
  plot(
    x = 0, y = 0, xlim = c(0, 100), ylim = c(0, 100), type = "n",
    xlab = "X-axis", ylab = "Y-axis",
    main = "Plot Will Appear Here Once Data is Selected"
  )
}

quality_plots <- function(selectedData, selectedFiles, importedData) {
  plot_functions <- list(
    "Average Quality" = rqcReadQualityPlot,
    "Cycle-specific Average Quality" = rqcCycleAverageQualityPlot,
    "Cycle-specific Base Call Proportion" = rqcCycleBaseCallsLinePlot,
    "Cycle-specific GC Content" = rqcCycleGCPlot,
    "Cycle-specific Quality Distribution" = rqcCycleQualityPlot,
    "Cycle-specific Quality Distribution - Boxplot" = rqcCycleQualityBoxPlot,
    "Per Read Mean Quality Distribution of Files" = rqcReadQualityBoxPlot,
    "Read Frequency" = rqcReadFrequencyPlot,
    "Read Length Distribution" = rqcReadWidthPlot
  )
  plot_functions[[selectedData]](importedData()[selectedFiles])
}

datatable <- function(selectedData, selectedFiles, importedData) {
  data_functions <- list(
    "Average Quality" = rqcReadQualityCalc,
    "Cycle-specific Average Quality" = rqcCycleAverageQualityCalc,
    "Cycle-specific Base Call Proportion" = rqcCycleBaseCallsCalc,
    "Cycle-specific GC Content" = rqcCycleGCCalc,
    "Cycle-specific Quality Distribution" = rqcCycleQualityCalc,
    "Cycle-specific Quality Distribution - Boxplot" = rqcCycleQualityBoxCalc,
    "Per Read Mean Quality Distribution of Files" = rqcReadQualityBoxCalc,
    "Read Frequency" = rqcReadFrequencyCalc,
    "Read Length Distribution" = rqcReadWidthCalc
  )
  resultList <- lapply(selectedFiles, function(file) {
    data_functions[[selectedData]](importedData()[file])
  })
  do.call(rbind, resultList)
}
