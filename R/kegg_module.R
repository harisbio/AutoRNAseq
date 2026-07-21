# ============================================================
# KEGG Pathway Analysis Module
# AutoRNAseq v1.0.0
# ============================================================

keggAnalysisUI <- function(id) {
    ns <- NS(id)
    tagList(
        box(
            width = NULL,
            title = span(icon("road"), " KEGG Pathway Analysis"),
            status = "primary",
            solidHeader = TRUE,
            style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
            sidebarLayout(
                sidebarPanel(
                    style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
                    h4("⚙️ Parameters", style = "color: #0092AC;"),
                    selectInput(ns("kegg_organism"), "Organism (KEGG code):",
                        choices = downstream_kegg_choices(),
                        selected = "hsa"
                    ),
                    selectInput(ns("gene_id_type"), "Gene ID Type:",
                        choices = c("SYMBOL", "ENSEMBL", "ENTREZID"), selected = "SYMBOL"
                    ),
                    selectInput(ns("org_db"), "OrgDb Package:",
                        choices = downstream_orgdb_choices(include_none = TRUE),
                        selected = "org.Hs.eg.db"
                    ),
                    numericInput(ns("pval_cutoff"), "p-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
                    numericInput(ns("lfc_cutoff"), "Min |log2FC|:", value = 1.0, min = 0, max = 10, step = 0.1),
                    hr(),
                    h4("📂 DEG Results Input", style = "color: #0092AC;"),
                    fileInput(ns("deg_file"), "Upload DEG CSV (.csv)"),
                    hr(),
                    div(
                        style = "display:flex; gap:10px;",
                        actionButton(ns("btn_run_kegg"), "Run KEGG Analysis",
                            icon = icon("play"),
                            style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                ),
                mainPanel(
                    tabBox(
                        width = 12,
                        tabPanel(
                            title = span(icon("table"), " Enrichment Table"),
                            value = "kegg-table",
                            box(
                                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                style = "background-color:#fff; border-radius:10px; padding:15px;",
                                DT::dataTableOutput(ns("kegg_table"))
                            ),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_kegg_csv"), "Download CSV",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("chart-bar"), " Dotplot"),
                            value = "kegg-dot",
                            plotOutput(ns("kegg_dotplot"), height = "500px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_kegg_dot"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("network-wired"), " Enrichment Map"),
                            value = "kegg-emap",
                            plotOutput(ns("kegg_emap"), height = "600px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_kegg_emap"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        )
                    ),
                    br(),
                    div(
                        style = "display:flex; justify-content:left;",
                        actionButton(ns("help_kegg"), "Help",
                            icon = icon("info-circle"),
                            style = "color:#fff; background-color:#F39C12; border-color:#E67E22; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                )
            )
        )
    )
}


keggAnalysisServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observeEvent(input$kegg_organism, {
            mapped_orgdb <- downstream_orgdb_from_kegg(input$kegg_organism)
            if (is.na(mapped_orgdb) || !nzchar(mapped_orgdb)) {
                mapped_orgdb <- "none"
            }
            updateSelectInput(session, "org_db", selected = mapped_orgdb)
        }, ignoreInit = TRUE)

        observeEvent(input$help_kegg, {
            show_step_help(
              "KEGG Help",
              c(
                "Upload a DEG CSV with a gene column plus log2FoldChange and padj.",
                "Choose the correct KEGG organism code for your species.",
                "If you select ENTREZID, the gene identifiers must already be mapped.",
                "Use a filtered DEG list first to verify the input format.",
                "KEGG results work best when the annotation and organism choice match the biology of your study."
              ),
              "The DEG CSV is usually small enough for a browser upload; the main requirement is correct identifier type and organism mapping."
            )
        })

        # ── DEG loading & gene list ──────────────────────────────────────────────
        deg_data <- reactive({
            req(input$deg_file)
            df <- read.csv(input$deg_file$datapath, header = TRUE, check.names = FALSE)
            colnames(df)[1] <- "gene"
            df
        })

        entrez_ids <- reactive({
            req(deg_data(), input$gene_id_type, input$org_db)
            df <- deg_data()
            sig <- df[!is.na(df$padj) & !is.na(df$log2FoldChange) &
                df$padj < input$pval_cutoff &
                abs(df$log2FoldChange) >= input$lfc_cutoff, ]
            genes <- as.character(sig$gene)
            genes <- genes[!is.na(genes) & genes != ""]

            if (input$gene_id_type == "ENTREZID") {
                return(genes)
            }

            pkg <- input$org_db
            if (pkg == "none") {
                showNotification("This organism has no local OrgDb mapping. Use ENTREZID as Gene ID Type.", type = "warning")
                return(NULL)
            }
            
            if (!requireNamespace(pkg, quietly = TRUE)) {
                showNotification(paste0("OrgDB package '", pkg, "' is not installed."), type = "error")
                return(NULL)
            }
            orgDB <- tryCatch(getExportedValue(pkg, pkg), error = function(e) NULL)
            
            if (is.null(orgDB)) {
                showNotification(paste0("OrgDB package '", pkg, "' could not be loaded."), type = "error")
                return(NULL)
            }


            mapped <- AnnotationDbi::mapIds(orgDB,
                keys = genes,
                column = "ENTREZID", keytype = input$gene_id_type,
                multiVals = "first"
            )
            as.character(mapped[!is.na(mapped)])
        })

        # ── KEGG enrichment ──────────────────────────────────────────────────────
        kegg_result <- reactiveVal(NULL)

        observeEvent(input$btn_run_kegg, {
            ids <- entrez_ids()
            req(ids)
            if (length(ids) == 0) {
                showNotification("No significant genes to analyse.", type = "warning")
                return()
            }
            withProgress(message = "Running KEGG Pathway Analysis...", value = 0.5, {
                tryCatch(
                    {
                        res <- clusterProfiler::enrichKEGG(
                            gene = ids,
                            organism = input$kegg_organism,
                            pvalueCutoff = input$pval_cutoff,
                            pAdjustMethod = "BH"
                        )
                        kegg_result(res)
                        showNotification("KEGG Analysis completed!", type = "message")
                    },
                    error = function(e) {
                        showNotification(paste("Error:", conditionMessage(e)), type = "error")
                    }
                )
            })
        })

        # ── outputs ──────────────────────────────────────────────────────────────
        output$kegg_table <- DT::renderDataTable({
            req(kegg_result())
            DT::datatable(as.data.frame(kegg_result()), options = list(scrollX = TRUE, pageLength = 10))
        })

        output$kegg_dotplot <- renderPlot({
            req(kegg_result())
            enrichplot::dotplot(kegg_result(), showCategory = 20) +
                ggplot2::ggtitle("KEGG Pathway Dotplot")
        })

        output$kegg_emap <- renderPlot({
            req(kegg_result())
            res2 <- enrichplot::pairwise_termsim(kegg_result())
            enrichplot::emapplot(res2) + ggplot2::ggtitle("KEGG Enrichment Map")
        })

        # ── downloads ─────────────────────────────────────────────────────────────
        output$dl_kegg_csv <- downloadHandler(
            filename = function() paste0("KEGG_enrichment_", Sys.Date(), ".csv"),
            content  = function(file) write.csv(as.data.frame(kegg_result()), file, row.names = FALSE)
        )
        output$dl_kegg_dot <- downloadHandler(
            filename = function() paste0("KEGG_dotplot_", Sys.Date(), ".png"),
            content = function(file) {
                png(file, width = 1200, height = 800, res = 120)
                print(enrichplot::dotplot(kegg_result(), showCategory = 20))
                dev.off()
            }
        )
        output$dl_kegg_emap <- downloadHandler(
            filename = function() paste0("KEGG_emapplot_", Sys.Date(), ".png"),
            content = function(file) {
                png(file, width = 1200, height = 800, res = 120)
                res2 <- enrichplot::pairwise_termsim(kegg_result())
                print(enrichplot::emapplot(res2))
                dev.off()
            }
        )
    })
}

