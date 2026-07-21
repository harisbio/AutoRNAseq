# ============================================================
# GSEA Module
# AutoRNAseq v2
# ============================================================

gseaAnalysisUI <- function(id) {
    ns <- NS(id)
    tagList(
        box(
            width = NULL,
            title = span(icon("wave-square"), " Gene Set Enrichment Analysis (GSEA)"),
            status = "primary",
            solidHeader = TRUE,
            style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
            sidebarLayout(
                sidebarPanel(
                    style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
                    h4("⚙️ Parameters", style = "color: #0092AC;"),
                    selectInput(ns("gsea_type"), "GSEA Method:",
                        choices = c("KEGG" = "kegg", "GO Biological Process" = "GO_BP"),
                        selected = "kegg"
                    ),
                    selectInput(ns("gsea_organism"), "Organism:",
                        choices = downstream_kegg_choices(),
                        selected = "hsa"
                    ),
                    selectInput(ns("organ_db"), "OrgDb Package:",
                        choices = downstream_orgdb_choices(include_none = TRUE),
                        selected = "org.Hs.eg.db"
                    ),
                    selectInput(ns("gene_id_type"), "Gene ID Type:",
                        choices = c("SYMBOL", "ENSEMBL", "ENTREZID"), selected = "SYMBOL"
                    ),
                    selectInput(ns("rank_by"), "Rank gene list by:",
                        choices = c(
                            "log2FoldChange" = "log2FoldChange",
                            "stat (Wald stat)" = "stat"
                        ),
                        selected = "log2FoldChange"
                    ),
                    numericInput(ns("pval_cutoff"), "p-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
                    hr(),
                    h4("📂 DEG Results Input", style = "color: #0092AC;"),
                    fileInput(ns("deg_file"), "Upload DEG CSV (.csv)"),
                    hr(),
                    div(
                        style = "display:flex; gap:10px;",
                        actionButton(ns("btn_run_gsea"), "Run GSEA",
                            icon = icon("play"),
                            style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                ),
                mainPanel(
                    tabBox(
                        width = 12,
                        tabPanel(
                            title = span(icon("table"), " GSEA Results Table"),
                            value = "gsea-table",
                            box(
                                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                style = "background-color:#fff; border-radius:10px; padding:15px;",
                                DT::dataTableOutput(ns("gsea_table"))
                            ),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_gsea_csv"), "Download CSV",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("chart-bar"), " GSEA Dotplot"),
                            value = "gsea-dot",
                            plotOutput(ns("gsea_dotplot"), height = "500px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_gsea_dot"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("wave-square"), " Enrichment Score Plot"),
                            value = "gsea-es",
                            uiOutput(ns("gsea_pathway_ui")),
                            plotOutput(ns("gsea_es_plot"), height = "500px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_gsea_es"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("list"), " Leading Edge Genes"),
                            value = "gsea-leading",
                            box(
                                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                style = "background-color:#fff; border-radius:10px; padding:15px;",
                                DT::dataTableOutput(ns("gsea_leading"))
                            )
                        )
                    ),
                    br(),
                    div(
                        style = "display:flex; justify-content:left;",
                        help_modal_button(ns("help_gsea"), "Help")
                    )
                )
            )
        )
    )
}


gseaAnalysisServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observeEvent(input$gsea_organism, {
            mapped_orgdb <- downstream_orgdb_from_kegg(input$gsea_organism)
            if (is.na(mapped_orgdb) || !nzchar(mapped_orgdb)) {
                mapped_orgdb <- "none"
            }
            updateSelectInput(session, "organ_db", selected = mapped_orgdb)
        }, ignoreInit = TRUE)

        observeEvent(input$help_gsea, {
            show_step_help(
                "GSEA Help",
                c(
                    "Upload a DEG CSV with a gene column and a ranking column such as log2FoldChange or stat.",
                    "Choose the GSEA method, organism, gene ID type, and ranking metric before running the analysis.",
                    "For GO-based GSEA, the selected OrgDb package must be installed and should match the organism.",
                    "For KEGG-based GSEA, select the correct organism code so pathway mapping stays biologically consistent.",
                    "Use a filtered or ranked DEG file first so you can confirm the settings on a smaller input set.",
                    "This module is for downstream enrichment analysis, not raw FASTQ or count-matrix input."
                ),
                extra = "If the module fails to load results, double-check that the gene IDs match the selected ID type and that the chosen organism has the required annotation support."
            )
        })

        # ── Load and rank DEG ────────────────────────────────────────────────────
        deg_data <- reactive({
            req(input$deg_file)
            df <- read.csv(input$deg_file$datapath, header = TRUE, check.names = FALSE)
            colnames(df)[1] <- "gene"
            df
        })

        ranked_gene_list <- reactive({
            req(deg_data(), input$rank_by, input$gene_id_type, input$organ_db)
            df <- deg_data()
            col <- input$rank_by
            genes <- as.character(df$gene)
            values <- as.numeric(df[[col]])

            # Map to ENTREZID if needed
            if (input$gene_id_type != "ENTREZID") {
                if (input$organ_db == "none") {
                    showNotification("This organism has no local OrgDb mapping. Use ENTREZID as Gene ID Type.", type = "warning")
                    return(NULL)
                }
                if (!requireNamespace(input$organ_db, quietly = TRUE)) {
                    showNotification(paste0("OrgDB package '", input$organ_db, "' is not installed."), type = "error")
                    return(NULL)
                }
                orgDB <- tryCatch(getExportedValue(input$organ_db, input$organ_db), error = function(e) NULL)
                if (is.null(orgDB)) {
                    showNotification("OrgDB not found.", type = "error")
                    return(NULL)
                }
                mapped <- AnnotationDbi::mapIds(orgDB,
                    keys = genes,
                    column = "ENTREZID", keytype = input$gene_id_type,
                    multiVals = "first"
                )
                names(values) <- as.character(mapped)
            } else {
                names(values) <- genes
            }

            # Remove NAs and sort descending
            values <- values[!is.na(names(values)) & !is.na(values)]
            values <- sort(values, decreasing = TRUE)
            values
        })

        # ── Run GSEA ─────────────────────────────────────────────────────────────
        gsea_result <- reactiveVal(NULL)

        observeEvent(input$btn_run_gsea, {
            ranked <- ranked_gene_list()
            req(ranked)
            withProgress(message = "Running GSEA...", value = 0.5, {
                tryCatch(
                    {
                        res <- if (input$gsea_type == "kegg") {
                            clusterProfiler::gseKEGG(
                                geneList = ranked,
                                organism = input$gsea_organism,
                                pvalueCutoff = input$pval_cutoff,
                                pAdjustMethod = "BH",
                                verbose = FALSE
                            )
                        } else {
                            pkg <- input$organ_db
                            if (!requireNamespace(pkg, quietly = TRUE)) {
                                showNotification(paste0("OrgDB package '", pkg, "' is not installed."), type = "error")
                                return(NULL)
                            }
                            orgDB <- tryCatch(getExportedValue(pkg, pkg), error = function(e) NULL)
                            
                            if (is.null(orgDB)) {
                                showNotification(paste0("OrgDB package '", pkg, "' could not be loaded."), type = "error")
                                return(NULL)
                            }
                            
                            clusterProfiler::gseGO(
                                geneList = ranked,
                                OrgDb = orgDB,
                                ont = "BP",
                                keyType = "ENTREZID",
                                pvalueCutoff = input$pval_cutoff,
                                pAdjustMethod = "BH",
                                verbose = FALSE
                            )
                        }
                        gsea_result(res)
                        showNotification("GSEA completed!", type = "message")
                    },
                    error = function(e) {
                        showNotification(paste("Error:", conditionMessage(e)), type = "error")
                    }
                )
            })
        })

        # Pathway selector for ES plot
        output$gsea_pathway_ui <- renderUI({
            req(gsea_result())
            pathways <- gsea_result()@result$Description
            selectInput(ns("selected_pathway"), "Select pathway for ES plot:",
                choices = pathways, selected = pathways[1]
            )
        })

        # ── Outputs ──────────────────────────────────────────────────────────────
        output$gsea_table <- DT::renderDataTable({
            req(gsea_result())
            df <- as.data.frame(gsea_result())
            df <- df[, c("Description", "NES", "pvalue", "p.adjust", "setSize", "core_enrichment")]
            DT::datatable(round_numeric(df), options = list(scrollX = TRUE, pageLength = 10))
        })

        output$gsea_dotplot <- renderPlot({
            req(gsea_result())
            enrichplot::dotplot(gsea_result(), showCategory = 20) +
                ggplot2::ggtitle("GSEA Dotplot")
        })

        output$gsea_es_plot <- renderPlot({
            req(gsea_result(), input$selected_pathway)
            idx <- which(gsea_result()@result$Description == input$selected_pathway)
            if (length(idx) == 0) {
                return()
            }
            enrichplot::gseaplot2(gsea_result(),
                geneSetID = idx,
                title = input$selected_pathway
            )
        })

        output$gsea_leading <- DT::renderDataTable({
            req(gsea_result())
            df <- as.data.frame(gsea_result())
            leading <- data.frame(
                Pathway = df$Description,
                NES = round(df$NES, 4),
                pvalue = round(df$pvalue, 6),
                p.adjust = round(df$p.adjust, 6),
                LeadingEdge = df$core_enrichment,
                stringsAsFactors = FALSE
            )
            DT::datatable(leading, options = list(scrollX = TRUE, pageLength = 10))
        })

        # ── Downloads ─────────────────────────────────────────────────────────────
        output$dl_gsea_csv <- downloadHandler(
            filename = function() paste0("GSEA_results_", Sys.Date(), ".csv"),
            content  = function(file) write.csv(as.data.frame(gsea_result()), file, row.names = FALSE)
        )
        output$dl_gsea_dot <- downloadHandler(
            filename = function() paste0("GSEA_dotplot_", Sys.Date(), ".png"),
            content = function(file) {
                png(file, width = 1200, height = 800, res = 120)
                print(enrichplot::dotplot(gsea_result(), showCategory = 20))
                dev.off()
            }
        )
        output$dl_gsea_es <- downloadHandler(
            filename = function() paste0("GSEA_ES_plot_", Sys.Date(), ".png"),
            content = function(file) {
                req(input$selected_pathway)
                idx <- which(gsea_result()@result$Description == input$selected_pathway)
                png(file, width = 1200, height = 600, res = 120)
                print(enrichplot::gseaplot2(gsea_result(), geneSetID = idx, title = input$selected_pathway))
                dev.off()
            }
        )
    })
}

# helper: round numeric cols
round_numeric <- function(df) {
    num_cols <- sapply(df, is.numeric)
    df[, num_cols] <- lapply(df[, num_cols, drop = FALSE], function(x) round(x, 4))
    df
}

