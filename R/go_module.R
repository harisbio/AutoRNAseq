# ============================================================
# GO Enrichment Module
# AutoRNAseq v2
# ============================================================

goEnrichmentUI <- function(id) {
    ns <- NS(id)
    tagList(
        box(
            width = NULL,
            title = span(icon("circle-nodes"), " Gene Ontology (GO) Enrichment Analysis"),
            status = "primary",
            solidHeader = TRUE,
            style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
            sidebarLayout(
                sidebarPanel(
                    style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
                    h4("⚙️ Parameters", style = "color: #0092AC;"),
                    selectInput(ns("gene_id_type"), "Gene ID Type:",
                        choices = c("SYMBOL", "ENSEMBL", "ENTREZID"), selected = "SYMBOL"
                    ),
                    selectInput(ns("organism"), "Organism:",
                        choices = downstream_kegg_choices(),
                        selected = "hsa"
                    ),
                    selectInput(ns("ontology"), "GO Ontology:",
                        choices = c(
                            "Biological Process (BP)" = "BP",
                            "Molecular Function (MF)" = "MF",
                            "Cellular Component (CC)" = "CC"
                        ),
                        selected = "BP"
                    ),
                    numericInput(ns("padj_cutoff"), "Adj. p-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
                    numericInput(ns("lfc_cutoff"), "Min |log2FC|:", value = 1.0, min = 0, max = 10, step = 0.1),
                    hr(),
                    h4("📂 DEG Results Input", style = "color: #0092AC;"),
                    helpText("Upload a CSV file from DESeq2 results (needs columns: log2FoldChange, padj, or use gene names as row names)."),
                    fileInput(ns("deg_file"), "Upload DEG CSV (.csv)"),
                    uiOutput(ns("go_mapping_notice")),
                    uiOutput(ns("go_mapping_upload_ui")),
                    hr(),
                    div(
                        style = "display:flex; gap:10px;",
                        actionButton(ns("btn_run_go"), "Run GO Enrichment",
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
                            value = "go-table",
                            box(
                                width = NULL, status = "primary", solidHeader = TRUE,
                                collapsible = TRUE,
                                style = "background-color:#fff; border-radius:10px; padding:15px;",
                                DT::dataTableOutput(ns("go_table"))
                            ),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_go_csv"), "Download CSV",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("chart-bar"), " Dotplot"),
                            value = "go-dot",
                            plotOutput(ns("go_dotplot"), height = "500px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_dotplot"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("chart-column"), " Barplot"),
                            value = "go-bar",
                            plotOutput(ns("go_barplot"), height = "500px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_barplot"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("diagram-project"), " Gene-Concept Network"),
                            value = "go-net",
                            plotOutput(ns("go_cnet"), height = "600px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_cnet"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("network-wired"), " Enrichment Map"),
                            value = "go-emap",
                            plotOutput(ns("go_emap"), height = "600px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_emap"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        )
                    ),
                    br(),
                    div(
                        style = "display:flex; justify-content:left;",
                        actionButton(ns("help_go"), "Help",
                            icon = icon("info-circle"),
                            style = "color:#fff; background-color:#F39C12; border-color:#E67E22; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                )
            )
        )
    )
}


goEnrichmentServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        selected_orgdb <- reactive({
            downstream_orgdb_from_kegg(input$organism)
        })

        output$go_mapping_notice <- renderUI({
            orgdb <- selected_orgdb()
            if (is.na(orgdb) || !nzchar(orgdb)) {
                tags$div(
                    style = "margin-top:10px; color:#F4D03F; font-weight:600;",
                    "This species does not have a built-in OrgDb package. Upload a custom GO mapping file to run GO enrichment."
                )
            } else {
                helpText(paste0("Detected OrgDb package: ", orgdb))
            }
        })

        output$go_mapping_upload_ui <- renderUI({
            orgdb <- selected_orgdb()
            if (is.na(orgdb) || !nzchar(orgdb)) {
                tagList(
                    fileInput(ns("go_map_file"), "Upload GO Mapping (.csv/.tsv/.txt)"),
                    helpText("Required columns: gene_id and go_id (GO:0000000 format). Optional column: term_name."),
                    downloadButton(
                        ns("dl_go_mapping_template"),
                        "Download GO Mapping Template",
                        style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                    )
                )
            } else {
                NULL
            }
        })

        # Help modal
        observeEvent(input$help_go, {
            show_step_help(
              "GO Enrichment Help",
              c(
                "Upload a DEG CSV with gene names and the columns log2FoldChange and padj.",
                "Gene IDs may be SYMBOL, ENSEMBL, or ENTREZID depending on the dataset and organism.",
                "If the organism has no built-in OrgDb, upload a GO mapping file with gene_id and go_id columns.",
                "Use a small, filtered DEG file first to confirm the format before running a larger list.",
                "The GO result tables and plots are meant for interpretation after differential expression, not raw read input."
              ),
              "Typical DEG files are small to medium in size, but the pathway annotations can expand quickly if the gene list is very large."
            )
        })

        # ── reactive: read and filter DEG file ──────────────────────────────────
        deg_data <- reactive({
            req(input$deg_file)
            df <- read.csv(input$deg_file$datapath, header = TRUE, check.names = FALSE)
            # First column treated as gene ID if not already row names
            if (!"gene" %in% tolower(colnames(df))) {
                colnames(df)[1] <- "gene"
            }
            df
        })

        gene_list_filtered <- reactive({
            req(deg_data())
            df <- deg_data()
            lfc <- input$lfc_cutoff
            pco <- input$padj_cutoff

            # Filter by thresholds
            sig <- df[!is.na(df$padj) & !is.na(df$log2FoldChange) &
                df$padj < pco & abs(df$log2FoldChange) >= lfc, ]

            genes <- as.character(sig[[1]]) # first column = gene ids
            genes[!is.na(genes) & genes != ""]
        })

        entrez_ids <- reactive({
            req(gene_list_filtered(), input$organism, input$gene_id_type)
            genes <- gene_list_filtered()
            orgdb_pkg <- selected_orgdb()

            if (input$gene_id_type == "ENTREZID" || is.na(orgdb_pkg) || !nzchar(orgdb_pkg)) {
                return(genes)
            }

            pkg <- orgdb_pkg
            if (!requireNamespace(pkg, quietly = TRUE)) {
                showNotification(paste0("OrgDB package '", pkg, "' is not installed in this container."), type = "error")
                return(NULL)
            }
            orgDB <- tryCatch(getExportedValue(pkg, pkg), error = function(e) NULL)
            if (is.null(orgDB)) {
                showNotification(paste0("OrgDB package '", pkg, "' could not be loaded. Please try again."), type = "error")
                return(NULL)

            }

            mapped <- AnnotationDbi::mapIds(orgDB,
                keys = genes,
                column = "ENTREZID",
                keytype = input$gene_id_type,
                multiVals = "first"
            )
            mapped <- mapped[!is.na(mapped)]
            as.character(mapped)
        })

        go_mapping <- reactive({
            req(input$go_map_file)
            file_name <- tolower(input$go_map_file$name)
            ext_tsv <- grepl("\\.(tsv|txt)$", file_name)
            df <- tryCatch(
                {
                    if (ext_tsv) {
                        read.delim(input$go_map_file$datapath, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
                    } else {
                        read.csv(input$go_map_file$datapath, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
                    }
                },
                error = function(e) NULL
            )
            if (is.null(df) || ncol(df) < 2) {
                showNotification("GO mapping file is invalid. Provide at least gene_id and go_id columns.", type = "error")
                return(NULL)
            }

            cn <- tolower(trimws(colnames(df)))
            gene_col <- which(cn %in% c("gene", "gene_id", "geneid", "symbol", "id"))[1]
            go_col <- which(cn %in% c("go", "go_id", "goid", "go_term"))[1]
            term_col <- which(cn %in% c("term_name", "term", "description", "name"))[1]
            if (is.na(gene_col)) gene_col <- 1
            if (is.na(go_col)) go_col <- 2

            raw_gene <- as.character(df[[gene_col]])
            raw_go <- as.character(df[[go_col]])
            split_go <- strsplit(raw_go, "[,;|\\s]+")
            pieces <- lapply(seq_along(raw_gene), function(i) {
                g <- trimws(raw_gene[i])
                gos <- unique(trimws(split_go[[i]]))
                gos <- gos[nzchar(gos) & grepl("^GO:\\d{7}$", gos)]
                if (!nzchar(g) || length(gos) == 0) {
                    return(NULL)
                }
                data.frame(term = gos, gene = g, stringsAsFactors = FALSE)
            })
            pieces <- pieces[!vapply(pieces, is.null, logical(1))]
            if (length(pieces) == 0) {
                showNotification("No valid GO IDs found (expected format: GO:0008150).", type = "error")
                return(NULL)
            }
            term2gene <- unique(do.call(rbind, pieces))

            term2name <- NULL
            if (!is.na(term_col)) {
                tmp <- data.frame(
                    term = as.character(df[[go_col]]),
                    name = as.character(df[[term_col]]),
                    stringsAsFactors = FALSE
                )
                tmp <- tmp[grepl("^GO:\\d{7}$", tmp$term) & nzchar(trimws(tmp$name)), , drop = FALSE]
                if (nrow(tmp) > 0) {
                    term2name <- unique(tmp)
                }
            }

            list(term2gene = term2gene, term2name = term2name)
        })

        # ── reactive: run GO ─────────────────────────────────────────────────────
        go_result <- reactiveVal(NULL)

        observeEvent(input$btn_run_go, {
            ids <- entrez_ids()
            req(ids)
            if (length(ids) == 0) {
                showNotification("No significant genes found with current thresholds.", type = "warning")
                return()
            }
            withProgress(message = "Running GO Enrichment...", value = 0.5, {
                tryCatch(
                    {
                        orgdb_pkg <- selected_orgdb()
                        res <- NULL
                        if (!is.na(orgdb_pkg) && nzchar(orgdb_pkg)) {
                            orgDB <- getExportedValue(orgdb_pkg, orgdb_pkg)
                            res <- clusterProfiler::enrichGO(
                                gene = ids,
                                OrgDb = orgDB,
                                keyType = "ENTREZID",
                                ont = input$ontology,
                                pAdjustMethod = "BH",
                                pvalueCutoff = input$padj_cutoff,
                                readable = TRUE
                            )
                        } else {
                            mapping <- go_mapping()
                            req(mapping)
                            res <- clusterProfiler::enricher(
                                gene = unique(ids),
                                TERM2GENE = mapping$term2gene,
                                TERM2NAME = mapping$term2name,
                                pvalueCutoff = input$padj_cutoff,
                                pAdjustMethod = "BH"
                            )
                        }
                        if (is.null(res) || nrow(as.data.frame(res)) == 0) {
                            showNotification("No enriched GO terms found for selected thresholds.", type = "warning")
                            go_result(NULL)
                            return()
                        }
                        go_result(res)
                        showNotification("GO Enrichment completed!", type = "message")
                    },
                    error = function(e) {
                        showNotification(paste("Error:", conditionMessage(e)), type = "error")
                    }
                )
            })
        })

        # ── outputs ───────────────────────────────────────────────────────────────
        output$go_table <- DT::renderDataTable({
            req(go_result())
            DT::datatable(as.data.frame(go_result()), options = list(scrollX = TRUE, pageLength = 10))
        })

        output$go_dotplot <- renderPlot({
            req(go_result())
            enrichplot::dotplot(go_result(), showCategory = 20) +
                ggplot2::ggtitle("GO Enrichment Dotplot")
        })

        output$go_barplot <- renderPlot({
            req(go_result())
            clusterProfiler::barplot(go_result(), showCategory = 20) +
                ggplot2::ggtitle("GO Enrichment Barplot")
        })

        output$go_cnet <- renderPlot({
            req(go_result())
            enrichplot::cnetplot(go_result(), showCategory = 5) +
                ggplot2::ggtitle("Gene-Concept Network")
        })

        output$go_emap <- renderPlot({
            req(go_result())
            res2 <- enrichplot::pairwise_termsim(go_result())
            enrichplot::emapplot(res2) +
                ggplot2::ggtitle("Enrichment Map")
        })

        # ── downloads ─────────────────────────────────────────────────────────────
        output$dl_go_mapping_template <- downloadHandler(
            filename = function() {
                paste0("GO_mapping_template_", Sys.Date(), ".csv")
            },
            content = function(file) {
                template <- data.frame(
                    gene_id = c("GeneA", "GeneB", "GeneC"),
                    go_id = c(
                        "GO:0008150;GO:0009987",
                        "GO:0003674",
                        "GO:0005575|GO:0005634"
                    ),
                    term_name = c(
                        "biological_process;cellular process",
                        "molecular_function",
                        "cellular_component;nucleus"
                    ),
                    stringsAsFactors = FALSE
                )
                write.csv(template, file, row.names = FALSE)
            }
        )

        output$dl_go_csv <- downloadHandler(
            filename = function() paste0("GO_enrichment_", Sys.Date(), ".csv"),
            content  = function(file) write.csv(as.data.frame(go_result()), file, row.names = FALSE)
        )

        make_plot_dl <- function(plot_expr, fname) {
            downloadHandler(
                filename = function() paste0(fname, "_", Sys.Date(), ".png"),
                content = function(file) {
                    png(file, width = 1200, height = 800, res = 120)
                    print(plot_expr())
                    dev.off()
                }
            )
        }

        output$dl_dotplot <- make_plot_dl(
            function() enrichplot::dotplot(go_result(), showCategory = 20), "GO_dotplot"
        )
        output$dl_barplot <- make_plot_dl(
            function() clusterProfiler::barplot(go_result(), showCategory = 20), "GO_barplot"
        )
        output$dl_cnet <- make_plot_dl(
            function() enrichplot::cnetplot(go_result(), showCategory = 5), "GO_cnetplot"
        )
        output$dl_emap <- make_plot_dl(function() {
            res2 <- enrichplot::pairwise_termsim(go_result())
            enrichplot::emapplot(res2)
        }, "GO_emapplot")
    })
}

