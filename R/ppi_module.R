# ============================================================
# PPI Network Analysis Module
# AutoRNAseq v2
# ============================================================

ppiNetworkUI <- function(id) {
    ns <- NS(id)
    tagList(
        box(
            width = NULL,
            title = span(icon("circle-nodes"), " Protein–Protein Interaction Network"),
            status = "primary",
            solidHeader = TRUE,
            style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
            sidebarLayout(
                sidebarPanel(
                    style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
                    h4("⚙️ STRINGdb Parameters", style = "color: #0092AC;"),
                    selectInput(ns("ppi_species"), "Species (NCBI Taxonomy ID):",
                        choices = downstream_ppi_choices(),
                        selected = "9606"
                    ),
                    sliderInput(ns("string_score"), "STRING Confidence Score (0–1000):",
                        min = 0, max = 1000, value = 400, step = 50
                    ),
                    numericInput(ns("top_n_genes"), "Max hub genes to highlight:", value = 10, min = 1, max = 50, step = 1),
                    hr(),
                    h4("📂 DEG Results Input", style = "color: #0092AC;"),
                    fileInput(ns("deg_file"), "Upload DEG CSV (.csv)"),
                    numericInput(ns("padj_cutoff"), "Adj. p-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
                    numericInput(ns("lfc_cutoff"), "Min |log2FC|:", value = 1.0, min = 0, max = 10, step = 0.1),
                    hr(),
                    div(
                        style = "display:flex; gap:10px;",
                        actionButton(ns("btn_run_ppi"), "Build PPI Network",
                            icon = icon("play"),
                            style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                ),
                mainPanel(
                    tabBox(
                        width = 12,
                        tabPanel(
                            title = span(icon("circle-nodes"), " Interactive Network"),
                            value = "ppi-network",
                            visNetwork::visNetworkOutput(ns("ppi_network"), height = "600px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_ppi_html"), "Download Network (HTML)",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("trophy"), " Hub Genes"),
                            value = "ppi-hubs",
                            box(
                                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                style = "background-color:#fff; border-radius:10px; padding:15px;",
                                DT::dataTableOutput(ns("hub_table"))
                            ),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_hub_csv"), "Download Hub Genes CSV",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        ),
                        tabPanel(
                            title = span(icon("chart-bar"), " Centrality Plot"),
                            value = "ppi-centrality",
                            plotOutput(ns("centrality_plot"), height = "500px"),
                            br(),
                            div(
                                style = "display:flex; gap:10px;",
                                downloadButton(ns("dl_centrality"), "Download Plot",
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                                )
                            )
                        )
                    ),
                    br(),
                    div(
                        style = "display:flex; justify-content:left;",
                        actionButton(ns("help_ppi"), "Help",
                            icon = icon("info-circle"),
                            style = "color:#fff; background-color:#F39C12; border-color:#E67E22; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                )
            )
        )
    )
}


ppiNetworkServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        observeEvent(input$help_ppi, {
            show_step_help(
              "PPI Network Help",
              c(
                "Upload a DEG CSV with gene names and the columns log2FoldChange and padj.",
                "STRINGdb uses the species taxonomy ID, so select the correct organism before building the network.",
                "A confidence score around 400 is a sensible starting point for most exploratory analyses.",
                "Use a smaller filtered DEG list first, because large networks can get slow or cluttered.",
                "Hub gene rankings and centrality plots are best interpreted together with DEG and pathway results."
              ),
              "If your dataset is large, the network itself can become visually dense even when the input file is small."
            )
        })

        # ── Load genes ───────────────────────────────────────────────────────────
        sig_genes <- reactive({
            req(input$deg_file)
            df <- read.csv(input$deg_file$datapath, header = TRUE, check.names = FALSE)
            colnames(df)[1] <- "gene"
            sig <- df[!is.na(df$padj) & !is.na(df$log2FoldChange) &
                df$padj < input$padj_cutoff &
                abs(df$log2FoldChange) >= input$lfc_cutoff, ]
            as.character(sig$gene)
        })

        # ── STRINGdb calls ───────────────────────────────────────────────────────
        ppi_data <- reactiveVal(NULL)

        observeEvent(input$btn_run_ppi, {
            genes <- sig_genes()
            req(genes)
            if (length(genes) == 0) {
                showNotification("No significant genes to build network with.", type = "warning")
                return()
            }

            withProgress(message = "Querying STRING database...", value = 0.3, {
                tryCatch(
                    {
                        string_db <- STRINGdb::STRINGdb$new(
                            version = "11.5",
                            species = as.integer(input$ppi_species),
                            score_threshold = as.integer(input$string_score),
                            input_directory = tempdir()
                        )

                        gene_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
                        mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
                        setProgress(0.6, message = "Building network...")

                        interactions <- string_db$get_interactions(mapped$STRING_id)
                        ppi_data(list(mapped = mapped, interactions = interactions, string_db = string_db))
                        showNotification("PPI Network built successfully!", type = "message")
                    },
                    error = function(e) {
                        showNotification(paste("STRINGdb Error:", conditionMessage(e)), type = "error")
                    }
                )
            })
        })

        # ── compute graph metrics ────────────────────────────────────────────────
        graph_metrics <- reactive({
            req(ppi_data())
            dat <- ppi_data()
            interactions <- dat$interactions
            mapped <- dat$mapped

            if (is.null(interactions) || nrow(interactions) == 0) {
                showNotification("No interactions found at this confidence threshold.", type = "warning")
                return(NULL)
            }

            g <- igraph::graph_from_data_frame(
                d         = interactions[, c("from", "to")],
                directed  = FALSE,
                vertices  = data.frame(name = mapped$STRING_id)
            )

            mapped_genes <- mapped[match(igraph::V(g)$name, mapped$STRING_id), ]
            igraph::V(g)$label <- mapped_genes$gene
            igraph::V(g)$degree <- igraph::degree(g)
            igraph::V(g)$between <- igraph::betweenness(g, normalized = TRUE)
            igraph::V(g)$close <- igraph::closeness(g, normalized = TRUE)
            g
        })

        # ── Interactive network via visNetwork ────────────────────────────────────
        output$ppi_network <- visNetwork::renderVisNetwork({
            req(graph_metrics())
            g <- graph_metrics()

            nodes <- data.frame(
                id = igraph::V(g)$name,
                label = igraph::V(g)$label,
                value = igraph::V(g)$degree,
                title = paste0(
                    "<b>", igraph::V(g)$label, "</b><br>",
                    "Degree: ", igraph::V(g)$degree, "<br>",
                    "Betweenness: ", round(igraph::V(g)$between, 4)
                ),
                color = colorRampPalette(c("#AED6F1", "#1A5276"))(100)[
                    cut(igraph::V(g)$degree, breaks = 100, labels = FALSE, include.lowest = TRUE)
                ],
                stringsAsFactors = FALSE
            )

            edges <- data.frame(
                from = igraph::ends(g, igraph::E(g))[, 1],
                to = igraph::ends(g, igraph::E(g))[, 2],
                stringsAsFactors = FALSE
            )

            visNetwork::visNetwork(nodes, edges, height = "580px", width = "100%") %>%
                visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
                visNetwork::visLayout(randomSeed = 42) %>%
                visNetwork::visPhysics(stabilization = TRUE) %>%
                visNetwork::visEdges(smooth = FALSE, color = list(color = "#cccccc")) %>%
                visNetwork::visNodes(scaling = list(min = 10, max = 40))
        })

        # ── Hub gene table ────────────────────────────────────────────────────────
        hub_df <- reactive({
            req(graph_metrics())
            g <- graph_metrics()
            data.frame(
                Gene = igraph::V(g)$label,
                Degree = igraph::V(g)$degree,
                Betweenness = round(igraph::V(g)$between, 4),
                Closeness = round(igraph::V(g)$close, 4),
                stringsAsFactors = FALSE
            ) %>% dplyr::arrange(dplyr::desc(Degree))
        })

        output$hub_table <- DT::renderDataTable({
            req(hub_df())
            DT::datatable(hub_df(), options = list(scrollX = TRUE, pageLength = 15))
        })

        # ── Centrality bar plot ────────────────────────────────────────────────────
        output$centrality_plot <- renderPlot({
            req(hub_df())
            df <- head(hub_df(), input$top_n_genes)
            ggplot2::ggplot(df, ggplot2::aes(x = reorder(Gene, Degree), y = Degree, fill = Betweenness)) +
                ggplot2::geom_col() +
                ggplot2::coord_flip() +
                ggplot2::scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
                ggplot2::labs(
                    title = paste("Top", input$top_n_genes, "Hub Genes by Degree Centrality"),
                    x = "Gene", y = "Degree", fill = "Betweenness"
                ) +
                ggplot2::theme_minimal(base_size = 13)
        })

        # ── Downloads ─────────────────────────────────────────────────────────────
        output$dl_hub_csv <- downloadHandler(
            filename = function() paste0("PPI_hub_genes_", Sys.Date(), ".csv"),
            content  = function(file) write.csv(hub_df(), file, row.names = FALSE)
        )
        output$dl_centrality <- downloadHandler(
            filename = function() paste0("PPI_centrality_", Sys.Date(), ".png"),
            content = function(file) {
                df <- head(hub_df(), input$top_n_genes)
                p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(Gene, Degree), y = Degree, fill = Betweenness)) +
                    ggplot2::geom_col() +
                    ggplot2::coord_flip() +
                    ggplot2::scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
                    ggplot2::labs(title = "Hub Genes — Degree Centrality", x = "Gene", y = "Degree") +
                    ggplot2::theme_minimal(base_size = 13)
                png(file, width = 1200, height = 800, res = 120)
                print(p)
                dev.off()
            }
        )
        output$dl_ppi_html <- downloadHandler(
            filename = function() paste0("PPI_network_", Sys.Date(), ".html"),
            content = function(file) {
                req(graph_metrics())
                g <- graph_metrics()
                nodes <- data.frame(
                    id = igraph::V(g)$name, label = igraph::V(g)$label,
                    value = igraph::V(g)$degree, stringsAsFactors = FALSE
                )
                edges <- data.frame(
                    from = igraph::ends(g, igraph::E(g))[, 1],
                    to = igraph::ends(g, igraph::E(g))[, 2], stringsAsFactors = FALSE
                )
                net <- visNetwork::visNetwork(nodes, edges) %>%
                    visNetwork::visOptions(highlightNearest = TRUE)
                visNetwork::visSave(net, file = file)
            }
        )
    })
}

