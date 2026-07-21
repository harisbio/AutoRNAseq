dif_gene_ex_analysisUI <- function(id) {
  ns <- NS(id)
  tagList(


    # 🔹 Main UI Box
    box(
      width = NULL,
      title = span(icon("dna"), " Differential Gene Expression Analysis"),
      status = "primary",
      solidHeader = TRUE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",

      tabBox(
        width = 12,

        # 📂 Step 1: Upload Experimental Design
        tabPanel(
          title = span(icon("upload"), " Upload Experimental Design"),
          value = "tab-metadata",
          sidebarLayout(
            sidebarPanel(
              style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
              h4("📋 Upload the metadata of your study"),
              fileInput(ns("phenotype_data"), "Import metadata (.csv)"),
              uiOutput(ns("desc_columns"))
            ),
            mainPanel(
              box(
                title = span(icon("table"), " Experimental Design Preview"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("desc_phenotype_data"))
              )
            )
          )
        ),

        # 📊 Step 2: Upload Count Matrix
        tabPanel(
          title = span(icon("upload"), " Upload Count Matrix"),
          value = "tab-countData",
          sidebarLayout(
            sidebarPanel(
              style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
              h4("📊 Upload your count matrix"),
              fileInput(ns("gene_expression_data"), "Import count data (.csv)")
            ),
            mainPanel(
              tabBox(
                width = 12,
                tabPanel(
                  title = span(icon("table"), " Data Table"),
                  value = "tab_countMatrix",
                  box(
                    title = span(icon("table"), " Count Matrix Preview"),
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                    DT::dataTableOutput(ns("count_data"))
                  )
                ),
                tabPanel(
                  title = span(icon("table"), " Count Data with Experimental Design"),
                  value = "tab_exp_count_data",
                  box(
                    title = span(icon("table"), " Experimental Design Data"),
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                    DT::dataTableOutput(ns("desc_expression_data"))
                  )
                ),

                tabPanel(
                  title = span(icon("image"), " Visualize Gene Expression"),
                  value = "tab_visual_count",
                  uiOutput(ns("visual_gene_expression"))
                )
              )
            )
          )
        ),

        # 🔬 Step 3: Run Differential Expression Analysis
        tabPanel(
          title = span(icon("play"), " Run Differential Expression Analysis"),
          value = "tab-DE",
          sidebarLayout(
            sidebarPanel(
              style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
              h4("🧪 Pre-filtering: Removing genes with low counts"),
              textInput(ns("filtering_factor"), "Enter the minimum count threshold", value = "10"),
              textInput(ns("input_alpha"), "Enter statistical significance parameter (alpha)", value = "0.05"),
              hr(),
              actionButton(ns("btn_run_DFEXR"), "Run Analysis", icon = icon("play"),
                           style = "color: #ffffff; background-color: #0092AC;
                           border-color: #007B9E; padding: 6px 12px; font-size: 14px;
                           border-radius: 5px; width: auto;")
            ),
            mainPanel(
              box(
                title = span(icon("table"), " Differential Expression Results"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("results"))
              ),
              br(),
              div(
                style = "display: flex; justify-content: left;",
                downloadButton(ns("btn_dwn_Degs"),
                               "Download Results",
                               style = "color: #ffffff; background-color: #0092AC;
               border-color: #007B9E; padding: 6px 12px; font-size: 14px;
               border-radius: 5px; width: auto;")

              )
              ,
              br(),
              div(
                style = "display: flex; justify-content: left;",
                help_modal_button(ns("help_de"), "Help")
              )
              # downloadButton(ns("btn_dwn_Degs"), "Download Results",
              #                style = "padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;")
            )
          )
        ),

        # 📈 Step 4: Summary Plots
        tabPanel(
          title = span(icon("image"), " Summary Plots"),
          value = "tab-plots",
          sidebarLayout(
            sidebarPanel(
              style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
              h4("📊 Select Plot Type", style = "color: #0092AC;"),
              selectInput(ns("select_plot_type"), "Choose a visualization:", choices = c(
                "None",
                "MA plot",
                "Shrunken LFC MA Plot",
                "Dispersion plot",
                "PCA plot",
                "Volcano plot",
                "Heatmap of Pairwise Sample Distances",
                "Heatmap of Z Scores for the Top Genes",
                "Summary Plot for DEGs"
              )),
              uiOutput(ns("plot_controls")),
              textOutput(ns("plot_desc"))
            ),
            mainPanel(

                    plotOutput(ns("plot_summary"), height = "430px"),
                  # 📥 Download Filtering report Button (Small)
                  div(
                    style = "display: flex; justify-content: left;",
                    downloadButton(ns("download_plot"), "Download Plot",
                                   style = "color: #ffffff; background-color: #0092AC;
                             border-color: #007B9E; padding: 6px 12px; font-size: 14px;
                             border-radius: 5px; width: auto;"))

            )
          )
        )
      )
    )
  )
}



dif_gene_ex_analysisServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Force max upload size dynamically at module load time
    options(shiny.maxRequestSize = 8192 * 1024^2)

    #-------------------------------------------------------------------------
    # Aesthetic Overhaul: Premium Theme
    #-------------------------------------------------------------------------
    theme_seq_premium <- function() {
      theme_minimal() +
      theme(
        text = element_text(family = "sans", color = "#2C3E50"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        plot.subtitle = element_text(size = 14, color = "#7F8C8D", hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid.major = element_line(color = "#ECF0F1", size = 0.5),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
        plot.margin = margin(30, 30, 30, 30)
      )
    }

    library(dplyr)
    library(tidyr)
    library(DESeq2)
    library(DT)
    library(shiny)
    library(ggplot2)
    library(EnhancedVolcano)
    library(pheatmap)
    library(RColorBrewer)
    library(apeglm)

    #-------------------------------------------------------------------------
    # Function to read the metadata
    #-------------------------------------------------------------------------
    observeEvent(input$help_de, {
      show_step_help(
        "Differential Expression Help",
        c(
          "Upload a metadata CSV with one row per sample and columns for sample names and conditions.",
          "Upload a count matrix CSV with genes as rows and sample names as columns.",
          "Keep the sample names identical between metadata and count matrix.",
          "Choose a minimum count threshold and alpha value before running DESeq2.",
          "A small pilot dataset is the best way to verify the file format before using a large study."
        ),
        "For smooth results, use raw count data from the quantification step rather than normalized values."
      )
    })

    get_metadata <- reactive({
      req(input$phenotype_data$datapath)  # Use req() to check if the file path is provided
      pheno_data <- read.csv(input$phenotype_data$datapath, header = TRUE)
      rownames(pheno_data) <- NULL
      print(rownames(pheno_data))
      return(pheno_data)
    })
    #-------------------------------------------------------------------------
    # Render metadata table
    #-------------------------------------------------------------------------
    output$desc_phenotype_data <- DT::renderDataTable({
      if (is.null(input$phenotype_data)) {
        return(data.frame(Message = "Upload the info on the experimental design as a .csv file"))
      }
      pheno_data <- get_metadata()
      rownames(pheno_data) <- NULL
      DT::datatable(pheno_data, options = list(scrollX = TRUE))
    })

    #-------------------------------------------------------------------------
    # Render UI for selecting metadata columns
    #-------------------------------------------------------------------------
    output$desc_columns <- renderUI({
      req(input$phenotype_data$datapath)  # Check if phenotype data is selected
      pheno_data <-  get_metadata()
      tagList(
        h4("Select the targeted columns for the analysis"),
        selectInput(
          ns("conditions_column"),
          "Choose the conditions column",
          choices = c("None", colnames(pheno_data))
        ),
        selectInput(
          ns("samples_column"),
          "Choose the samples names column",
          choices = c("None", colnames(pheno_data))
        )
      )

    })

    #-------------------------------------------------------------------------
    # Modified metadata based on selected columns
    #-------------------------------------------------------------------------
    modified_metadata <- reactive({
      req(input$conditions_column,
          input$samples_column,
          get_metadata())
      if (input$conditions_column == "None" || input$samples_column == "None") {
        return(NULL)
      }
      pheno_data <-  get_metadata()
      rownames(pheno_data) <- NULL
      if (!(input$conditions_column %in% colnames(pheno_data)) || !(input$samples_column %in% colnames(pheno_data))) {
        return(NULL)
      }
      metadata <- pheno_data %>%
        select(input$conditions_column, input$samples_column) %>%
        dplyr::rename(conditions = colnames(.[1]), samples = colnames(.[2]))
      print(metadata)
      return(metadata)
    })




    #-------------------------------------------------------------------------
    # Function to read count matrix data
    #-------------------------------------------------------------------------
    get_gene_data <- reactive({
      req(input$gene_expression_data$datapath)  # Check if gene expression data is selected
      gene_data <- read.csv(input$gene_expression_data$datapath, header = TRUE)
      colnames(gene_data)[1] <- "geneID"
      colnames(gene_data) <- gsub("\\.bam$", "", colnames(gene_data), ignore.case = TRUE)
      return(gene_data)
    })
    #-------------------------------------------------------------------------
    # Render count matrix table
    #-------------------------------------------------------------------------
    output$count_data <- DT::renderDataTable({
      if (is.null(input$gene_expression_data)) {
        return(data.frame(Message = "Upload your count matrix as a .csv file"))
      }

      gene_data <- get_gene_data()


      DT::datatable(gene_data, options = list(scrollX = TRUE))
    })


    #-------------------------------------------------------------------------
    # Join gene expression data with modified metadata
    #-------------------------------------------------------------------------
    join_data <- reactive({
      req(get_gene_data())  # Ensure gene data exists
      req(modified_metadata())  # Ensure metadata exists

      gene_data <- get_gene_data()
      meta <- modified_metadata()
      req(!is.null(meta))

      # Check if gene_data has at least one row
      req(nrow(gene_data) > 0)

      gene_data <- gene_data %>%
        gather(key = 'samples', value = 'count', -geneID) %>%
        dplyr::left_join(., meta, by = c("samples" = "samples"))

      return(gene_data)
    })



    #-------------------------------------------------------------------------
    # Render count matrix data table with the metadata
    #-------------------------------------------------------------------------
    output$desc_expression_data <- DT::renderDataTable({
      if (is.null(input$gene_expression_data) ||
          is.null(input$phenotype_data) ||
          input$conditions_column == 'None' ||
          input$samples_column == "None" || is.null(get_metadata())) {
        return(
          data.frame(Message = "Check your metadata file or your count matrix or select the targeted columns")
        )
      }
      DT::datatable(join_data(), options = list(scrollX = TRUE))
    })

    #-------------------------------------------------------------------------
    # Render UI for count matrix data visualization
    #-------------------------------------------------------------------------
    output$visual_gene_expression <- renderUI({
      list(
        shiny::fluidRow(
          shiny::column(6, textInput(ns("geneName"), "Enter the name of gene")),
          shiny::column(6, radioButtons(
            ns("choose_plot_type"),
            "Choose the type of visualization plot",
            choices = c("Bar plot")
          ))
        ),
           plotOutput(ns("visualization"), height = "430px"))


    })

    #-------------------------------------------------------------------------
    # Render selected gene expression data visualization
    #-------------------------------------------------------------------------
    output$visualization <- renderPlot({
      # Existing check logic remains same...

      # Check required inputs
      if (is.null(input$gene_expression_data) ||
          is.null(input$phenotype_data) ||
          is.null(input$choose_plot_type) ||
          is.null(modified_metadata())) {
        return(plot(
          1, type = "n", xlab = "", ylab = "", xlim = c(0, 5), ylim = c(0, 5),
          main = "Plot Will Appear Here Once Data is Selected"
        ))
      }
      
      plot_type <- input$choose_plot_type
      data <- join_data()
      if (is.null(data)) {
        return(plot(1, type = "n", xlab = "", ylab = "", main = "Please configure metadata columns correctly."))
      }
      
      if (plot_type == "Bar plot") {
        if (is.null(input$geneName) || trimws(input$geneName) == "") {
           return(plot(1, type="n", xlab="", ylab="", main="Please type a specific gene name in the text box."))
        }
        
        sub_data <- data %>% filter(geneID == input$geneName)
        if (nrow(sub_data) == 0) {
           return(plot(1, type="n", xlab="", ylab="", main=paste("Gene '", input$geneName, "' not found in expression data.", sep="")))
        }
        
        return(ggplot(sub_data, aes(x = samples, y = count, fill = conditions)) +
          geom_col(width = 0.6) +
          theme_seq_premium() +
          scale_fill_manual(values = c("firebrick", "steelblue", "forestgreen", "darkorchid")) +
          labs(title = paste("Expression Analysis: ", input$geneName),
               subtitle = "Normalized cross-sample counts",
               x = "Sample ID", y = "Normalized Counts"))
      }
    }, width = 700, height = 430) 


    #-------------------------------------------------------------------------
    # Run DESq2: Differential expression analysis
    #-------------------------------------------------------------------------
    run_DGEA <- reactive({
      # Get user input
      refColumn <- input$conditions_column
      sampleColumn <- input$samples_column
      data <- get_gene_data()
      metadata <- get_metadata()
      req(!is.null(data), !is.null(metadata))
      if (refColumn == "None" || sampleColumn == "None") {
        showNotification("Please select valid condition and sample columns in metadata.", type = "error")
        return(NULL)
      }
      if (!(refColumn %in% colnames(metadata)) || !(sampleColumn %in% colnames(metadata))) {
        showNotification("Selected metadata columns are invalid.", type = "error")
        return(NULL)
      }

      alpha_val <- suppressWarnings(as.numeric(input$input_alpha))
      filter_val <- suppressWarnings(as.numeric(input$filtering_factor))
      if (is.na(alpha_val) || alpha_val <= 0 || alpha_val >= 1) {
        showNotification("Alpha must be a number between 0 and 1.", type = "error")
        return(NULL)
      }
      if (is.na(filter_val) || filter_val < 0) {
        showNotification("Minimum count threshold must be a non-negative number.", type = "error")
        return(NULL)
      }

      print(data)
      print(rownames(data))

      # Ensure first column is set as row names and remove it
      rownames(data) <- data[[1]]
      data <- data[, -1, drop = FALSE]
      data <- data[, sapply(data, is.numeric), drop = FALSE]
      if (ncol(data) < 2) {
        showNotification("Count matrix must include at least two numeric sample columns.", type = "error")
        return(NULL)
      }

      metadata <- metadata[, c(sampleColumn, refColumn), drop = FALSE]
      colnames(metadata) <- c("sample", "condition")
      metadata$sample <- as.character(metadata$sample)
      metadata$condition <- as.character(metadata$condition)
      metadata <- metadata[!is.na(metadata$sample) & nzchar(metadata$sample), , drop = FALSE]

      common_samples <- intersect(colnames(data), metadata$sample)
      if (length(common_samples) < 2) {
        showNotification("Metadata sample names do not match count matrix columns.", type = "error")
        return(NULL)
      }
      data <- data[, common_samples, drop = FALSE]
      metadata <- metadata[match(common_samples, metadata$sample), , drop = FALSE]
      rownames(metadata) <- metadata$sample
      metadata$sample <- NULL

      # Convert count data to a matrix of integers
      data <- as.matrix(data)
      storage.mode(data) <- "integer"

      # Convert metadata reference column to factor
      if (!is.factor(metadata$condition)) {
        metadata$condition <- as.factor(metadata$condition)
      }
      if (nlevels(metadata$condition) < 2) {
        showNotification("At least two condition groups are required for DESeq2.", type = "error")
        return(NULL)
      }

      # Ensure metadata is a data frame
      metadata <- as.data.frame(metadata)

      # Construct DESeq2 dataset
      dds <- DESeqDataSetFromMatrix(
        countData = data,
        colData = metadata,
        design = ~ condition
      )

      # Filter genes with at least X reads
      keep <- rowSums(counts(dds)) >= filter_val
      dds <- dds[keep, ]
      if (nrow(dds) == 0) {
        showNotification("No genes remain after filtering. Reduce the filtering threshold.", type = "warning")
        return(NULL)
      }

      # Run DESeq2 analysis
      dds <- DESeq(dds)
      res <- results(dds, alpha = alpha_val)

      # Convert S4 results object to a data frame before filtering
      res_df <- as.data.frame(res)

      # Ensure complete cases before filtering
      res_df <- res_df[complete.cases(res_df), ]

      # Round numeric values for better readability
      numeric_cols <- sapply(res_df, is.numeric)
      res_df[, numeric_cols] <- round(res_df[, numeric_cols], 4)

      return(list(res_df, dds, res, metadata))
    })

    reactive_DGE <- reactiveVal(NULL)


    observeEvent(input$btn_run_DFEXR, {
      results <- tryCatch(
        run_DGEA(),
        error = function(e) {
          showNotification(paste("Differential analysis error:", conditionMessage(e)), type = "error")
          NULL
        }
      )

      if (is.null(results)) {
        showNotification("Differential gene expression analysis did not run. Please check that all inputs are correctly provided.", type = "error")
      } else {
        showNotification("Differential gene expression analysis completed successfully.", type = "message")
      }
      reactive_DGE(results)
    })

    # Libraries now loaded at the top of moduleServer

    cal_z_scores <- function(x) {
      (x - mean(x)) / sd(x)
    }



    #-------------------------------------------------------------------------
    #Render the table results
    #-------------------------------------------------------------------------
    output$results <- DT::renderDataTable({
      if (is.null(input$gene_expression_data)) {
        return(
          data.frame(Message = " Setup your data for the analysis (Experimental design/count matrix)")
        )
      }
      resultData <- reactive_DGE()
      req(resultData)
      DT::datatable(as.data.frame(resultData[[1]]) , options = list(scrollX = TRUE))
    })

    #-------------------------------------------------------------------------
    # UI controls for plots
    #-------------------------------------------------------------------------
    output$plot_controls <- renderUI({
      req(input$select_plot_type)
      plot_type <- input$select_plot_type
      
      if (plot_type == "Volcano plot") {
        tagList(
          numericInput(ns("volcano_lfc"), "Log2 Fold Change Cutoff:", value = 1.0, step = 0.1),
          numericInput(ns("volcano_pval"), "P-value Cutoff (padj):", value = 0.05, step = 0.01)
        )
      } else if (plot_type == "Heatmap of Pairwise Sample Distances") {
        numericInput(ns("heatmap_sample_top_n"), "Number of top variable genes to use:", value = 500, min = 10, step = 10)
      } else if (plot_type == "Heatmap of Z Scores for the Top Genes") {
        numericInput(ns("heatmap_top_n"), "Number of top genes:", value = 20, min = 5, step = 5)
      } else if (plot_type == "Summary Plot for DEGs") {
        tagList(
          numericInput(ns("summary_lfc"), "Log2 Fold Change Cutoff:", value = 1.0, step = 0.1),
          numericInput(ns("summary_pval"), "P-value Cutoff (padj):", value = 0.05, step = 0.01)
        )
      } else {
        NULL
      }
    })

    #-------------------------------------------------------------------------
    #Render summary plots
    #-------------------------------------------------------------------------
    output$plot_summary <- renderPlot({
      # Check if selected_file is not NULL or empty
      selectedPlot <- input$select_plot_type

      if (is.null(input$gene_expression_data) ||
          is.null(input$phenotype_data) || selectedPlot == "None") {
        return(plot(
          1,
          type = "n",
          xlab = "",
          ylab = "",
          xlim = c(0, 5),
          ylim = c(0, 5),
          main = "Plot Will Appear Here Once Data is Selected"
        ))
      }
      resultData <- reactive_DGE()
      colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      req(resultData)
      switch(
        selectedPlot,
        "MA plot" = {
          plotMA(resultData[[3]], cex = 0.7, ylim = c(-10, 10))
          abline(h = c(-1, 1), col = "red", lwd = 3)
        },
        "Shrunken LFC MA Plot" = {
          res_shrunk <- lfcShrink(resultData[[2]], coef = resultsNames(resultData[[2]])[[2]], type = "apeglm")
          plotMA(res_shrunk, cex = 0.7, ylim = c(-10, 10))
          abline(h = c(-1, 1), col = "red", lwd = 3)
        },
        "Dispersion plot" = plotDispEsts(resultData[[2]], main = "Dispersion Plot"),
        "PCA plot" = {
          rld <- rlogTransformation(resultData[[2]], blind = FALSE)
          plotPCA(rld, intgroup = "condition") +
            geom_text(aes(label = name), size = 3, vjust = 1.5, color = "#2C3E50") + 
            theme_seq_premium() +
            scale_color_brewer(palette = "Set1") +
            ggtitle("Sample Variation Analysis (PCA)")
        },
        "Volcano plot" = {
          # Aesthetic settings
          lfc_cut <- input$volcano_lfc
          pval_cut <- input$volcano_pval
          if(is.null(lfc_cut)) lfc_cut <- 1.0
          if(is.null(pval_cut)) pval_cut <- 0.05
          
          res_to_plot <- resultData[[1]]
          
          EnhancedVolcano(
            res_to_plot,
            lab = rownames(res_to_plot),
            x = "log2FoldChange",
            y = "padj",
            title = "Differential Expression Profile",
            subtitle = paste("Comparison:", resultsNames(resultData[[2]])[[2]]),
            pCutoff = pval_cut,
            FCcutoff = lfc_cut,
            pointSize = 2.5,
            labSize = 4.0,
            col = c("grey30", "forestgreen", "royalblue", "red2"),
            colAlpha = 0.6,
            legendPosition = "right",
            legendLabSize = 12,
            legendIconSize = 4.0,
            drawConnectors = TRUE,
            widthConnectors = 0.5
          )
        },
        "Heatmap of Pairwise Sample Distances" = {
          rld <- rlogTransformation(resultData[[2]], blind = FALSE)
          top_n <- input$heatmap_sample_top_n
          if(is.null(top_n)) top_n <- 500
          
          rv <- apply(assay(rld), 1, var)
          select <- order(rv, decreasing = TRUE)[1:min(top_n, length(rv))]
          
          dists <- dist(t(assay(rld)[select, ]))
          mat <- as.matrix(dists)
          
          # Premium Blue-White-Red Palette
          p_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
          
          pheatmap(
            mat,
            clustering_distance_rows = dists,
            clustering_distance_cols = dists,
            col = p_colors,
            main = "Sample Cluster Analysis (Distance Heatmap)",
            border_color = "white",
            cellwidth = 40, cellheight = 40,
            fontsize = 12
          )
        },
        "Heatmap of Z Scores for the Top Genes" = {
          top_n <- input$heatmap_top_n
          if(is.null(top_n)) top_n <- 20
          res_df <- resultData[[1]]
          res_df <- res_df[order(res_df$padj), ]
          n_genes <- min(top_n, nrow(res_df))
          top_genes <- rownames(res_df)[1:n_genes]
          
          mat <- counts(resultData[[2]], normalized = TRUE)
          mat <- mat[top_genes, , drop=FALSE]
          mat_z <- t(apply(mat, 1, cal_z_scores))
          
          meta <- resultData[[4]]
          
          # Intelligent scaling: Fix height at 800 for UI, adjust fonts
          # If n_genes > 150, hide names for a 'fingerprint' aesthetic
          show_rn <- n_genes <= 120
          fs <- ifelse(n_genes > 50, max(1, 400 / n_genes), 10)
          
          # Modern colors
          p_colors <- colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"))(100)
          
          pheatmap(
            mat_z,
            cluster_rows = TRUE,
            show_rownames = show_rn,
            cluster_cols = TRUE,
            annotation_col = meta,
            fontsize_row = fs,
            color = p_colors,
            main = paste("Expression Heatmap: Top", n_genes, "Genes"),
            border_color = NA
          )
        },
        "Summary Plot for DEGs" = {
          res_df <- resultData[[1]]
          lfc_cut <- input$summary_lfc
          pval_cut <- input$summary_pval
          if(is.null(lfc_cut)) lfc_cut <- 1.0
          if(is.null(pval_cut)) pval_cut <- 0.05
          
          up <- sum(res_df$log2FoldChange > lfc_cut & res_df$padj < pval_cut, na.rm = TRUE)
          down <- sum(res_df$log2FoldChange < -lfc_cut & res_df$padj < pval_cut, na.rm = TRUE)
          
          df_bar <- data.frame(
            Category = factor(c("Upregulated", "Downregulated"), levels = c("Upregulated", "Downregulated")),
            Count = c(up, down)
          )
          
          ggplot(df_bar, aes(x = Category, y = Count, fill = Category)) +
            geom_bar(stat = "identity", color="white", alpha=0.9, width = 0.35) +  # Sleek, narrow bars
            geom_text(aes(label = Count), vjust = -1, size = 6, face="bold", color="#2C3E50") +
            theme_seq_premium() +
            scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
            labs(title = "Differential Expression Summary",
                 subtitle = paste("Thresholds: |LFC| >", lfc_cut, ", padj <", pval_cut),
                 x = "", y = "Gene Count") +
            theme(legend.position = "none")
        }
      )
    }, width = 700, height = 430) 


    #-------------------------------------------------------------------------
    # DOWNLOAD THE PLOT :
    #-------------------------------------------------------------------------
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("Plot-", gsub(" ", "_", input$select_plot_type), "-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        req(input$select_plot_type)
        selectedPlot <- input$select_plot_type
        resultData <- reactive_DGE()
        colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
        
        # Default landscape dimensions (10x8 inch @ 300 DPI)
        img_w <- 10 * 300
        img_h <- 8 * 300
        
        if (!is.null(resultData) && selectedPlot != "None") {
          
          switch(
            selectedPlot,
            "MA plot" = {
              png(file, width = img_w, height = img_h, res = 300)
              plotMA(resultData[[3]], cex = 0.7, ylim = c(-10, 10))
              abline(h = c(-1, 1), col = "red", lwd = 3)
            },
            "Shrunken LFC MA Plot" = {
              png(file, width = img_w, height = img_h, res = 300)
              res_shrunk <- lfcShrink(resultData[[2]], coef = resultsNames(resultData[[2]])[[2]], type = "apeglm")
              plotMA(res_shrunk, cex = 0.7, ylim = c(-10, 10))
              abline(h = c(-1, 1), col = "red", lwd = 3)
            },
            "Dispersion plot" = {
              png(file, width = img_w, height = img_h, res = 300)
              plotDispEsts(resultData[[2]], main = "Dispersion Plot")
            },
            "PCA plot" = {
              png(file, width = img_w, height = img_h, res = 300)
              rld <- rlogTransformation(resultData[[2]], blind = FALSE)
              p <- plotPCA(rld, intgroup = "condition") +
                geom_text(aes(label = name), size = 3, vjust = 1.5, color = "#2C3E50") + 
                theme_seq_premium() +
                scale_color_brewer(palette = "Set1") +
                ggtitle("Variation Analysis (PCA)")
              print(p)
            },
            "Volcano plot" = {
              png(file, width = img_w, height = img_h, res = 300)
              lfc_cut <- input$volcano_lfc
              pval_cut <- input$volcano_pval
              if(is.null(lfc_cut)) lfc_cut <- 1.0
              if(is.null(pval_cut)) pval_cut <- 0.05
              p <- EnhancedVolcano(
                resultData[[1]],
                lab = rownames(resultData[[1]]),
                x = "log2FoldChange",
                y = "padj",
                title = "Differential Expression Profile",
                subtitle = paste("Comparison:", resultsNames(resultData[[2]])[[2]]),
                pCutoff = pval_cut,
                FCcutoff = lfc_cut,
                pointSize = 2.5,
                labSize = 4.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 0.6,
                legendPosition = "right",
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5
              )
              print(p)
            },
            "Heatmap of Pairwise Sample Distances" = {
              rld <- rlogTransformation(resultData[[2]], blind = FALSE)
              top_n <- input$heatmap_sample_top_n
              if(is.null(top_n)) top_n <- 500
              
              rv <- apply(assay(rld), 1, var)
              select <- order(rv, decreasing = TRUE)[1:min(top_n, length(rv))]
              
              dists <- dist(t(assay(rld)[select, ]))
              mat <- as.matrix(dists)
              
              p_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
              
              png(file, width = img_w, height = img_h, res = 300)
              pheatmap(
                mat,
                clustering_distance_rows = dists,
                clustering_distance_cols = dists,
                col = p_colors,
                header_rel = 0.2,
                main = "Sample Cluster Analysis",
                border_color = "white"
              )
            },
            "Heatmap of Z Scores for the Top Genes" = {
              top_n <- input$heatmap_top_n
              if(is.null(top_n)) top_n <- 20
              res_df <- resultData[[1]]
              res_df <- res_df[order(res_df$padj), ]
              n_genes <- min(top_n, nrow(res_df))
              top_genes <- rownames(res_df)[1:n_genes]
              
              mat <- counts(resultData[[2]], normalized = TRUE)
              mat <- mat[top_genes, , drop=FALSE]
              mat_z <- t(apply(mat, 1, cal_z_scores))
              
              meta <- resultData[[4]]
              
              show_rn <- n_genes <= 120
              # Precision Font Scaling: Use a conservative coefficient to prevent overlap in the fixed high-res area
              fs <- min(8, (img_h * 0.12) / n_genes) 
              if (n_genes > 120) show_rn <- FALSE
              
              p_colors <- colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"))(100)
              
              png(file, width = img_w, height = img_h, res = 300)
              pheatmap(
                mat_z,
                cluster_rows = TRUE,
                show_rownames = show_rn,
                cluster_cols = TRUE,
                annotation_col = meta,
                fontsize_row = fs,
                color = p_colors,
                main = paste("Top", n_genes, "Genes Heatmap"),
                border_color = NA
              )
            },
            "Summary Plot for DEGs" = {
              res_df <- resultData[[1]]
              lfc_cut <- input$summary_lfc
              pval_cut <- input$summary_pval
              if(is.null(lfc_cut)) lfc_cut <- 1.0
              if(is.null(pval_cut)) pval_cut <- 0.05
              
              up <- sum(res_df$log2FoldChange > lfc_cut & res_df$padj < pval_cut, na.rm = TRUE)
              down <- sum(res_df$log2FoldChange < -lfc_cut & res_df$padj < pval_cut, na.rm = TRUE)
              
              df_bar <- data.frame(
                Category = factor(c("Upregulated", "Downregulated"), levels = c("Upregulated", "Downregulated")),
                Count = c(up, down)
              )
              
              png(file, width = img_w, height = img_h, res = 300)
              p <- ggplot(df_bar, aes(x = Category, y = Count, fill = Category)) +
                geom_bar(stat = "identity", color="white", alpha=0.9, width = 0.35) +
                geom_text(aes(label = Count), vjust = -1, size = 6, face="bold", color="#2C3E50") +
                theme_seq_premium() +
                scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
                labs(title = "Differential Expression Summary",
                     subtitle = paste("Thresholds: |LFC| >", lfc_cut, ", padj <", pval_cut),
                     x = "", y = "Gene Count") +
                theme(legend.position = "none")
              print(p)
            }
          )
          dev.off()
        } else {
          png(file, width = img_w, height = img_h, res = 300)
          plot(1, type = "n", main = "No data available to plot")
          dev.off()
        }
      }
    )

    output$plot_desc <- renderText({
      selectedPlot <- input$select_plot_type
      if (selectedPlot == "None") {
        return("Setup your input data!")
      }
      switch(
        selectedPlot,
        "MA plot" = "it visualizes and indentifies gene expression changes from two different conditions (e.g.untreated vs. treated) in terms of log fold change (M)
                       on Y-axis and log of the mean of normalized expression counts of two conditions on X-axis.
                       Generally, genes with lower mean expression values will have higher variable log fold changes.
                       The blue dots are representing the differentially expressed genes and they have adjusted p-values of less than alpha parameter.
                       The triangles in the edge of the plot are genes with higher fold changes and the directions of the triangle are the direction of the fold change.
                       Genes on the right part of the plot means that these genes have high mean of normalized count
                       and high fold changes (They are very interesting to look into). " ,
        "Shrunken LFC MA Plot" = "The shrunken MA plot is similar to the MA plot, but it visualizes shrunken log2 fold changes, which remove noise associated
                               with log2 fold changes from low-count genes without requiring arbitrary filtering thresholds.",
        "Dispersion plot" = "Displays the dispersion estimates for the differential expression analysis.",
        "PCA plot" = "Generates a PCA plot based on the rlog-transformed counts, highlighting differences between conditions." ,
        "Volcano plot" = "Displays a volcano plot for differential expression analysis, highlighting genes with significant log2 fold changes and p-values.",
        "Heatmap of Pairwise Sample Distances" = "Displays a heatmap of the sample distances and clustering based on the log-transformed normalized counts.",
        "Heatmap of Z Scores for the Top Genes" = "Displays a heatmap of Z scores for the top selected genes",
        "Summary Plot for DEGs" = "Displays a bar plot summarizing the total, upregulated, and downregulated differentially expressed genes based on the configured cutoffs."
      )

    })
    #-------------------------------------------------------------------------
    # DOWNLOAD THE RESULT TABLE :
    #-------------------------------------------------------------------------
    output$btn_dwn_Degs <- downloadHandler(
      filename = function() {
        paste("DEGs-Analysis-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        resultData <- reactive_DGE()
        write.csv(as.matrix(resultData[[1]]),
                  file,
                  quote = F,
                  row.names = T)
      }
    )


  })
}
