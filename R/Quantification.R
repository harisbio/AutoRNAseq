quantificationUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      width = NULL,
      title = span(icon("calculator"), " Feature Counting with featureCounts"),
      status = "primary",
      solidHeader = TRUE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",

      sidebarLayout(
        sidebarPanel(
          style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",

          # 📂 Upload Annotation File
          h4("📂 Upload Annotation File", style = "color: #0092AC;"),
          fileInput(ns("dw_GTF"),
                    "Select an Annotation File in gtf format",
                    accept = c(".gtf", ".gtf.gz")),

          # # 🔍 Feature & Attribute Selection
          # h4("🔍 Select Feature & Attribute", style = "color: #0092AC;"),
          # selectInput(ns("select_feature_type"), "Select Feature Type:", choices = NULL),
          # selectInput(ns("select_attribute_type"), "Select Attribute Type:", choices = NULL),

          # 📂 BAM File Directory
          h4("📁 Mapping Files", style = "color: #0092AC;"),
          fileInput(
            ns("bam_upload"),
            "Upload BAM/SAM files",
            multiple = TRUE,
            accept = c(".bam", ".sam")
          ),
          tags$p("Preferred option: upload one or more BAM/SAM files directly.",
                 style = "color: #D6EAF8; font-weight:bold;"),
          br(),
          textInputIcon(
            ns("BAMFiles"),
            label = "Or provide a folder path for Mapping Files (BAM/SAM)",
            placeholder = "e.g., C:/path/to/mapping_files",
            icon = icon("folder-open")
          ),
          selectInput(
            ns("library_layout"),
            "Library Layout",
            choices = c("Paired-end" = "paired", "Single-end" = "single"),
            selected = "paired"
          ),
          tags$p("If using a folder path, the specified folder should contain only mapping files.",
                 style = "color:blue; font-weight:bold;"),
          hr(),

          # 🚀 Run Read Counting Button (Small)
          div(
            style = "display: flex; justify-content: left;",
            actionButton(ns("btn_count"),
                         "Run Read Counting",
                         icon = icon("play"),
                         style = "color: #ffffff; background-color: #0092AC;
                         border-color: #007B9E; padding: 6px 12px; font-size: 14px;
                         border-radius: 5px; width: auto;")
          ),
          br(),
          div(
            style = "display: flex; justify-content: left;",
            help_modal_button(ns("help_quant"), "Help")
          )
        ),

        mainPanel(
          tabBox(
            width = 12,

            # 📄 Annotation File Preview Tab
            tabPanel(
              title = span(icon("file"), " Annotation File Preview"),
              value = "tab-annotation",
              box(
                title = span(icon("database"), " Annotation File Overview"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("AnnotationFile"))
              )
            ),

            # 📊 Flattened Annotation Features Tab
            tabPanel(
              title = span(icon("align-left"), " Flattened Annotation Features"),
              value = "tab-flatten",
              box(
                title = span(icon("list"), " Processed Annotation Features"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("flatten_features"))
              )
            ),

            # 📊 FeatureCounts Results Tab
            tabPanel(
              title = span(icon("table"), " FeatureCounts Results"),
              value = "tab-featurecounts",
              box(
                title = span(icon("chart-bar"), " Count Matrix"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("featureCount"))
              ),
              # 📥 Download Count Matrix Button (Small)
              div(
                style = "display: flex; justify-content: left;",
                downloadButton(ns("dwn_count_matrix"),
                               "Download Count Matrix",
                               style = "color: #ffffff; background-color: #0092AC;
               border-color: #007B9E; padding: 6px 12px; font-size: 14px;
               border-radius: 5px; width: auto;")

              )
            ),
            tabPanel(
              title = span(icon("chart-area"), " Quantification Visualizations"),
              value = "tab-quant-viz",
              selectInput(
                ns("count_viz_type"),
                "Visualization Type",
                choices = c(
                  "Density Plot (Log10 Counts)" = "density",
                  "Library Size per Sample" = "library_size",
                  "Boxplot of Log2 Counts" = "boxplot",
                  "Sample PCA (Log2 Counts)" = "pca"
                ),
                selected = "density"
              ),
              box(
                title = span(icon("chart-bar"), " Quantification Plots"),
                width = NULL,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                plotOutput(ns("count_viz_plot"), height = "430px")
              )
            )
          ),

          br(),

          # ❓ Help Button (Small)
          div(
            style = "display: flex; justify-content: left;",
            actionButton(ns("help"),
                         "Help",
                         icon = icon("info-circle"),
                         style = "color: #ffffff; background-color: #F39C12;
                         border-color: #E67E22; padding: 6px 12px; font-size: 14px;
                         border-radius: 5px; width: auto;")
          )
        )
      )
    )
  )
}


# Quantification Server module
quantificationServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # Help modal for detailed instructions
      observeEvent(input$help_quant, {
        show_step_help(
          "Quantification Help",
          c(
            "Upload one annotation file in GTF, GTF.GZ, or GFF3 format.",
            "Use the same genome build that was used for alignment.",
            "Provide BAM/SAM files either by uploading them or by pointing to a folder path.",
            "Choose paired-end or single-end to match your alignment output.",
            "Start with a small test set before large batches so file and annotation mismatches are easier to spot."
          ),
          "featureCounts expects consistent read naming and a compatible annotation. If the first sample fails, check file paths, annotation format, and library layout."
        )
      })

      dataGTF <- reactive({
        if (is.null(input$dw_GTF$datapath)) {
          return(NULL)
        }
        GTF_File <- gsub("\\\\", "/", input$dw_GTF$datapath)
        # Import as GFF3 (change format to "gtf" if appropriate)
        Annotation <- rtracklayer::import(GTF_File, format = "gtf")
        return( as.data.frame(Annotation))
      })

      # Legacy selectors for feature/attribute columns were removed from the UI.
      # Keep the workflow focused on required inputs only (annotation + BAM/SAM files).

      # Render the annotation file preview
      output$AnnotationFile <- DT::renderDataTable({
        if (is.null(input$dw_GTF$datapath) || input$dw_GTF$datapath == "" ) {
          return(as.matrix("Please import an annotation file."))
        }
        DT::datatable(dataGTF(),
                      options = list(scrollX = TRUE, pageLength = 5),
                      caption = "Annotation File Details")
      })

      # Render the flattened annotation features table with error checking
      output$flatten_features <- DT::renderDataTable({
        if (is.null(input$dw_GTF$datapath) || input$dw_GTF$datapath == "") {
          return(as.matrix("Please import an annotation file."))
        }
        file <- gsub("\\\\", "/", input$dw_GTF$datapath)

          flattened <- flattenGTF(GTFfile = file)

          DT::datatable(as.data.frame(flattened),
                        options = list(scrollX = TRUE, pageLength = 5),
                        caption = "Flattened Features from Annotation")
      })

      mapping_files <- reactive({
        uploaded <- input$bam_upload
        if (!is.null(uploaded) && nrow(uploaded) > 0) {
          file_paths <- uploaded$datapath
          file_names <- uploaded$name
          valid_idx <- grepl("\\.(bam|sam)$", file_names, ignore.case = TRUE)
          file_paths <- file_paths[valid_idx]
          file_names <- file_names[valid_idx]
          if (length(file_paths) == 0) {
            showNotification("Uploaded files do not include valid BAM/SAM files.", type = "error")
            return(NULL)
          }
          names(file_paths) <- file_names
          return(file_paths)
        }

        bam_folder <- trimws(input$BAMFiles)
        if (is.null(bam_folder) || bam_folder == "") {
          return(NULL)
        }
        if (!dir.exists(bam_folder)) {
          showNotification("The provided mapping files folder does not exist.", type = "error")
          return(NULL)
        }

        file_paths <- list.files(
          path = bam_folder,
          pattern = "\\.(bam|sam)$",
          full.names = TRUE,
          ignore.case = TRUE
        )
        if (length(file_paths) == 0) {
          showNotification("No BAM/SAM files were found in the specified folder.", type = "error")
          return(NULL)
        }
        file_names <- basename(file_paths)
        names(file_paths) <- file_names
        file_paths
      })

      # Run featureCounts when "Run Read Counting" is clicked
      runCount <- reactive({
        req(input$dw_GTF$datapath)
        bam_files <- mapping_files()
        req(bam_files)
        if (length(bam_files) == 0) {
          showNotification("No BAM/SAM files are available for quantification.", type = "error")
          return(NULL)
        }
        AnnotationFile <- input$dw_GTF$datapath

        showNotification("Running read counting...", type = "message")
        countMatrix <- tryCatch(
          {
            fc_SE <- featureCounts(
              unname(bam_files),
              isGTFAnnotationFile = TRUE,
              annot.ext = AnnotationFile,
              isPairedEnd = identical(input$library_layout, "paired")
            )
            out <- as.data.frame(fc_SE$counts)
            output_names <- names(bam_files)
            output_names <- stringr::str_remove(output_names, "\\.(bam|sam)$")
            output_names <- stringr::str_remove(output_names, "-subread$")
            if (length(output_names) == ncol(out)) {
              colnames(out) <- output_names
            }
            out
          },
          error = function(e) {
            showNotification(paste("Read counting failed:", conditionMessage(e)), type = "error")
            NULL
          }
        )
        if (!is.null(countMatrix)) {
          showNotification("Read counting completed successfully.", type = "message")
        }
        return(countMatrix)
      })

      reactive_count <- reactiveVal(NULL)

      observeEvent(input$btn_count, {
        count_result <- runCount()
        if (is.null(count_result)) {
          showNotification("FeatureCounts did not run. Please verify annotation file, BAM/SAM inputs, and library layout.", type = "error")
        } else {
          showNotification("Feature count completed successfully.", type = "message")
        }
        reactive_count(count_result)
      })

      draw_quant_placeholder <- function(msg) {
        plot(
          x = 0, y = 0, type = "n",
          xlim = c(0, 1), ylim = c(0, 1),
          xlab = "", ylab = "", axes = FALSE,
          main = msg
        )
      }

      count_for_plot <- reactive({
        df <- reactive_count()
        if (is.null(df) || nrow(df) == 0) return(NULL)
        num_df <- as.data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))), stringsAsFactors = FALSE)
        num_df <- num_df[, colSums(!is.na(num_df)) > 0, drop = FALSE]
        if (ncol(num_df) == 0) return(NULL)
        num_df
      })

      output$count_viz_plot <- renderPlot({
        cdf <- count_for_plot()
        if (is.null(cdf) || nrow(cdf) == 0) {
          draw_quant_placeholder("Run featureCounts to view quantification visualizations.")
          return()
        }

        viz <- input$count_viz_type

        if (identical(viz, "density")) {
          vals <- log10(cdf + 1)
          long_df <- stack(vals)
          names(long_df) <- c("value", "sample")
          p <- ggplot2::ggplot(long_df, ggplot2::aes(x = value, color = sample)) +
            ggplot2::geom_density(linewidth = 0.9) +
            ggplot2::labs(title = "Density Plot of Log10 Counts", x = "log10(Count + 1)", y = "Density") +
            ggplot2::theme_minimal(base_size = 13)
          print(p)
          return()
        }

        if (identical(viz, "library_size")) {
          lib <- colSums(cdf, na.rm = TRUE)
          lib_df <- data.frame(sample = names(lib), library_size = as.numeric(lib), stringsAsFactors = FALSE)
          p <- ggplot2::ggplot(lib_df, ggplot2::aes(x = sample, y = library_size, fill = library_size)) +
            ggplot2::geom_col(alpha = 0.95) +
            ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey40") +
            ggplot2::labs(title = "Library Size per Sample", x = "Sample", y = "Total Assigned Reads") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        if (identical(viz, "boxplot")) {
          vals <- log2(cdf + 1)
          long_df <- stack(vals)
          names(long_df) <- c("value", "sample")
          p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = value, fill = sample)) +
            ggplot2::geom_boxplot(alpha = 0.9, outlier.alpha = 0.5) +
            ggplot2::labs(title = "Boxplot of Log2 Counts", x = "Sample", y = "log2(Count + 1)") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
          print(p)
          return()
        }

        mat <- log2(as.matrix(cdf) + 1)
        if (ncol(mat) < 2) {
          draw_quant_placeholder("At least two samples are required for PCA.")
          return()
        }
        pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
        pca_df <- data.frame(
          sample = rownames(pca$x),
          PC1 = pca$x[, 1],
          PC2 = pca$x[, 2],
          stringsAsFactors = FALSE
        )
        p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = sample)) +
          ggplot2::geom_point(size = 3, color = "#00B4D8") +
          ggplot2::geom_text(vjust = -0.7, size = 4, color = "#EEEEEE") +
          ggplot2::labs(title = "Sample PCA (Log2 Counts)", x = "PC1", y = "PC2") +
          ggplot2::theme_minimal(base_size = 13)
        print(p)
      })


      # Render the featureCounts count matrix table, limited to 15 rows per page
      output$featureCount <- DT::renderDataTable({
        if ((is.null(input$bam_upload) || nrow(input$bam_upload) == 0) &&
            (is.null(input$BAMFiles) || trimws(input$BAMFiles) == "")) {
          return(as.matrix("Please upload BAM/SAM files or provide a folder path that contains them."))
        }
        req(reactive_count())
        DT::datatable(reactive_count(),
                      options = list(pageLength = 5, scrollX = TRUE),
                      caption = "FeatureCounts: Count Matrix of Mapped Reads")
      })

      # Download handler for the count matrix
      output$dwn_count_matrix <- downloadHandler(
        filename = function() {
          paste("count-Matrix-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(reactive_count(), file, row.names = FALSE)
        }
      )
    }
  )
}
