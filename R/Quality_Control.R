# Quality control UI
qualityControlUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      width = NULL,
      title = span(icon("check-circle"), " Quality Control Assessment"),
      status = "primary",
      solidHeader = TRUE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
      sidebarLayout(
        sidebarPanel(
          style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
          h4("Select FASTQ Files", style = "color: #0092AC;"),
          selectInput(ns("select_fastq_file"), span(icon("file-upload"), " Choose FASTQ Files:"), choices = c("None"), multiple = TRUE),
          h4("Quality Control Tool", style = "color: #0092AC;"),
          selectInput(
            ns("qc_tool"),
            "Choose QC Tool",
            choices = c("Select QC Tool" = "none", "Rqc" = "rqc", "FastQC" = "fastqc", "MultiQC (FastQC outputs)" = "multiqc"),
            selected = "none"
          ),
          textInput(ns("qc_output_dir"), "QC Output Directory (optional)", placeholder = "Leave blank to auto-create inside FASTQ folder"),
          uiOutput(ns("qc_metric_ui")),
          actionButton(
            ns("run_qc_tool"),
            "Run Quality Control",
            icon = icon("play"),
            style = "color: #ffffff; background-color: #0092AC; border-color: #007B9E; padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;"
          ),
          hr(),
          downloadButton(
            ns("download_data"),
            "Download QC Report",
            style = "color: #ffffff; background-color: #0092AC; border-color: #007B9E; padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;"
          )
        ),
        mainPanel(
          tabBox(
            width = 12,
            tabPanel(
              title = span(icon("chart-bar"), " QC Visualizations"),
              uiOutput(ns("qc_viz_ui")),
              plotOutput(ns("quality_plot"), height = "420px"),
              br(),
              downloadButton(
                ns("download_plot"),
                "Download Plot",
                style = "color: #ffffff; background-color: #0092AC; border-color: #007B9E; padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;"
              )
            ),
            tabPanel(
              title = span(icon("table"), " QC Data Summary"),
              box(
                title = span(icon("database"), " Data Summary of QC Analysis"),
                width = NULL,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("quality_data"))
              ),
              br(),
              downloadButton(
                ns("download_data_csv"),
                "Download CSV",
                style = "color: #ffffff; background-color: #0092AC; border-color: #007B9E; padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;"
              )
            ),
            tabPanel(
              title = span(icon("tools"), " QC Run Details"),
              box(
                title = span(icon("info-circle"), " Tool Execution Status"),
                width = NULL,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                textOutput(ns("qc_tool_status")),
                br(),
                DT::dataTableOutput(ns("qc_tool_results"))
              ),
              box(
                title = span(icon("terminal"), " Tool Logs"),
                width = NULL,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                verbatimTextOutput(ns("qc_tool_log"))
              )
            )
          ),
          br(),
          actionButton(
            ns("help"),
            "Help",
            icon = icon("info-circle"),
            style = "color: #ffffff; background-color: #F39C12; border-color: #E67E22; padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;"
          )
        )
      )
    )
  )
}

# Server Module
qualityControlServer <- function(id, importedData) {
  moduleServer(id, function(input, output, session) {
    qc_tool_results <- reactiveVal(data.frame())
    qc_tool_log <- reactiveVal("No QC run yet.")
    qc_tool_status <- reactiveVal("Select FASTQ files, choose a QC tool, then run quality control.")
    qc_run_ready <- reactiveVal(FALSE)
    qc_last_tool <- reactiveVal("none")
    fastqc_module_df <- reactiveVal(data.frame())
    fastqc_basic_df <- reactiveVal(data.frame())
    multiqc_df <- reactiveVal(data.frame())

    select_files(session, "select_fastq_file", importedData, "Select Files")

    output$qc_metric_ui <- renderUI({
      if (!identical(input$qc_tool, "rqc")) return(NULL)
      tagList(
        h4("Rqc Metric Selection", style = "color: #0092AC;"),
        selectInput(
          session$ns("data_type"),
          span(icon("chart-line"), " Choose Rqc Metric:"),
          choices = c(
            "None",
            "Average Quality",
            "Cycle-specific Average Quality",
            "Cycle-specific Base Call Proportion",
            "Cycle-specific GC Content",
            "Cycle-specific Quality Distribution",
            "Cycle-specific Quality Distribution - Boxplot",
            "Per Read Mean Quality Distribution of Files",
            "Read Frequency",
            "Read Length Distribution"
          )
        )
      )
    })

    output$qc_viz_ui <- renderUI({
      if (identical(input$qc_tool, "rqc")) {
        return(selectInput(session$ns("qc_viz_type"), "Rqc Visualization Type", choices = c(
          "Rqc Metric Plot" = "rqc_native",
          "Metric Profile (Mean by Feature)" = "rqc_profile",
          "Metric Heatmap" = "rqc_heatmap"
        ), selected = "rqc_native"))
      }
      if (identical(input$qc_tool, "fastqc")) {
        return(selectInput(session$ns("qc_viz_type"), "FastQC Visualization Type", choices = c(
          "FastQC Module Status Distribution" = "fastqc_status_dist",
          "FastQC Module Status by Sample" = "fastqc_status_by_sample",
          "Total Sequences by Sample" = "fastqc_total_sequences",
          "GC Content by Sample" = "fastqc_gc"
        ), selected = "fastqc_status_dist"))
      }
      if (identical(input$qc_tool, "multiqc")) {
        return(selectInput(session$ns("qc_viz_type"), "MultiQC Visualization Type", choices = c(
          "Read Count by Sample" = "multiqc_reads",
          "GC Content by Sample" = "multiqc_gc",
          "Average Read Length by Sample" = "multiqc_length",
          "MultiQC Numeric Heatmap" = "multiqc_heatmap"
        ), selected = "multiqc_reads"))
      }
      helpText("Choose a QC tool to see tool-specific visualization options.")
    })

    draw_placeholder <- function(msg) {
      plot(x = 0, y = 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE, main = msg)
    }

    first_existing_col <- function(df, candidates) {
      for (nm in candidates) if (nm %in% colnames(df)) return(nm)
      NULL
    }

    resolve_input_paths <- function(selected_files) {
      selected <- selected_files[!selected_files %in% "None"]
      if (length(selected) == 0 || is.null(importedData())) return(NULL)
      path_info <- perFileInformation(importedData())$path
      if (is.null(path_info) || length(path_info) == 0) return(NULL)
      base_path <- path_info[1]
      full_paths <- file.path(base_path, selected)
      list(base_path = base_path, selected = selected, full_paths = full_paths)
    }

    parse_fastqc_zip <- function(zip_path) {
      if (!file.exists(zip_path)) return(list(module = data.frame(), basic = data.frame()))
      listing <- tryCatch(unzip(zip_path, list = TRUE)$Name, error = function(e) character(0))
      if (length(listing) == 0) return(list(module = data.frame(), basic = data.frame()))

      summary_entry <- listing[grepl("summary\\.txt$", listing)][1]
      data_entry <- listing[grepl("fastqc_data\\.txt$", listing)][1]
      module_df <- data.frame()
      basic_df <- data.frame()

      if (!is.na(summary_entry)) {
        sm <- tryCatch(readLines(unz(zip_path, summary_entry), warn = FALSE), error = function(e) character(0))
        if (length(sm) > 0) {
          module_rows <- lapply(strsplit(sm, "\t"), function(x) if (length(x) >= 2) data.frame(status = x[1], module = x[2], stringsAsFactors = FALSE))
          module_rows <- Filter(Negate(is.null), module_rows)
          if (length(module_rows) > 0) module_df <- do.call(rbind, module_rows)
        }
      }

      if (!is.na(data_entry)) {
        dl <- tryCatch(readLines(unz(zip_path, data_entry), warn = FALSE), error = function(e) character(0))
        if (length(dl) > 0) {
          s <- grep("^>>Basic Statistics", dl)
          e <- grep("^>>END_MODULE", dl)
          if (length(s) > 0 && length(e) > 0) {
            e <- e[e > s[1]][1]
            if (!is.na(e) && (s[1] + 2) <= (e - 1)) {
              sec <- dl[(s[1] + 2):(e - 1)]
              kv <- strsplit(sec, "\t")
              kv <- kv[vapply(kv, function(x) length(x) >= 2, logical(1))]
              if (length(kv) > 0) {
                vals <- setNames(lapply(kv, function(x) x[2]), vapply(kv, function(x) x[1], character(1)))
                basic_df <- as.data.frame(vals, stringsAsFactors = FALSE)
              }
            }
          }
        }
      }
      list(module = module_df, basic = basic_df)
    }

    path_inside_package <- function(path) {
      pkg_root <- normalizePath(system.file(package = "AutoRNAseq"), winslash = "/", mustWork = FALSE)
      target <- normalizePath(path, winslash = "/", mustWork = FALSE)
      nzchar(pkg_root) && startsWith(target, pkg_root)
    }

    qc_output_dir <- function(base_path, tool) {
      safe_tool <- gsub("[^A-Za-z0-9_-]+", "_", tool)
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      if (path_inside_package(base_path)) {
        root <- Sys.getenv(
          "AUTORNASEQ_QC_OUTPUT_PATH",
          unset = file.path(tempdir(), "AutoRNAseq", "QC_tools_output")
        )
      } else {
        root <- file.path(base_path, "QC_tools_output")
      }
      file.path(root, paste0(safe_tool, "_", stamp))
    }

    command_status <- function(result) {
      status_code <- attr(result, "status")
      if (!is.null(status_code) && !is.na(status_code) && status_code != 0) {
        return("Failed")
      }
      if (any(grepl("^ERROR:", result))) {
        return("Failed")
      }
      "Completed"
    }

    multiqc_env <- function() {
      ubuntu_python_lib <- "/usr/lib/python3/dist-packages"
      if (dir.exists(ubuntu_python_lib)) {
        return(paste0("PYTHONPATH=", ubuntu_python_lib))
      }
      character(0)
    }

    bind_result_rows <- function(...) {
      rows <- list(...)
      rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
      if (length(rows) == 0) {
        return(data.frame())
      }
      all_cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
      aligned <- lapply(rows, function(df) {
        missing_cols <- setdiff(all_cols, names(df))
        for (col in missing_cols) {
          df[[col]] <- NA_character_
        }
        df[, all_cols, drop = FALSE]
      })
      do.call(rbind, aligned)
    }

    run_fastqc <- function(input_paths, outdir) {
      rows <- list(); logs <- c(); module_rows <- list(); basic_rows <- list()
      for (i in seq_along(input_paths)) {
        file_path <- input_paths[i]
        result <- tryCatch(
          suppressWarnings(system2("fastqc", args = c(shQuote(file_path), "--outdir", shQuote(outdir)), stdout = TRUE, stderr = TRUE)),
          error = function(e) paste("ERROR:", conditionMessage(e))
        )
        sample_name <- tools::file_path_sans_ext(basename(file_path))
        zip_path <- file.path(outdir, paste0(sample_name, "_fastqc.zip"))
        html_path <- file.path(outdir, paste0(sample_name, "_fastqc.html"))
        parsed <- parse_fastqc_zip(zip_path)
        if (nrow(parsed$module) > 0) { tmp <- parsed$module; tmp$sample <- basename(file_path); module_rows[[length(module_rows) + 1]] <- tmp }
        if (nrow(parsed$basic) > 0) { tmp <- parsed$basic; tmp$sample <- basename(file_path); basic_rows[[length(basic_rows) + 1]] <- tmp }
        rows[[i]] <- data.frame(
          tool = "FastQC", input_file = basename(file_path),
          status = command_status(result),
          zip_report = zip_path, html_report = html_path, output_directory = outdir, stringsAsFactors = FALSE
        )
        if (identical(rows[[i]]$status, "Failed")) {
          status_code <- attr(result, "status")
          if (!is.null(status_code) && !is.na(status_code)) {
            result <- c(paste0("Exit status: ", status_code), result)
          }
        }
        logs <- c(logs, paste0("FastQC -> ", basename(file_path)), result, "")
      }
      list(
        df = do.call(bind_result_rows, rows), log = logs,
        module_df = if (length(module_rows) > 0) do.call(rbind, module_rows) else data.frame(),
        basic_df = if (length(basic_rows) > 0) do.call(rbind, basic_rows) else data.frame()
      )
    }

    parse_multiqc <- function(outdir) {
      mq_file <- file.path(outdir, "multiqc_data", "multiqc_fastqc.txt")
      if (!file.exists(mq_file)) return(data.frame())
      df <- tryCatch(utils::read.delim(mq_file, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame())
      if (nrow(df) == 0) return(df)
      sample_col <- first_existing_col(df, c("Sample", "Filename", "sample", "Name"))
      if (!is.null(sample_col)) df$sample <- as.character(df[[sample_col]]) else df$sample <- paste0("Sample_", seq_len(nrow(df)))
      for (nm in colnames(df)) {
        suppressWarnings(num <- as.numeric(df[[nm]]))
        if (!all(is.na(num))) df[[nm]] <- num
      }
      df
    }

    observeEvent(input$run_qc_tool, {
      qc_run_ready(FALSE); qc_last_tool("none")
      fastqc_module_df(data.frame()); fastqc_basic_df(data.frame()); multiqc_df(data.frame())

      paths <- resolve_input_paths(input$select_fastq_file)
      if (is.null(paths)) { qc_tool_status("No valid FASTQ files selected."); qc_tool_log("Please select one or more FASTQ files first."); qc_tool_results(data.frame()); return() }
      tool <- input$qc_tool
      if (is.null(tool) || tool == "none") { qc_tool_status("Please select a QC tool first."); qc_tool_log("No QC tool selected."); qc_tool_results(data.frame()); return() }

      if (tool == "rqc") {
        qc_tool_results(data.frame(tool = "Rqc", input_file = paste(paths$selected, collapse = ", "), status = "Completed", output_directory = "In-app visualizations", stringsAsFactors = FALSE))
        qc_tool_log("Rqc completed in-app."); qc_tool_status("Rqc run completed. Use Rqc visualizations."); qc_last_tool("rqc"); qc_run_ready(TRUE)
        return()
      }

      outdir <- trimws(input$qc_output_dir)
      if (outdir == "") outdir <- qc_output_dir(paths$base_path, tool)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      if (Sys.which("fastqc") == "") { qc_tool_status("FastQC is not installed on this system."); qc_tool_log("Install FastQC and re-run."); qc_tool_results(data.frame()); return() }

      withProgress(message = "Running selected QC tool...", value = 0.1, {
        f <- run_fastqc(paths$full_paths, outdir)
        fastqc_module_df(f$module_df); fastqc_basic_df(f$basic_df); incProgress(0.6)

        if (tool == "fastqc") {
          qc_tool_results(f$df); qc_tool_log(paste(f$log, collapse = "\n"))
          qc_tool_status(paste0("FastQC completed. Output directory: ", outdir)); qc_last_tool("fastqc"); qc_run_ready(TRUE); incProgress(0.3)
          return()
        }

        if (Sys.which("multiqc") == "") {
          qc_tool_status("MultiQC is not installed on this system."); qc_tool_results(f$df); qc_tool_log(paste(c(f$log, "MultiQC is required for this option."), collapse = "\n"))
          return()
        }

        mlog <- tryCatch(
          suppressWarnings(system2("multiqc", args = c(shQuote(outdir), "-o", shQuote(outdir), "--force"), stdout = TRUE, stderr = TRUE, env = multiqc_env())),
          error = function(e) paste("ERROR:", conditionMessage(e))
        )
        mdf <- parse_multiqc(outdir)
        multiqc_df(mdf)
        multiqc_status <- command_status(mlog)
        if (identical(multiqc_status, "Failed")) {
          status_code <- attr(mlog, "status")
          if (!is.null(status_code) && !is.na(status_code)) {
            mlog <- c(paste0("Exit status: ", status_code), mlog)
          }
        }
        mrow <- data.frame(tool = "MultiQC", input_file = basename(outdir), status = multiqc_status, output_directory = outdir, stringsAsFactors = FALSE)
        qc_tool_results(bind_result_rows(f$df, mrow)); qc_tool_log(paste(c(f$log, "MultiQC execution:", mlog), collapse = "\n"))
        if (identical(multiqc_status, "Completed")) {
          qc_tool_status(paste0("MultiQC completed on FastQC outputs. Output directory: ", outdir))
          qc_last_tool("multiqc")
          qc_run_ready(TRUE)
        } else {
          qc_tool_status(paste0("MultiQC failed, but the Shiny session is still stable. Check Tool Logs. FastQC outputs are in: ", outdir))
          qc_last_tool("fastqc")
          qc_run_ready(TRUE)
        }
        incProgress(0.3)
      })
    })

    draw_rqc_plot <- function(viz) {
      if (is.null(input$select_fastq_file) || "None" %in% input$select_fastq_file || is.null(importedData())) return(draw_placeholder("Select FASTQ files first."))
      metric <- input$data_type
      if (is.null(metric) || metric == "None") return(draw_placeholder("Choose an Rqc metric to plot."))
      if (identical(viz, "rqc_native")) return(print(quality_plots(metric, input$select_fastq_file, importedData)))
      raw_df <- datatable(metric, input$select_fastq_file, importedData)
      num_df <- raw_df[, vapply(raw_df, is.numeric, logical(1)), drop = FALSE]
      num_df <- num_df[, colSums(!is.na(num_df)) > 0, drop = FALSE]
      if (ncol(num_df) == 0) return(draw_placeholder("No numeric Rqc data available for this view."))
      if (identical(viz, "rqc_profile")) {
        p <- data.frame(feature = names(colMeans(num_df, na.rm = TRUE)), mean_value = as.numeric(colMeans(num_df, na.rm = TRUE)), stringsAsFactors = FALSE)
        return(print(ggplot2::ggplot(p, ggplot2::aes(x = feature, y = mean_value, group = 1)) + ggplot2::geom_line(color = "#00B4D8", linewidth = 1) + ggplot2::geom_point(color = "#90E0EF", size = 2.2) + ggplot2::labs(title = "Rqc Metric Profile", x = "Feature", y = "Mean Value") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
      }
      mat <- as.matrix(num_df); if (nrow(mat) > 200) mat <- mat[seq_len(200), , drop = FALSE]
      ld <- data.frame(row_id = rep(seq_len(nrow(mat)), times = ncol(mat)), feature = rep(colnames(mat), each = nrow(mat)), value = as.vector(mat), stringsAsFactors = FALSE)
      print(ggplot2::ggplot(ld, ggplot2::aes(x = feature, y = row_id, fill = value)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey40") + ggplot2::labs(title = "Rqc Metric Heatmap", x = "Feature", y = "Observation", fill = "Value") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
    }

    draw_fastqc_plot <- function(viz) {
      mdf <- fastqc_module_df(); bdf <- fastqc_basic_df()
      if (identical(viz, "fastqc_status_dist")) {
        if (nrow(mdf) == 0) return(draw_placeholder("No FastQC module summary found. Run FastQC first."))
        s <- as.data.frame(table(mdf$status), stringsAsFactors = FALSE); names(s) <- c("status", "count")
        return(print(ggplot2::ggplot(s, ggplot2::aes(x = status, y = count, fill = status)) + ggplot2::geom_col(width = 0.7, alpha = 0.9) + ggplot2::labs(title = "FastQC Module Status Distribution", x = "Status", y = "Count") + ggplot2::theme_minimal(base_size = 13)))
      }
      if (identical(viz, "fastqc_status_by_sample")) {
        if (nrow(mdf) == 0) return(draw_placeholder("No FastQC module summary found. Run FastQC first."))
        s <- as.data.frame(table(mdf$sample, mdf$status), stringsAsFactors = FALSE); names(s) <- c("sample", "status", "count")
        return(print(ggplot2::ggplot(s, ggplot2::aes(x = sample, y = count, fill = status)) + ggplot2::geom_col(position = "stack", alpha = 0.9) + ggplot2::labs(title = "FastQC Status by Sample", x = "Sample", y = "Module Count") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
      }
      if (nrow(bdf) == 0) return(draw_placeholder("No FastQC basic statistics found. Run FastQC first."))
      seq_col <- first_existing_col(bdf, c("Total Sequences", "Total.Sequences"))
      gc_col <- first_existing_col(bdf, c("%GC", "X.GC", "GC"))
      if (identical(viz, "fastqc_total_sequences")) {
        if (is.null(seq_col)) return(draw_placeholder("FastQC Total Sequences column not found."))
        bdf$seqs <- as.numeric(bdf[[seq_col]])
        return(print(ggplot2::ggplot(bdf, ggplot2::aes(x = sample, y = seqs, fill = sample)) + ggplot2::geom_col(alpha = 0.9, show.legend = FALSE) + ggplot2::labs(title = "FastQC Total Sequences by Sample", x = "Sample", y = "Total Sequences") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
      }
      if (is.null(gc_col)) return(draw_placeholder("FastQC GC column not found."))
      bdf$gc <- as.numeric(bdf[[gc_col]])
      print(ggplot2::ggplot(bdf, ggplot2::aes(x = sample, y = gc, fill = gc)) + ggplot2::geom_col(alpha = 0.9) + ggplot2::scale_fill_viridis_c(option = "B", na.value = "grey40") + ggplot2::labs(title = "FastQC GC Content by Sample", x = "Sample", y = "GC (%)") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
    }

    draw_multiqc_plot <- function(viz) {
      mdf <- multiqc_df()
      if (nrow(mdf) == 0) return(draw_placeholder("No MultiQC data found. Run MultiQC first."))
      sample_col <- first_existing_col(mdf, c("sample", "Sample", "Filename")); if (is.null(sample_col)) { mdf$sample <- paste0("Sample_", seq_len(nrow(mdf))); sample_col <- "sample" }
      reads_col <- first_existing_col(mdf, c("Total Sequences", "total_sequences", "total_sequences_fastqc"))
      gc_col <- first_existing_col(mdf, c("%GC", "gc_content", "percent_gc"))
      len_col <- first_existing_col(mdf, c("avg_sequence_length", "Sequence length", "sequence_length"))

      if (identical(viz, "multiqc_reads")) {
        if (is.null(reads_col)) return(draw_placeholder("MultiQC read count column not available."))
        mdf$reads <- suppressWarnings(as.numeric(mdf[[reads_col]]))
        return(print(ggplot2::ggplot(mdf, ggplot2::aes(x = .data[[sample_col]], y = reads, fill = .data[[sample_col]])) + ggplot2::geom_col(alpha = 0.9, show.legend = FALSE) + ggplot2::labs(title = "MultiQC Read Count by Sample", x = "Sample", y = "Reads") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
      }
      if (identical(viz, "multiqc_gc")) {
        if (is.null(gc_col)) return(draw_placeholder("MultiQC GC column not available."))
        mdf$gc <- suppressWarnings(as.numeric(mdf[[gc_col]]))
        return(print(ggplot2::ggplot(mdf, ggplot2::aes(x = .data[[sample_col]], y = gc, fill = gc)) + ggplot2::geom_col(alpha = 0.9) + ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey40") + ggplot2::labs(title = "MultiQC GC Content by Sample", x = "Sample", y = "GC (%)") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
      }
      if (identical(viz, "multiqc_length")) {
        if (is.null(len_col)) return(draw_placeholder("MultiQC read length column not available."))
        mdf$len <- suppressWarnings(as.numeric(mdf[[len_col]]))
        return(print(ggplot2::ggplot(mdf, ggplot2::aes(x = .data[[sample_col]], y = len, fill = .data[[sample_col]])) + ggplot2::geom_col(alpha = 0.9, show.legend = FALSE) + ggplot2::labs(title = "MultiQC Average Read Length by Sample", x = "Sample", y = "Read Length") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
      }
      num_df <- mdf[, vapply(mdf, is.numeric, logical(1)), drop = FALSE]; num_df <- num_df[, colSums(!is.na(num_df)) > 0, drop = FALSE]
      if (ncol(num_df) == 0) return(draw_placeholder("No numeric MultiQC columns available for heatmap."))
      if (ncol(num_df) > 20) num_df <- num_df[, seq_len(20), drop = FALSE]
      mat <- as.matrix(num_df)
      ld <- data.frame(row_id = rep(seq_len(nrow(mat)), times = ncol(mat)), feature = rep(colnames(mat), each = nrow(mat)), value = as.vector(mat), stringsAsFactors = FALSE)
      print(ggplot2::ggplot(ld, ggplot2::aes(x = feature, y = row_id, fill = value)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "D", na.value = "grey40") + ggplot2::labs(title = "MultiQC Numeric Heatmap", x = "Metric", y = "Sample Index", fill = "Value") + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
    }

    draw_qc_plot <- function() {
      if (!isTRUE(qc_run_ready())) return(draw_placeholder("Run quality control to view tool-specific visualizations."))
      if (identical(qc_last_tool(), "rqc")) return(draw_rqc_plot(input$qc_viz_type))
      if (identical(qc_last_tool(), "fastqc")) return(draw_fastqc_plot(input$qc_viz_type))
      if (identical(qc_last_tool(), "multiqc")) return(draw_multiqc_plot(input$qc_viz_type))
      draw_placeholder("No QC result available.")
    }

    output$quality_plot <- renderPlot({ draw_qc_plot() })

    current_qc_table <- reactive({
      if (!isTRUE(qc_run_ready())) return(data.frame(Message = "Run quality control first."))
      if (identical(qc_last_tool(), "rqc")) {
        if (is.null(input$data_type) || input$data_type == "None" || is.null(importedData()) || is.null(input$select_fastq_file) || "None" %in% input$select_fastq_file) {
          return(data.frame(Message = "For Rqc, choose files and a metric."))
        }
        return(datatable(input$data_type, input$select_fastq_file, importedData))
      }
      if (identical(qc_last_tool(), "fastqc")) {
        if (nrow(fastqc_basic_df()) > 0) return(fastqc_basic_df())
        if (nrow(fastqc_module_df()) > 0) return(fastqc_module_df())
        return(data.frame(Message = "FastQC summary is not available."))
      }
      if (identical(qc_last_tool(), "multiqc")) {
        if (nrow(multiqc_df()) > 0) return(multiqc_df())
        return(data.frame(Message = "MultiQC summary is not available."))
      }
      data.frame(Message = "No QC data available.")
    })

    output$quality_data <- DT::renderDataTable({ DT::datatable(current_qc_table(), options = list(scrollX = TRUE, pageLength = 8)) })
    output$qc_tool_status <- renderText({ qc_tool_status() })
    output$qc_tool_log <- renderText({ qc_tool_log() })
    output$qc_tool_results <- DT::renderDataTable({ DT::datatable(qc_tool_results(), options = list(scrollX = TRUE, pageLength = 8)) })

    output$download_plot <- downloadHandler(
      filename = function() paste("quality_plot_", Sys.Date(), ".png", sep = ""),
      content = function(file) { png(file, width = 1200, height = 650); draw_qc_plot(); dev.off() }
    )
    output$download_data_csv <- downloadHandler(
      filename = function() paste("quality_data_", Sys.Date(), ".csv", sep = ""),
      content = function(file) write.csv(current_qc_table(), file, row.names = FALSE)
    )
    output$download_data <- downloadHandler(
      filename = function() paste("quality_data_", Sys.Date(), ".csv", sep = ""),
      content = function(file) write.csv(current_qc_table(), file, row.names = FALSE)
    )

    observeEvent(input$help, {
      show_step_help(
        "Quality Control Help",
        c(
          "Select one or more FASTQ or FASTQ.GZ files that were imported through Data Setup.",
          "Choose Rqc for built-in read quality summaries, FastQC for per-file QC reports, or MultiQC to summarize FastQC outputs.",
          "Paired-end data should keep both mates available and named consistently so the QC tools can match them correctly.",
          "This module is designed for raw sequencing reads, not BAM or count matrices.",
          "For a smooth first run, start with one or two small FASTQ files before using a full project.",
          "Very large FASTQ files may take time and disk space, especially when FastQC or MultiQC generates reports."
        ),
        extra = "If the tool does not appear to react, check that the selected files exist in the imported dataset and that the QC software is available inside the container."
      )
    })
  })
}
