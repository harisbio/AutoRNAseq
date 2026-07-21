
Filter_TrimUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      width = NULL,
      title = span(icon("filter"), " Filtering & Trimming Reads"),
      status = "primary",
      solidHeader = TRUE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
      sidebarLayout(
        sidebarPanel(
          style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
          h4("Filtering & Trimming Tool Selection", style = "color: #0092AC;"),
          selectInput(
            ns("trimming_tool"),
            "Select Trimming Engine",
            choices = c(
              "Select Trimming Tool" = "none",
              "QuasR preprocessReads" = "quasr",
              "fastp" = "fastp",
              "Trimmomatic" = "trimmomatic"
            ),
            selected = "none"
          ),
          textInput(
            ns("trim_output_dir"),
            "Tool Output Directory (optional)",
            placeholder = "Leave blank to auto-create inside FASTQ folder"
          ),
          hr(),
          h4("Set Filtering Parameters", style = "color: #0092AC;"),
          checkboxInput(ns("pair_single"), "Paired-End Sequencing", FALSE),
          textInput(ns("truncateStartBases"), "Bases to Remove (Start):", placeholder = "Start Bases"),
          textInput(ns("truncateEndBases"), "Bases to Remove (End):", placeholder = "End Bases"),
          textInput(ns("minLength"), "Minimum Sequence Length:", placeholder = "Min Length"),
          uiOutput(ns("adapters_remover")),
          textInput(ns("nBases"), "Max 'N' Bases Allowed:", placeholder = "Max Ns"),
          textInput(ns("complexity"), "Minimum Sequence Complexity:", placeholder = "Complexity"),
          uiOutput(ns("trimming_tool_options")),
          hr(),
          actionButton(
            ns("btn_filter"),
            "Run Filtering & Trimming",
            icon = icon("play"),
            style = "color: #ffffff; background-color: #0092AC;
                                 border-color: #007B9E; padding: 6px 12px; font-size: 14px;
                                 border-radius: 5px; width: auto;"
          ),
          br(),
          div(
            style = "display: flex; justify-content: left;",
            help_modal_button(ns("help_trim"), "Help")
          )
        ),
        mainPanel(
          tabBox(
            width = 12,
            tabPanel(
              title = span(icon("chart-bar"), " Preprocessing Results"),
              value = "tab-preprocessing",
              h4("Select FASTQ Files", style = "color: #0092AC;"),
              selectInput(ns("inputfile"), "Choose FASTQ Files:", choices = c(), multiple = TRUE),
              uiOutput(ns("filemate")),
              helpText("Multiple files can be processed, but all sequence file vectors must have identical lengths."),
              hr(),
              box(
                title = span(icon("database"), " Filtered & Trimmed Data Summary"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("tb_filter_trim"))
              ),
              div(
                style = "display: flex; justify-content: left;",
                downloadButton(
                  ns("dwn_filter_report"),
                  "Download",
                  style = "color: #ffffff; background-color: #0092AC;
                                           border-color: #007B9E; padding: 6px 12px; font-size: 14px;
                                           border-radius: 5px; width: auto;"
                )
              )
            ),
            tabPanel(
              title = span(icon("terminal"), " Tool Logs"),
              value = "tab-trim-logs",
              box(
                title = span(icon("info-circle"), " Execution Status"),
                width = NULL,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                textOutput(ns("trim_tool_status"))
              ),
              box(
                title = span(icon("terminal"), " Command Logs"),
                width = NULL,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                verbatimTextOutput(ns("trim_tool_log"))
              )
            ),
            tabPanel(
              title = span(icon("chart-column"), " Trimming Visualizations"),
              value = "tab-trim-viz",
              uiOutput(ns("trim_viz_ui")),
              box(
                title = span(icon("chart-bar"), " Trimming Output Visualizations"),
                width = NULL,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                plotOutput(ns("trim_plot"), height = "420px")
              )
            )
          ),
        )
      )
    )
  )
}


Filter_TrimServer <- function(id, importedData) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- NS(id)

      filtered_data <- reactiveVal(NULL)
      trim_tool_status <- reactiveVal("No trimming run started.")
      trim_tool_log <- reactiveVal("No command log available yet.")
      trim_plot_data <- reactiveVal(data.frame())
      trim_fastp_summary <- reactiveVal(data.frame())
      trim_last_tool <- reactiveVal("none")

      select_files(session, "inputfile", importedData, "Select Files")

      observeEvent(input$pair_single, {
        if (input$pair_single) {
          output$filemate <- renderUI({
            selectInput(ns("filenameMate"), "Sequence pairs from paired-end experiments", choices = c(), multiple = TRUE)
          })
          select_files(session, "filenameMate", importedData, "Sequence pairs from paired-end experiments")
        } else {
          output$filemate <- renderUI({ NULL })
        }
      })

      output$adapters_remover <- renderUI({
        if (input$trimming_tool %in% c("quasr", "fastp")) {
          list(
            textInput(ns("Lpattern"), "Left Adapter (5'-end):", placeholder = "Left Adapter"),
            textInput(ns("Rpattern"), "Right Adapter (3'-end):", placeholder = "Right Adapter")
          )
        } else {
          NULL
        }
      })

      output$trimming_tool_options <- renderUI({
        if (input$trimming_tool == "trimmomatic") {
          textInput(
            ns("trimmomatic_adapter_fa"),
            "Trimmomatic Adapter FASTA (optional)",
            placeholder = "Absolute path to adapters.fa"
          )
        } else {
          NULL
        }
      })

      output$trim_viz_ui <- renderUI({
        tool <- input$trimming_tool
        if (is.null(tool) || tool == "none") {
          return(helpText("Choose a trimming tool to see tool-specific visualization options."))
        }

        common_choices <- c(
          "Before vs After (Grouped Bars)" = "bars",
          "Size Reduction (%) by Sample" = "reduction_pct",
          "Sample-wise Before/After Trend" = "slope",
          "Retention Ratio (After/Before)" = "ratio",
          "Absolute Size Change (MB)" = "delta_mb",
          "Tool Run Status Summary" = "status_summary"
        )

        if (tool == "fastp") {
          choices <- c(
            common_choices,
            "fastp: Reads Before vs After" = "fastp_reads",
            "fastp: Q30 and GC Improvement" = "fastp_quality"
          )
          return(selectInput(ns("trim_viz_type"), "fastp Visualization Type", choices = choices, selected = "bars"))
        }

        label <- if (tool == "trimmomatic") "Trimmomatic Visualization Type" else "QuasR Visualization Type"
        selectInput(ns("trim_viz_type"), label, choices = common_choices, selected = "bars")
      })

      resolve_paths <- function() {
        input_file <- input$inputfile
        if (is.null(input_file) || length(input_file) == 0) {
          return(NULL)
        }
        path_info <- perFileInformation(importedData())$path
        if (is.null(path_info) || length(path_info) == 0) {
          return(NULL)
        }

        base_path <- path_info[1]
        input_files_path <- vapply(input_file, function(file) file.path(base_path, file), FUN.VALUE = character(1))
        mate_files_path <- character(0)

        if (isTRUE(input$pair_single)) {
          filename_mate <- input$filenameMate
          if (is.null(filename_mate) || length(filename_mate) == 0) {
            return(NULL)
          }
          mate_files_path <- vapply(filename_mate, function(file) file.path(base_path, file), FUN.VALUE = character(1))
          if (length(mate_files_path) != length(input_files_path)) {
            return(NULL)
          }
        }

        list(base_path = base_path, input_files_path = input_files_path, mate_files_path = mate_files_path)
      }

      path_inside_package <- function(path) {
        pkg_root <- normalizePath(system.file(package = "AutoRNAseq"), winslash = "/", mustWork = FALSE)
        target <- normalizePath(path, winslash = "/", mustWork = FALSE)
        nzchar(pkg_root) && startsWith(target, pkg_root)
      }

      run_output_dir <- function(base_path, tool) {
        safe_tool <- gsub("[^A-Za-z0-9_-]+", "_", tool)
        stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        if (path_inside_package(base_path)) {
          root <- Sys.getenv(
            "AUTORNASEQ_TRIM_OUTPUT_PATH",
            unset = file.path(tempdir(), "AutoRNAseq", "FilteredOutput")
          )
        } else {
          root <- file.path(base_path, "FilteredOutput")
        }
        file.path(root, paste0(safe_tool, "_", stamp))
      }

      unique_output_path <- function(path) {
        if (!file.exists(path)) {
          return(path)
        }
        dir_name <- dirname(path)
        ext <- tools::file_ext(path)
        stem <- tools::file_path_sans_ext(basename(path))
        if (tolower(ext) == "gz" && grepl("\\.(fastq|fq)$", stem, ignore.case = TRUE)) {
          second_ext <- tools::file_ext(stem)
          stem <- tools::file_path_sans_ext(stem)
          ext <- paste0(second_ext, ".", ext)
        }
        suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
        file.path(dir_name, paste0(stem, "-", suffix, if (nzchar(ext)) paste0(".", ext) else ""))
      }

      numeric_or_null <- function(value) {
        parsed <- suppressWarnings(as.numeric(value))
        if (length(parsed) == 0 || is.na(parsed[1])) {
          return(NULL)
        }
        parsed[1]
      }

      numeric_or_default <- function(value, default) {
        parsed <- numeric_or_null(value)
        if (is.null(parsed)) {
          return(default)
        }
        parsed
      }

      text_or_empty <- function(value) {
        if (is.null(value) || length(value) == 0) {
          return("")
        }
        trimws(as.character(value[1]))
      }

      fastq_output_name <- function(path, suffix, extension = ".fastq.gz") {
        name <- basename(path)
        stem <- sub("\\.(fastq|fq)(\\.gz)?$", "", name, ignore.case = TRUE)
        paste0(stem, suffix, extension)
      }

      tool_status <- function(result) {
        status_code <- attr(result, "status")
        if (!is.null(status_code) && !is.na(status_code) && status_code != 0) {
          return("Failed")
        }
        if (any(grepl("^ERROR:", result))) {
          return("Failed")
        }
        "Completed"
      }

      run_fastp <- function(paths, output_dir) {
        if (Sys.which("fastp") == "") {
          stop("fastp is not installed on this system.")
        }

        rows <- list()
        logs <- c()
        in_files <- paths$input_files_path
        mate_files <- paths$mate_files_path

        for (i in seq_along(in_files)) {
          in1 <- in_files[i]
          sample_stem <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(in1), ignore.case = TRUE)
          out1 <- unique_output_path(file.path(output_dir, fastq_output_name(in1, "-fastp")))
          json_report <- unique_output_path(file.path(output_dir, paste0(sample_stem, "-fastp.json")))
          html_report <- unique_output_path(file.path(output_dir, paste0(sample_stem, "-fastp.html")))

          args <- c("-i", shQuote(in1), "-o", shQuote(out1), "--json", shQuote(json_report), "--html", shQuote(html_report))

          if (length(mate_files) > 0) {
            in2 <- mate_files[i]
            out2 <- unique_output_path(file.path(output_dir, fastq_output_name(in2, "-fastp")))
            args <- c(args, "-I", shQuote(in2), "-O", shQuote(out2))
          }

          trim_front <- numeric_or_null(input$truncateStartBases)
          trim_tail <- numeric_or_null(input$truncateEndBases)
          min_length <- numeric_or_null(input$minLength)
          max_n <- numeric_or_null(input$nBases)

          if (!is.null(trim_front) && trim_front > 0) {
            args <- c(args, "--trim_front1", as.character(trim_front))
          }
          if (!is.null(trim_tail) && trim_tail > 0) {
            args <- c(args, "--trim_tail1", as.character(trim_tail))
          }
          if (!is.null(min_length) && min_length > 0) {
            args <- c(args, "--length_required", as.character(min_length))
          }
          if (!is.null(max_n) && max_n >= 0) {
            args <- c(args, "--n_base_limit", as.character(max_n))
          }
          left_adapter <- text_or_empty(input$Lpattern)
          right_adapter <- text_or_empty(input$Rpattern)
          if (nzchar(left_adapter)) {
            args <- c(args, "--adapter_sequence", shQuote(left_adapter))
          }
          if (nzchar(right_adapter) && length(mate_files) > 0) {
            args <- c(args, "--adapter_sequence_r2", shQuote(right_adapter))
          }

          result <- tryCatch(
            suppressWarnings(system2("fastp", args = args, stdout = TRUE, stderr = TRUE)),
            error = function(e) paste("ERROR:", conditionMessage(e))
          )
          status <- tool_status(result)

          rows[[i]] <- data.frame(
            tool = "fastp",
            input_file = basename(in1),
            output_file = basename(out1),
            input_path = in1,
            output_path = out1,
            json_report = json_report,
            html_report = html_report,
            status = status,
            output_directory = output_dir,
            stringsAsFactors = FALSE
          )
          logs <- c(logs, paste0("fastp -> ", basename(in1)), result, "")
        }

        list(df = do.call(rbind, rows), log = logs)
      }

      run_trimmomatic <- function(paths, output_dir) {
        rows <- list()
        logs <- c()
        in_files <- paths$input_files_path
        mate_files <- paths$mate_files_path
        is_paired <- length(mate_files) > 0

        trimmomatic_cmd <- if (is_paired) "TrimmomaticPE" else "TrimmomaticSE"
        if (Sys.which(trimmomatic_cmd) == "") {
          stop(paste0("Trimmomatic executable '", trimmomatic_cmd, "' not found on PATH."))
        }

        for (i in seq_along(in_files)) {
          in1 <- in_files[i]
          out1 <- unique_output_path(file.path(output_dir, fastq_output_name(in1, "-trimmomatic")))

          if (is_paired) {
            in2 <- mate_files[i]
            out2 <- unique_output_path(file.path(output_dir, fastq_output_name(in2, "-trimmomatic")))
            out1_unpaired <- unique_output_path(file.path(output_dir, fastq_output_name(in1, "-unpaired")))
            out2_unpaired <- unique_output_path(file.path(output_dir, fastq_output_name(in2, "-unpaired")))
            args <- c("-threads", "1", shQuote(in1), shQuote(in2), shQuote(out1), shQuote(out1_unpaired), shQuote(out2), shQuote(out2_unpaired))
          } else {
            args <- c("-threads", "1", shQuote(in1), shQuote(out1))
          }

          trimmers <- character(0)
          trim_front <- numeric_or_null(input$truncateStartBases)
          min_length <- numeric_or_default(input$minLength, 14)

          if (!is.null(trim_front) && trim_front > 0) {
            trimmers <- c(trimmers, paste0("HEADCROP:", trim_front))
          }
          if (!is.null(min_length) && min_length > 0) {
            trimmers <- c(trimmers, paste0("MINLEN:", min_length))
          }
          adapter_file <- text_or_empty(input$trimmomatic_adapter_fa)
          if (nzchar(adapter_file)) {
            trimmers <- c(trimmers, paste0("ILLUMINACLIP:", adapter_file, ":2:30:10"))
          }

          args <- c(args, trimmers)

          result <- tryCatch(
            suppressWarnings(system2(trimmomatic_cmd, args = args, stdout = TRUE, stderr = TRUE)),
            error = function(e) paste("ERROR:", conditionMessage(e))
          )
          status <- tool_status(result)
          if (identical(status, "Failed")) {
            status_code <- attr(result, "status")
            if (!is.null(status_code) && !is.na(status_code)) {
              result <- c(paste0("Exit status: ", status_code), result)
            }
          }

          rows[[i]] <- data.frame(
            tool = "Trimmomatic",
            input_file = basename(in1),
            output_file = basename(out1),
            output_path = out1,
            status = status,
            output_directory = output_dir,
            stringsAsFactors = FALSE
          )
          logs <- c(logs, paste0("Trimmomatic -> ", basename(in1)), result, "")
        }

        list(df = do.call(rbind, rows), log = logs)
      }

      run_quasr <- function(paths, output_dir) {
        output_files <- vapply(
          file.path(output_dir, vapply(paths$input_files_path, fastq_output_name, suffix = "-Filtered", FUN.VALUE = character(1))),
          unique_output_path,
          FUN.VALUE = character(1)
        )
        left_adapter <- text_or_empty(input$Lpattern)
        right_adapter <- text_or_empty(input$Rpattern)
        preprocess_args <- list(
          filename = paths$input_files_path,
          outputFilename = output_files,
          nBases = numeric_or_default(input$nBases, 2),
          truncateStartBases = numeric_or_null(input$truncateStartBases),
          truncateEndBases = numeric_or_null(input$truncateEndBases),
          Lpattern = left_adapter,
          Rpattern = right_adapter,
          complexity = numeric_or_null(input$complexity),
          minLength = numeric_or_default(input$minLength, 14)
        )

        if (length(paths$mate_files_path) > 0) {
          output_mate_files <- vapply(
            file.path(output_dir, vapply(paths$mate_files_path, fastq_output_name, suffix = "-FilteredMate", FUN.VALUE = character(1))),
            unique_output_path,
            FUN.VALUE = character(1)
          )
          preprocess_args$filenameMate <- paths$mate_files_path
          preprocess_args$outputFilenameMate <- output_mate_files
        }

        result <- tryCatch(
          do.call(preprocessReads, preprocess_args),
          error = function(e) {
            stop(paste("QuasR preprocessReads failed:", conditionMessage(e)), call. = FALSE)
          }
        )
        missing_outputs <- output_files[!file.exists(output_files)]
        if (length(missing_outputs) > 0) {
          warning(
            "QuasR completed without creating expected output file(s): ",
            paste(basename(missing_outputs), collapse = ", "),
            call. = FALSE
          )
        }
        standardized <- data.frame(
          tool = "QuasR preprocessReads",
          input_file = basename(paths$input_files_path),
          output_file = basename(output_files),
          output_path = output_files,
          status = "Completed",
          output_directory = output_dir,
          stringsAsFactors = FALSE
        )
        list(df = standardized, log = c("QuasR preprocessReads completed.", capture.output(print(result))))
      }

      build_trim_plot_df <- function(paths, result_df) {
        if (is.null(result_df) || nrow(result_df) == 0 || !("output_path" %in% colnames(result_df))) {
          return(data.frame())
        }
        in_paths <- paths$input_files_path
        out_paths <- result_df$output_path
        n <- min(length(in_paths), length(out_paths))
        if (n == 0) {
          return(data.frame())
        }
        in_size <- file.info(in_paths[seq_len(n)])$size / (1024^2)
        out_size <- file.info(out_paths[seq_len(n)])$size / (1024^2)
        data.frame(
          sample = basename(in_paths[seq_len(n)]),
          before_mb = round(in_size, 3),
          after_mb = round(out_size, 3),
          stringsAsFactors = FALSE
        )
      }

      parse_fastp_reports <- function(result_df) {
        if (is.null(result_df) || nrow(result_df) == 0 || !("json_report" %in% colnames(result_df))) {
          return(data.frame())
        }
        if (!requireNamespace("jsonlite", quietly = TRUE)) {
          return(data.frame())
        }

        rows <- list()
        idx <- 1
        for (i in seq_len(nrow(result_df))) {
          jp <- result_df$json_report[i]
          if (is.na(jp) || !nzchar(jp) || !file.exists(jp)) {
            next
          }
          js <- tryCatch(jsonlite::fromJSON(jp), error = function(e) NULL)
          if (is.null(js) || is.null(js$summary)) {
            next
          }
          bf <- js$summary$before_filtering
          af <- js$summary$after_filtering
          rows[[idx]] <- data.frame(
            sample = result_df$input_file[i],
            reads_before = as.numeric(bf$total_reads %||% NA_real_),
            reads_after = as.numeric(af$total_reads %||% NA_real_),
            bases_before = as.numeric(bf$total_bases %||% NA_real_),
            bases_after = as.numeric(af$total_bases %||% NA_real_),
            q30_before = as.numeric(bf$q30_rate %||% NA_real_),
            q30_after = as.numeric(af$q30_rate %||% NA_real_),
            gc_before = as.numeric(bf$gc_content %||% NA_real_),
            gc_after = as.numeric(af$gc_content %||% NA_real_),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
        if (length(rows) == 0) return(data.frame())
        do.call(rbind, rows)
      }

      `%||%` <- function(x, y) if (is.null(x)) y else x

      observeEvent(input$help_trim, {
        show_step_help(
          "Filtering & Trimming Help",
          c(
            "Select one or more FASTQ/FASTQ.GZ files imported from Data Setup.",
            "Choose QuasR, fastp, or Trimmomatic depending on the preprocessing you want.",
            "If paired-end is enabled, make sure the mate files match sample-by-sample.",
            "Leave the output directory blank to create a default folder next to the input files.",
            "For very large files, start with one or two samples to validate the settings before batch processing."
          ),
          "Trimming is safest when the sample names, mate pairs, and file naming convention are consistent."
        )
      })

      observeEvent(input$btn_filter, {
        trim_fastp_summary(data.frame())
        trim_last_tool("none")
        if (is.null(importedData())) {
          showModal(
            modalDialog(
              title = "Data Import Required",
              "Please import your FASTQ files before proceeding.",
              easyClose = TRUE
            )
          )
          filtered_data(NULL)
          return()
        }

        paths <- resolve_paths()
        if (is.null(paths)) {
          showModal(modalDialog(title = "Input Error", "Please check selected files and paired-end mates.", easyClose = TRUE))
          return()
        }

        if (is.null(input$trimming_tool) || input$trimming_tool == "none") {
          showModal(modalDialog(title = "Tool Selection Required", "Please select a trimming tool first.", easyClose = TRUE))
          return()
        }

        output_dir <- trimws(input$trim_output_dir)
        if (output_dir == "") {
          output_dir <- run_output_dir(paths$base_path, input$trimming_tool)
        }
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        withProgress(message = "Processing data", value = 0.2, {
          res <- tryCatch({
            if (input$trimming_tool == "fastp") {
              run_fastp(paths, output_dir)
            } else if (input$trimming_tool == "trimmomatic") {
              run_trimmomatic(paths, output_dir)
            } else {
              run_quasr(paths, output_dir)
            }
          }, warning = function(w) {
            list(df = data.frame(status = "Warning", detail = conditionMessage(w), stringsAsFactors = FALSE), log = paste("WARNING:", conditionMessage(w)))
          }, error = function(e) {
            list(df = data.frame(status = "Error", detail = conditionMessage(e), stringsAsFactors = FALSE), log = paste("ERROR:", conditionMessage(e)))
          })

          filtered_data(res$df)
          trim_plot_data(build_trim_plot_df(paths, res$df))
          if (identical(input$trimming_tool, "fastp")) {
            trim_fastp_summary(parse_fastp_reports(res$df))
          } else {
            trim_fastp_summary(data.frame())
          }
          trim_last_tool(input$trimming_tool)
          trim_tool_log(paste(res$log, collapse = "\n"))
          trim_tool_status(paste0("Trimming engine '", input$trimming_tool, "' finished. Output directory: ", output_dir))
          incProgress(0.8)
        })
      })

      output$tb_filter_trim <- DT::renderDataTable({
        if (is.null(importedData())) {
          return(DT::datatable(data.frame(Message = "Import your FASTQ files")))
        }
        if (!is.null(filtered_data())) {
          return(DT::datatable(filtered_data(), options = list(scrollX = TRUE)))
        }
        DT::datatable(data.frame(Message = "Run filtering/trimming to view results."))
      })

      output$trim_tool_status <- renderText({
        trim_tool_status()
      })

      output$trim_tool_log <- renderText({
        trim_tool_log()
      })

      output$dwn_filter_report <- downloadHandler(
        filename = function() {
          paste("filter_report", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(filtered_data(), file, row.names = FALSE)
        }
      )

      draw_trim_placeholder <- function(msg) {
        plot(
          x = 0, y = 0, type = "n",
          xlim = c(0, 1), ylim = c(0, 1),
          xlab = "", ylab = "", axes = FALSE,
          main = msg
        )
      }

      output$trim_plot <- renderPlot({
        viz <- input$trim_viz_type
        df <- trim_plot_data()
        res_df <- filtered_data()
        fastp_df <- trim_fastp_summary()
        selected_tool <- trim_last_tool()

        if (is.null(viz) || selected_tool == "none") {
          draw_trim_placeholder("Select a trimming tool, run it, then choose a visualization.")
          return()
        }

        if (identical(viz, "status_summary")) {
          if (is.null(res_df) || nrow(res_df) == 0 || !("status" %in% names(res_df))) {
            draw_trim_placeholder("Run a trimming tool to view status summary.")
            return()
          }
          status_df <- as.data.frame(table(res_df$status), stringsAsFactors = FALSE)
          names(status_df) <- c("status", "count")
          p <- ggplot2::ggplot(status_df, ggplot2::aes(x = status, y = count, fill = status)) +
            ggplot2::geom_col(width = 0.7, alpha = 0.9) +
            ggplot2::labs(title = "Trimming Tool Status Summary", x = "Status", y = "Count") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(legend.position = "none")
          print(p)
          return()
        }

        if (identical(viz, "fastp_reads")) {
          if (!identical(selected_tool, "fastp")) {
            draw_trim_placeholder("This visualization is available after running fastp.")
            return()
          }
          if (is.null(fastp_df) || nrow(fastp_df) == 0) {
            draw_trim_placeholder("fastp JSON summary not found. Run fastp to view read-retention plots.")
            return()
          }
          long_df <- rbind(
            data.frame(sample = fastp_df$sample, stage = "Before", reads = fastp_df$reads_before, stringsAsFactors = FALSE),
            data.frame(sample = fastp_df$sample, stage = "After", reads = fastp_df$reads_after, stringsAsFactors = FALSE)
          )
          p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = reads, fill = stage)) +
            ggplot2::geom_col(position = "dodge", alpha = 0.9) +
            ggplot2::labs(title = "fastp Read Counts Before vs After", x = "Sample", y = "Reads") +
            ggplot2::scale_fill_manual(values = c("Before" = "#E67E22", "After" = "#2A9D8F")) +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        if (identical(viz, "fastp_quality")) {
          if (!identical(selected_tool, "fastp")) {
            draw_trim_placeholder("This visualization is available after running fastp.")
            return()
          }
          if (is.null(fastp_df) || nrow(fastp_df) == 0) {
            draw_trim_placeholder("fastp JSON summary not found. Run fastp to view Q30/GC improvements.")
            return()
          }
          quality_df <- rbind(
            data.frame(sample = fastp_df$sample, metric = "Q30", before = fastp_df$q30_before, after = fastp_df$q30_after, stringsAsFactors = FALSE),
            data.frame(sample = fastp_df$sample, metric = "GC", before = fastp_df$gc_before, after = fastp_df$gc_after, stringsAsFactors = FALSE)
          )
          long_df <- rbind(
            data.frame(sample = quality_df$sample, metric = quality_df$metric, stage = "Before", value = quality_df$before, stringsAsFactors = FALSE),
            data.frame(sample = quality_df$sample, metric = quality_df$metric, stage = "After", value = quality_df$after, stringsAsFactors = FALSE)
          )
          p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = value, fill = stage)) +
            ggplot2::geom_col(position = "dodge", alpha = 0.9) +
            ggplot2::facet_wrap(~ metric, scales = "free_y") +
            ggplot2::labs(title = "fastp Quality Metrics Before vs After", x = "Sample", y = "Metric Value") +
            ggplot2::scale_fill_manual(values = c("Before" = "#8D99AE", "After" = "#2A9D8F")) +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        if (is.null(df) || nrow(df) == 0) {
          draw_trim_placeholder("Run a trimming tool to view trimming visualizations.")
          return()
        }

        df$reduction_pct <- ifelse(df$before_mb > 0, ((df$before_mb - df$after_mb) / df$before_mb) * 100, NA_real_)
        df$retention_ratio <- ifelse(df$before_mb > 0, df$after_mb / df$before_mb, NA_real_)
        df$delta_mb <- df$before_mb - df$after_mb

        if (identical(viz, "bars")) {
          long_df <- rbind(
            data.frame(sample = df$sample, stage = "Before", size_mb = df$before_mb, stringsAsFactors = FALSE),
            data.frame(sample = df$sample, stage = "After", size_mb = df$after_mb, stringsAsFactors = FALSE)
          )
          p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = size_mb, fill = stage)) +
            ggplot2::geom_col(position = "dodge", alpha = 0.9) +
            ggplot2::labs(x = "Sample", y = "File Size (MB)", title = "Before vs After Trimming") +
            ggplot2::scale_fill_manual(values = c("Before" = "#E67E22", "After" = "#0092AC")) +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        if (identical(viz, "reduction_pct")) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = reduction_pct, fill = reduction_pct)) +
            ggplot2::geom_col(alpha = 0.95) +
            ggplot2::scale_fill_viridis_c(option = "B", na.value = "grey40") +
            ggplot2::labs(title = "Size Reduction Percentage by Sample", x = "Sample", y = "Reduction (%)") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        if (identical(viz, "slope")) {
          long_df <- rbind(
            data.frame(sample = df$sample, stage = "Before", size_mb = df$before_mb, stringsAsFactors = FALSE),
            data.frame(sample = df$sample, stage = "After", size_mb = df$after_mb, stringsAsFactors = FALSE)
          )
          long_df$stage <- factor(long_df$stage, levels = c("Before", "After"))
          p <- ggplot2::ggplot(long_df, ggplot2::aes(x = stage, y = size_mb, group = sample, color = sample)) +
            ggplot2::geom_line(linewidth = 1, alpha = 0.85) +
            ggplot2::geom_point(size = 2.6) +
            ggplot2::labs(title = "Sample-wise Trimming Trend", x = "", y = "File Size (MB)") +
            ggplot2::theme_minimal(base_size = 13)
          print(p)
          return()
        }

        if (identical(viz, "ratio")) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = retention_ratio, fill = retention_ratio)) +
            ggplot2::geom_col(alpha = 0.95) +
            ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#EEEEEE") +
            ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey40") +
            ggplot2::labs(title = "Retention Ratio After/Before", x = "Sample", y = "Retention Ratio") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        if (identical(viz, "delta_mb")) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = delta_mb, fill = delta_mb)) +
            ggplot2::geom_col(alpha = 0.95) +
            ggplot2::scale_fill_viridis_c(option = "A", na.value = "grey40") +
            ggplot2::labs(title = "Absolute Size Change by Sample", x = "Sample", y = "Size Change (MB)") +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
          print(p)
          return()
        }

        draw_trim_placeholder("Selected trimming visualization is not available.")
      })
    }
  )
}
