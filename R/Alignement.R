alignmentUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      width = NULL,
      title = span(icon("dna"), " Mapping Reads to Reference Genome"),
      status = "primary",
      solidHeader = TRUE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
      sidebarLayout(
        sidebarPanel(
          style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",
          h4("Upload Reference Genome", style = "color: #0092AC;"),
          fileInput(ns("reference_file"), "Upload a reference file (FASTA format)"),
          h4("Aligner Selection", style = "color: #0092AC;"),
          selectInput(
            ns("aligner_tool"),
            "Choose Alignment Tool",
            choices = c(
              "Rsubread align" = "rsubread",
              "HISAT2" = "hisat2",
              "Bowtie2" = "bowtie2",
              "STAR" = "star"
            ),
            selected = "rsubread"
          ),
          h4("Select Read Files", style = "color: #0092AC;"),
          fileInput(
            ns("readfile1"), 
            "Upload Read Files (FASTQ format)", 
            multiple = TRUE,
            accept = c(".fastq", ".fq", ".gz")
          ),
          checkboxInput(ns("pair_single"), "Paired-End Sequencing", FALSE),
          uiOutput(ns("read2")),
          h4("Output File Format", style = "color: #0092AC;"),
          selectInput(
            ns("outputFormat"),
            "Choose the Output Format:",
            choices = c("BAM", "SAM"),
            selected = "BAM"
          ),
          textInput(
            ns("align_output_dir"),
            "Alignment Output Directory (optional)",
            placeholder = "Leave blank to create alignementOutput inside reads folder"
          ),
          hr(),
          actionButton(
            ns("btn_lunch_alignement"),
            "Align Reads",
            icon = icon("play"),
            style = "color: #ffffff; background-color: #0092AC; border-color: #007B9E; padding: 6px 12px; font-size: 14px; border-radius: 5px; width: auto;"
          ),
          br(),
          div(
            style = "display: flex; justify-content: left;",
            help_modal_button(ns("help_align"), "Help")
          )
        ),
        mainPanel(
          tabBox(
            width = 12,
            tabPanel(
              title = span(icon("file"), " Reference File"),
              box(
                title = span(icon("database"), " Your Reference File"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("ref_file"))
              )
            ),
            tabPanel(
              title = span(icon("chart-line"), " Alignment Results"),
              box(
                title = span(icon("chart-bar"), " Alignment Results"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("alignement_results"))
              )
            ),
            tabPanel(
              title = span(icon("chart-column"), " Alignment Visualizations"),
              selectInput(
                ns("alignment_viz_type"),
                "Visualization Type",
                choices = c(
                  "Alignment Status Summary" = "status_summary",
                  "Output File Size by Sample" = "output_size",
                  "Alignment Rate by Sample" = "rate_by_sample",
                  "Mapped vs Unmapped (Reads/Rate)" = "mapped_unmapped",
                  "Tool-wise Mean Alignment Rate" = "tool_rate"
                ),
                selected = "status_summary"
              ),
              box(
                title = span(icon("chart-bar"), " Alignment Plots"),
                width = NULL,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                plotOutput(ns("alignment_plot"), height = "430px")
              )
            ),
            tabPanel(
              title = span(icon("terminal"), " Alignment Logs"),
              box(
                title = span(icon("info-circle"), " Execution Status"),
                width = NULL,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                textOutput(ns("alignment_status"))
              ),
              box(
                title = span(icon("terminal"), " Command Logs"),
                width = NULL,
                status = "warning",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                verbatimTextOutput(ns("alignment_log"))
              )
            )
          ),
        )
      )
    )
  )
}


alignmentServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- NS(id)
    reactive_align <- reactiveVal(NULL)
    alignment_status <- reactiveVal("No alignment run started.")
    alignment_log <- reactiveVal("No logs available yet.")
    
    readfile1_val <- reactiveVal(NULL)
    readfile2_val <- reactiveVal(NULL)

    observeEvent(input$help_align, {
      show_step_help(
        "Alignment Help",
        c(
          "Upload a reference FASTA file that matches the genome build used for the RNA-seq experiment.",
          "Choose the alignment tool you want to use and make sure it is installed in the container.",
          "Point the read directory to the folder containing FASTQ/FASTQ.GZ files.",
          "If paired-end data is used, enable paired mode so the second read folder can be supplied.",
          "Alignment outputs can be large; test one sample first if your data is substantial."
        ),
        "Keep read names, reference genome, and sample folders consistent. That is the biggest predictor of a successful alignment run."
      )
    })

    observeEvent(input$pair_single, {
      if (isTRUE(input$pair_single)) {
        output$read2 <- renderUI({
          tagList(
            fileInput(
              ns("readfile2"), 
              "Upload Second Read Files (for paired-end, FASTQ format)", 
              multiple = TRUE,
              accept = c(".fastq", ".fq", ".gz")
            )
          )
        })
      } else {
        output$read2 <- renderUI({ NULL })
      }
    })



    normalize_dir <- function(path) {
      gsub("\\\\", "/", trimws(path))
    }

    list_fastq <- function(dir_path) {
      sort(list.files(
        path = dir_path,
        pattern = "\\.(fq|fastq)(\\.gz)?$",
        full.names = TRUE,
        ignore.case = TRUE
      ))
    }

    sample_name_from_fastq <- function(path) {
      x <- basename(path)
      x <- gsub("\\.(fastq|fq)(\\.gz)?$", "", x, ignore.case = TRUE)
      x
    }

    cmd_failed <- function(cmd_output) {
      status <- attr(cmd_output, "status")
      if (!is.null(status) && !is.na(status) && status != 0) {
        return(TRUE)
      }
      any(grepl("ERROR", cmd_output, fixed = TRUE))
    }

    extract_overall_rate <- function(log_lines) {
      txt <- paste(log_lines, collapse = "\n")
      m <- regmatches(txt, regexpr("[0-9]+\\.?[0-9]*% overall alignment rate", txt, ignore.case = TRUE))
      if (length(m) == 0 || is.na(m) || m == "") return(NA_real_)
      as.numeric(sub("%.*$", "", m))
    }

    extract_star_metric <- function(log_final_path, metric_label) {
      if (!file.exists(log_final_path)) return(NA_real_)
      lines <- tryCatch(readLines(log_final_path, warn = FALSE), error = function(e) character(0))
      if (length(lines) == 0) return(NA_real_)
      idx <- grep(paste0("^\\s*", gsub("([\\^\\$\\|\\(\\)\\[\\]\\*\\+\\?\\.\\\\])", "\\\\\\1", metric_label), "\\s*\\|"), lines)
      if (length(idx) == 0) return(NA_real_)
      val <- sub("^.*\\|\\s*", "", lines[idx[1]])
      val <- gsub("%", "", trimws(val))
      suppressWarnings(as.numeric(val))
    }

    normalize_alignment_df <- function(df, tool_name) {
      if (is.null(df) || nrow(df) == 0) return(df)
      out <- df
      if (!("tool" %in% names(out))) out$tool <- tool_name
      if (!("sample" %in% names(out))) {
        s_col <- names(out)[grepl("sample|file|name", names(out), ignore.case = TRUE)][1]
        if (!is.na(s_col) && nzchar(s_col)) out$sample <- as.character(out[[s_col]])
      }
      if (!("status" %in% names(out))) out$status <- "Completed"
      if (!("output_file" %in% names(out))) out$output_file <- NA_character_

      reads_col <- names(out)[grepl("nread|total.*read|input.*read", names(out), ignore.case = TRUE)][1]
      mapped_col <- names(out)[grepl("nmapped|mapped", names(out), ignore.case = TRUE)][1]
      if (!is.na(reads_col) && !is.na(mapped_col) && reads_col != "" && mapped_col != "") {
        suppressWarnings(total_reads <- as.numeric(out[[reads_col]]))
        suppressWarnings(mapped_reads <- as.numeric(out[[mapped_col]]))
        out$total_reads <- total_reads
        out$mapped_reads <- mapped_reads
        out$unmapped_reads <- pmax(total_reads - mapped_reads, 0, na.rm = TRUE)
        out$alignment_rate <- ifelse(total_reads > 0, (mapped_reads / total_reads) * 100, NA_real_)
      }

      out
    }

    maybe_convert_sam_to_bam <- function(sam_path, wanted_format, log_acc) {
      if (!identical(toupper(wanted_format), "BAM")) {
        return(list(path = sam_path, log = log_acc, status = "Completed"))
      }

      if (Sys.which("samtools") == "") {
        return(list(path = sam_path, log = c(log_acc, "samtools not found; keeping SAM output instead of BAM."), status = "Completed with warning"))
      }

      bam_path <- sub("\\.sam$", ".bam", sam_path, ignore.case = TRUE)
      cmd <- c("view", "-bS", shQuote(sam_path), "-o", shQuote(bam_path))
      conv <- tryCatch(system2("samtools", args = cmd, stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
      if (cmd_failed(conv)) {
        return(list(path = sam_path, log = c(log_acc, "samtools conversion failed:", conv), status = "Completed with warning"))
      }
      return(list(path = bam_path, log = c(log_acc, "samtools conversion to BAM completed.", conv), status = "Completed"))
    }

    run_rsubread <- function(reference_file, readfiles1, readfiles2, paired, output_dir, output_format) {
      logs <- c("Building Rsubread index...")
      buildindex(basename = file.path(output_dir, "reference_index"), reference = reference_file)
      logs <- c(logs, "Rsubread index built.")

      input_format <- if (grepl("\\.gz$", readfiles1[1], ignore.case = TRUE)) "gzFASTQ" else "FASTQ"
      output_format_tag <- toupper(output_format)
      out_ext <- paste0(".", tolower(output_format))
      out_files <- file.path(output_dir, paste0(sample_name_from_fastq(readfiles1), out_ext))

      align_args <- list(
        index = file.path(output_dir, "reference_index"),
        readfile1 = readfiles1,
        input_format = input_format,
        output_format = output_format_tag,
        output_file = out_files,
        phredOffset = 33
      )
      if (paired) {
        align_args$readfile2 <- readfiles2
      }

      stat <- do.call(align, align_args)
      stat_df <- as.data.frame(stat)
      stat_df$tool <- "Rsubread align"
      stat_df$sample <- sample_name_from_fastq(readfiles1)
      stat_df$output_file <- out_files
      stat_df <- normalize_alignment_df(stat_df, "Rsubread align")
      list(df = stat_df, logs = logs, status = "Completed")
    }

    run_hisat2 <- function(reference_file, readfiles1, readfiles2, paired, output_dir, output_format) {
      if (Sys.which("hisat2-build") == "" || Sys.which("hisat2") == "") {
        stop("HISAT2 is not installed on this system.")
      }
      logs <- c("Building HISAT2 index...")
      idx_base <- file.path(output_dir, "hisat2_index", "genome")
      dir.create(dirname(idx_base), recursive = TRUE, showWarnings = FALSE)
      bld <- tryCatch(system2("hisat2-build", args = c(shQuote(reference_file), shQuote(idx_base)), stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
      logs <- c(logs, bld)

      rows <- list()
      for (i in seq_along(readfiles1)) {
        sample <- sample_name_from_fastq(readfiles1[i])
        sam_out <- file.path(output_dir, paste0(sample, ".sam"))
        args <- c("-x", shQuote(idx_base), "-S", shQuote(sam_out))
        if (paired) {
          args <- c(args, "-1", shQuote(readfiles1[i]), "-2", shQuote(readfiles2[i]))
        } else {
          args <- c(args, "-U", shQuote(readfiles1[i]))
        }
        run <- tryCatch(system2("hisat2", args = args, stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
        status <- if (cmd_failed(run)) "Failed" else "Completed"
        conv <- maybe_convert_sam_to_bam(sam_out, output_format, run)
        ov_rate <- extract_overall_rate(run)
        rows[[i]] <- data.frame(
          tool = "HISAT2",
          sample = sample,
          read1 = basename(readfiles1[i]),
          read2 = if (paired) basename(readfiles2[i]) else NA_character_,
          status = if (status == "Failed") "Failed" else conv$status,
          alignment_rate = ov_rate,
          output_file = conv$path,
          stringsAsFactors = FALSE
        )
        logs <- c(logs, paste0("HISAT2 -> ", sample), conv$log, "")
      }
      list(df = normalize_alignment_df(do.call(rbind, rows), "HISAT2"), logs = logs, status = "Completed")
    }

    run_bowtie2 <- function(reference_file, readfiles1, readfiles2, paired, output_dir, output_format) {
      if (Sys.which("bowtie2-build") == "" || Sys.which("bowtie2") == "") {
        stop("Bowtie2 is not installed on this system.")
      }
      logs <- c("Building Bowtie2 index...")
      idx_base <- file.path(output_dir, "bowtie2_index", "genome")
      dir.create(dirname(idx_base), recursive = TRUE, showWarnings = FALSE)
      bld <- tryCatch(system2("bowtie2-build", args = c(shQuote(reference_file), shQuote(idx_base)), stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
      logs <- c(logs, bld)

      rows <- list()
      for (i in seq_along(readfiles1)) {
        sample <- sample_name_from_fastq(readfiles1[i])
        sam_out <- file.path(output_dir, paste0(sample, ".sam"))
        args <- c("-x", shQuote(idx_base), "-S", shQuote(sam_out))
        if (paired) {
          args <- c(args, "-1", shQuote(readfiles1[i]), "-2", shQuote(readfiles2[i]))
        } else {
          args <- c(args, "-U", shQuote(readfiles1[i]))
        }
        run <- tryCatch(system2("bowtie2", args = args, stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
        status <- if (cmd_failed(run)) "Failed" else "Completed"
        conv <- maybe_convert_sam_to_bam(sam_out, output_format, run)
        ov_rate <- extract_overall_rate(run)
        rows[[i]] <- data.frame(
          tool = "Bowtie2",
          sample = sample,
          read1 = basename(readfiles1[i]),
          read2 = if (paired) basename(readfiles2[i]) else NA_character_,
          status = if (status == "Failed") "Failed" else conv$status,
          alignment_rate = ov_rate,
          output_file = conv$path,
          stringsAsFactors = FALSE
        )
        logs <- c(logs, paste0("Bowtie2 -> ", sample), conv$log, "")
      }
      list(df = normalize_alignment_df(do.call(rbind, rows), "Bowtie2"), logs = logs, status = "Completed")
    }

    run_star <- function(reference_file, readfiles1, readfiles2, paired, output_dir, output_format) {
      if (Sys.which("STAR") == "") {
        stop("STAR is not installed on this system.")
      }
      logs <- c("Building STAR genome index...")
      genome_dir <- file.path(output_dir, "star_genome")
      dir.create(genome_dir, recursive = TRUE, showWarnings = FALSE)
      bld_args <- c("--runThreadN", "1", "--runMode", "genomeGenerate", "--genomeDir", shQuote(genome_dir), "--genomeFastaFiles", shQuote(reference_file))
      bld <- tryCatch(system2("STAR", args = bld_args, stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
      logs <- c(logs, bld)

      rows <- list()
      for (i in seq_along(readfiles1)) {
        sample <- sample_name_from_fastq(readfiles1[i])
        out_prefix <- file.path(output_dir, paste0(sample, "_"))
        args <- c("--runThreadN", "1", "--genomeDir", shQuote(genome_dir), "--outFileNamePrefix", shQuote(out_prefix))

        r_args <- c(shQuote(readfiles1[i]))
        if (paired) r_args <- c(r_args, shQuote(readfiles2[i]))
        args <- c(args, "--readFilesIn", r_args)
        if (grepl("\\.gz$", readfiles1[i], ignore.case = TRUE)) {
          args <- c(args, "--readFilesCommand", "zcat")
        }

        if (identical(toupper(output_format), "BAM")) {
          args <- c(args, "--outSAMtype", "BAM", "Unsorted")
        } else {
          args <- c(args, "--outSAMtype", "SAM")
        }

        run <- tryCatch(system2("STAR", args = args, stdout = TRUE, stderr = TRUE), error = function(e) paste("ERROR:", conditionMessage(e)))
        status <- if (cmd_failed(run)) "Failed" else "Completed"
        star_log_final <- file.path(output_dir, paste0(sample, "_Log.final.out"))
        uniq_rate <- extract_star_metric(star_log_final, "Uniquely mapped reads %")
        in_reads <- extract_star_metric(star_log_final, "Number of input reads")
        multi_reads <- extract_star_metric(star_log_final, "Number of reads mapped to multiple loci")
        uniq_reads <- extract_star_metric(star_log_final, "Uniquely mapped reads number")
        mapped_reads <- ifelse(is.na(uniq_reads) && is.na(multi_reads), NA_real_, sum(uniq_reads, multi_reads, na.rm = TRUE))
        out_file <- if (identical(toupper(output_format), "BAM")) {
          file.path(output_dir, paste0(sample, "_Aligned.out.bam"))
        } else {
          file.path(output_dir, paste0(sample, "_Aligned.out.sam"))
        }
        rows[[i]] <- data.frame(
          tool = "STAR",
          sample = sample,
          read1 = basename(readfiles1[i]),
          read2 = if (paired) basename(readfiles2[i]) else NA_character_,
          status = status,
          total_reads = in_reads,
          mapped_reads = mapped_reads,
          alignment_rate = uniq_rate,
          output_file = out_file,
          stringsAsFactors = FALSE
        )
        logs <- c(logs, paste0("STAR -> ", sample), run, "")
      }
      list(df = normalize_alignment_df(do.call(rbind, rows), "STAR"), logs = logs, status = "Completed")
    }

    run_alignment <- reactive({
      if (is.null(input$reference_file$datapath) || is.null(input$readfile1)) {
        return(list(df = NULL, status = "Reference FASTA and read files are required.", logs = "Missing inputs."))
      }

      readfolder1 <- file.path(tempdir(), "align_reads1")
      if (dir.exists(readfolder1)) unlink(readfolder1, recursive = TRUE)
      dir.create(readfolder1, recursive = TRUE)
      
      for (i in seq_len(nrow(input$readfile1))) {
        file.copy(input$readfile1$datapath[i], file.path(readfolder1, input$readfile1$name[i]))
      }
      
      readfiles1 <- list_fastq(readfolder1)
      if (length(readfiles1) == 0) {
        return(list(df = NULL, status = "No valid FASTQ files found in uploaded read files.", logs = "Invalid format."))
      }

      paired <- isTRUE(input$pair_single)
      readfiles2 <- character(0)
      if (paired) {
        if (is.null(input$readfile2)) {
          return(list(df = NULL, status = "Paired-end selected but second read files are missing.", logs = ""))
        }
        
        readfolder2 <- file.path(tempdir(), "align_reads2")
        if (dir.exists(readfolder2)) unlink(readfolder2, recursive = TRUE)
        dir.create(readfolder2, recursive = TRUE)
        
        for (i in seq_len(nrow(input$readfile2))) {
          file.copy(input$readfile2$datapath[i], file.path(readfolder2, input$readfile2$name[i]))
        }
        
        readfiles2 <- list_fastq(readfolder2)
        if (length(readfiles2) == 0) {
          return(list(df = NULL, status = "No valid FASTQ files found in second read uploads.", logs = "Invalid format."))
        }
        if (length(readfiles1) != length(readfiles2)) {
          return(list(df = NULL, status = "Paired-end uploads must have equal FASTQ file counts.", logs = paste(length(readfiles1), "vs", length(readfiles2))))
        }
      }

      output_dir <- trimws(input$align_output_dir)
      if (output_dir == "") {
        # Fallback since we don't have a source directory on host anymore
        if (dir.exists("/data/fastq")) {
          output_dir <- "/data/fastq/alignementOutput"
        } else {
          output_dir <- file.path(getwd(), "alignementOutput")
        }
      }
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

      tool <- input$aligner_tool
      reference_file <- input$reference_file$datapath
      out_fmt <- input$outputFormat

      tryCatch({
        if (tool == "rsubread") {
          run_rsubread(reference_file, readfiles1, readfiles2, paired, output_dir, out_fmt)
        } else if (tool == "hisat2") {
          run_hisat2(reference_file, readfiles1, readfiles2, paired, output_dir, out_fmt)
        } else if (tool == "bowtie2") {
          run_bowtie2(reference_file, readfiles1, readfiles2, paired, output_dir, out_fmt)
        } else {
          run_star(reference_file, readfiles1, readfiles2, paired, output_dir, out_fmt)
        }
      }, error = function(e) {
        list(df = NULL, status = paste("Alignment failed:", conditionMessage(e)), logs = paste("ERROR:", conditionMessage(e)))
      })
    })

    refData <- reactive({
      req(input$reference_file)
      finfo <- file.info(input$reference_file$datapath)
      data.frame(
        File = input$reference_file$name,
        Size_MB = round(finfo$size / (1024^2), 2),
        Status = "Uploaded Successfully",
        stringsAsFactors = FALSE
      )
    })

    output$ref_file <- DT::renderDataTable({
      if (is.null(input$reference_file)) {
        data.frame(Message = "Please import the reference file (FASTA format).")
      } else {
        DT::datatable(refData(), options = list(pageLength = 5, scrollX = TRUE), caption = "Reference File Details")
      }
    })

    observeEvent(input$btn_lunch_alignement, {
      withProgress(message = "Running alignment...", value = 0.2, {
        align_result <- run_alignment()
        incProgress(0.7)
        if (is.null(align_result$df)) {
          showNotification(align_result$status, type = "error")
          alignment_status(align_result$status)
        } else {
          showNotification("Alignment completed successfully.", type = "message")
          alignment_status(paste0("Alignment completed using ", input$aligner_tool, "."))
        }
        alignment_log(paste(align_result$logs, collapse = "\n"))
        reactive_align(align_result$df)
        incProgress(0.1)
      })
    })

    draw_alignment_placeholder <- function(msg) {
      plot(
        x = 0, y = 0, type = "n",
        xlim = c(0, 1), ylim = c(0, 1),
        xlab = "", ylab = "", axes = FALSE,
        main = msg
      )
    }

    prep_alignment_plot_df <- reactive({
      df <- reactive_align()
      if (is.null(df) || nrow(df) == 0) return(data.frame())
      out <- as.data.frame(df, stringsAsFactors = FALSE)
      if (!("sample" %in% names(out))) out$sample <- paste0("Sample_", seq_len(nrow(out)))
      if (!("status" %in% names(out))) out$status <- "Completed"
      if (!("tool" %in% names(out))) out$tool <- input$aligner_tool

      if (!("output_file" %in% names(out))) out$output_file <- NA_character_
      out$output_size_mb <- suppressWarnings(file.info(out$output_file)$size / (1024^2))

      if ("alignment_rate" %in% names(out)) {
        suppressWarnings(out$alignment_rate <- as.numeric(out$alignment_rate))
      } else {
        out$alignment_rate <- NA_real_
      }
      if ("total_reads" %in% names(out)) suppressWarnings(out$total_reads <- as.numeric(out$total_reads))
      if ("mapped_reads" %in% names(out)) suppressWarnings(out$mapped_reads <- as.numeric(out$mapped_reads))
      if ("unmapped_reads" %in% names(out)) suppressWarnings(out$unmapped_reads <- as.numeric(out$unmapped_reads))

      if (!("total_reads" %in% names(out))) out$total_reads <- NA_real_
      if (!("mapped_reads" %in% names(out))) out$mapped_reads <- NA_real_
      if (!("unmapped_reads" %in% names(out))) out$unmapped_reads <- NA_real_

      if (all(is.na(out$mapped_reads)) && !all(is.na(out$alignment_rate))) {
        out$mapped_reads <- out$alignment_rate
        out$unmapped_reads <- 100 - out$alignment_rate
      } else if (!all(is.na(out$total_reads)) && !all(is.na(out$mapped_reads)) && all(is.na(out$unmapped_reads))) {
        out$unmapped_reads <- pmax(out$total_reads - out$mapped_reads, 0, na.rm = TRUE)
      }

      if (all(is.na(out$alignment_rate)) && !all(is.na(out$total_reads)) && !all(is.na(out$mapped_reads))) {
        out$alignment_rate <- ifelse(out$total_reads > 0, (out$mapped_reads / out$total_reads) * 100, NA_real_)
      }

      out
    })

    output$alignment_plot <- renderPlot({
      df <- prep_alignment_plot_df()
      viz <- input$alignment_viz_type
      if (is.null(df) || nrow(df) == 0) {
        draw_alignment_placeholder("Run alignment first to view visualizations.")
        return()
      }

      if (identical(viz, "status_summary")) {
        s <- as.data.frame(table(df$status), stringsAsFactors = FALSE)
        names(s) <- c("status", "count")
        p <- ggplot2::ggplot(s, ggplot2::aes(x = status, y = count, fill = status)) +
          ggplot2::geom_col(alpha = 0.92, width = 0.7) +
          ggplot2::labs(title = "Alignment Status Summary", x = "Status", y = "Count") +
          ggplot2::theme_minimal(base_size = 13) +
          ggplot2::theme(legend.position = "none")
        print(p)
        return()
      }

      if (identical(viz, "output_size")) {
        if (all(is.na(df$output_size_mb))) {
          draw_alignment_placeholder("Output files not found yet to calculate file size.")
          return()
        }
        p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = output_size_mb, fill = tool)) +
          ggplot2::geom_col(alpha = 0.92) +
          ggplot2::labs(title = "Output File Size by Sample", x = "Sample", y = "Output Size (MB)") +
          ggplot2::theme_minimal(base_size = 13) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        print(p)
        return()
      }

      if (identical(viz, "rate_by_sample")) {
        if (all(is.na(df$alignment_rate))) {
          draw_alignment_placeholder("Alignment rate is not available for this result.")
          return()
        }
        p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = alignment_rate, fill = alignment_rate)) +
          ggplot2::geom_col(alpha = 0.95) +
          ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey40") +
          ggplot2::labs(title = "Alignment Rate by Sample", x = "Sample", y = "Alignment Rate (%)") +
          ggplot2::theme_minimal(base_size = 13) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        print(p)
        return()
      }

      if (identical(viz, "mapped_unmapped")) {
        if (all(is.na(df$mapped_reads)) || all(is.na(df$unmapped_reads))) {
          draw_alignment_placeholder("Mapped/unmapped values are not available for this result.")
          return()
        }
        long_df <- rbind(
          data.frame(sample = df$sample, class = "Mapped", value = df$mapped_reads, stringsAsFactors = FALSE),
          data.frame(sample = df$sample, class = "Unmapped", value = df$unmapped_reads, stringsAsFactors = FALSE)
        )
        p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = value, fill = class)) +
          ggplot2::geom_col(position = "dodge", alpha = 0.92) +
          ggplot2::labs(title = "Mapped vs Unmapped by Sample", x = "Sample", y = "Reads (or % when only rate is available)") +
          ggplot2::scale_fill_manual(values = c("Mapped" = "#2A9D8F", "Unmapped" = "#E76F51")) +
          ggplot2::theme_minimal(base_size = 13) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        print(p)
        return()
      }

      if (identical(viz, "tool_rate")) {
        if (all(is.na(df$alignment_rate))) {
          draw_alignment_placeholder("Alignment rate is not available for tool-wise comparison.")
          return()
        }
        sum_df <- stats::aggregate(alignment_rate ~ tool, data = df, FUN = function(x) mean(x, na.rm = TRUE))
        p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = tool, y = alignment_rate, fill = tool)) +
          ggplot2::geom_col(alpha = 0.92, show.legend = FALSE) +
          ggplot2::labs(title = "Mean Alignment Rate by Tool", x = "Tool", y = "Mean Alignment Rate (%)") +
          ggplot2::theme_minimal(base_size = 13)
        print(p)
        return()
      }

      draw_alignment_placeholder("Selected visualization is not available.")
    })

    output$alignment_status <- renderText({
      alignment_status()
    })

    output$alignment_log <- renderText({
      alignment_log()
    })

    output$alignement_results <- DT::renderDataTable({
      if (is.null(reactive_align())) {
        data.frame(Message = "Alignment results will appear here once the process is complete.")
      } else {
        DT::datatable(reactive_align(), options = list(pageLength = 15, scrollX = TRUE), caption = "Alignment Summary Results")
      }
    })
  })
}
