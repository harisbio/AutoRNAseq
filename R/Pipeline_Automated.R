automatedPipelineUI <- function(id) {
  ns <- NS(id)
  tagList(
    shinydashboard::box(
      width = NULL,
      title = span(icon("robot"), " One-Click Automated RNA-Seq Pipeline"),
      status = "success",
      solidHeader = TRUE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
      
      shiny::fluidRow(
        shiny::column(width = 4,
          shinydashboard::box(
            width = NULL, title = span(icon("database"), " 1. Data Input"), status = "primary", solidHeader = TRUE,
            fileInput(ns("input_metadata"), "Upload Study Metadata (.csv)"),
            selectInput(ns("select_group_column"), "Choose Condition Column:", choices = NULL),
            selectInput(ns("select_sample_column"), "Choose Sample Name Column:", choices = NULL),
            shiny::textAreaInput(ns("sra_accessions"), "Enter SRA Accessions (one per line):", placeholder = "SRR1234567\nSRR1234568", rows = 4),
            selectInput(ns("select_pair_single"), "Library Layout:", choices = c("paired-end", "single-end")),
            numericInput(ns("numberOfWorkers"), "Parallel Workers (CPU cores):", value = 2, min = 1)
          )
        ),
        shiny::column(width = 4,
          shinydashboard::box(
            width = NULL, title = span(icon("cogs"), " 2. Preprocessing & Alignment"), status = "primary", solidHeader = TRUE,
            h5("Quality Control"),
            selectInput(ns("qc_tool"), "QC Tool:", choices = c("None" = "none", "FastQC" = "fastqc")),
            hr(),
            h5("Trimming"),
            selectInput(ns("trimming_tool"), "Trimming Tool:", choices = c("None" = "none", "fastp" = "fastp", "Trimmomatic" = "trimmomatic")),
            textInput(ns("trim_start"), "Bases to Remove (Start):", value = "0"),
            textInput(ns("trim_minlen"), "Minimum Sequence Length:", value = "20"),
            hr(),
            h5("Alignment"),
            fileInput(ns("reference_fasta"), "Upload Reference Genome (.fasta)", accept = c(".fasta", ".fa", ".fna")),
            selectInput(ns("aligner_tool"), "Alignment Tool:", choices = c("HISAT2" = "hisat2", "Rsubread" = "rsubread", "Bowtie2" = "bowtie2", "STAR" = "star"))
          )
        ),
        shiny::column(width = 4,
          shinydashboard::box(
            width = NULL, title = span(icon("chart-bar"), " 3. Quantification & DE"), status = "primary", solidHeader = TRUE,
            h5("Quantification"),
            fileInput(ns("annotation_gtf"), "Upload Annotation (.gtf)", accept = c(".gtf", ".gtf.gz", ".gff3")),
            helpText("Uses featureCounts under the hood."),
            hr(),
            h5("Differential Expression (DESeq2)"),
            textInput(ns("deseq_min_counts"), "Minimum Count Threshold:", value = "10"),
            textInput(ns("deseq_alpha"), "Significance Alpha:", value = "0.05")
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(width = 12,
          shinydashboard::box(
            width = NULL, title = span(icon("play-circle"), " Execute"), status = "warning", solidHeader = TRUE,
            actionButton(ns("btn_run_pipeline"), "Run Automated Pipeline", icon = icon("rocket"), 
                         style = "color: #fff; background-color: #E74C3C; border-color: #C0392B; padding: 10px 20px; font-size: 18px; border-radius: 8px; width: 100%;"),
            hr(),
            verbatimTextOutput(ns("pipeline_logs")),
            br(),
            help_modal_button(ns("help_pipeline"), "Help")
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(width = 12,
          shinydashboard::tabBox(
            width = 12, title = span(icon("poll"), " Final Results (DESeq2)"),
            tabPanel("DE Results Table", DT::dataTableOutput(ns("final_deseq_table"))),
            tabPanel("Volcano Plot", plotOutput(ns("final_volcano_plot"), height="500px"))
          )
        )
      )
    )
  )
}

automatedPipelineServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Increase upload limits locally for this tab just in case
    options(shiny.maxRequestSize = 8192 * 1024^2)
    
    metadata <- reactiveVal(NULL)
    logs <- reactiveVal("Pipeline ready. Configure settings and click 'Run Automated Pipeline'.")
    deseq_res <- reactiveVal(NULL)
    deseq_plot <- reactiveVal(NULL)
    
    append_log <- function(msg) {
      cur <- logs()
      logs(paste(cur, paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg), sep = "\n"))
    }
    
    observeEvent(input$input_metadata, {
      req(input$input_metadata)
      df <- tryCatch(read.csv(input$input_metadata$datapath, header = TRUE, check.names = FALSE), error=function(e) NULL)
      if (!is.null(df) && nrow(df) > 0) {
        metadata(df)
        if (ncol(df) > 0) {
          updateSelectInput(session, "select_group_column", choices = names(df), selected = names(df)[1])
          updateSelectInput(session, "select_sample_column", choices = names(df), selected = names(df)[1])
        }
      } else {
        showNotification("Metadata file could not be read or is empty.", type = "error")
      }
    })

    observeEvent(input$help_pipeline, {
      show_step_help(
        "Automated Pipeline Help",
        c(
          "Upload a study metadata CSV with sample names and condition/group columns.",
          "Provide SRA accessions only if you want the app to download raw reads automatically.",
          "Upload a reference FASTA before alignment and an annotation GTF/GFF before quantification.",
          "Use a smaller test dataset first to verify the end-to-end configuration.",
          "Large runs can take significant disk space and time because the pipeline performs download, QC, trimming, alignment, quantification, and DE."
        ),
        "For the smoothest run, ensure sample names, FASTQ names, and metadata labels are consistent before clicking Run Automated Pipeline."
      )
    })
    
    observeEvent(input$btn_run_pipeline, {
      req(metadata())
      if (is.null(input$select_group_column) || is.null(input$select_sample_column)) {
        showNotification("Please choose both the condition column and sample name column before running the pipeline.", type = "error")
        return()
      }
      if (is.null(input$reference_fasta) || is.null(input$reference_fasta$datapath)) {
        showNotification("Please upload a reference FASTA file before alignment.", type = "error")
        return()
      }
      if (is.null(input$annotation_gtf) || is.null(input$annotation_gtf$datapath)) {
        showNotification("Please upload an annotation GTF/GFF file before quantification.", type = "error")
        return()
      }
      append_log("--- Starting Automated Pipeline ---")
      
      out_dir <- Sys.getenv("FASTQ_MOUNT_PATH", "/data/fastq")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      workers <- input$numberOfWorkers
      is_paired <- input$select_pair_single == "paired-end"
      
      withProgress(message = "Automated Pipeline Running", value = 0, {
        tryCatch({
        
        # ---------------------------------------------------------
        # 1. Fetch SRA Data
        # ---------------------------------------------------------
        incProgress(0.1, detail = "Phase 1: Fetching SRA Data")
        append_log("Phase 1: Fetching SRA Data")
        
        txt <- input$sra_accessions
        accs <- character(0)
        if (!is.null(txt) && nzchar(txt)) {
          accs <- unlist(strsplit(txt, "[\r\n, ]+"))
          accs <- unique(accs[nzchar(accs)])
        }
        
        if (length(accs) == 0) {
          stop("No SRA accessions provided.")
        }
        
        for (i in seq_along(accs)) {
          acc <- accs[i]
          # Simple check if already exists
          existing <- list.files(out_dir, pattern = paste0("^", acc, "(_[1-2])?\\.(fq|fastq)(\\.gz)?$"))
          if (length(existing) > 0) {
             append_log(paste(acc, "already exists. Skipping download."))
             next
          }
          
          append_log(paste("Downloading via prefetch:", acc))
          cmd_pref <- paste("prefetch", acc, "--output-directory", out_dir)
          system(cmd_pref, ignore.stdout = TRUE, ignore.stderr = TRUE)
          
          append_log(paste("Converting to FASTQ:", acc))
          sra_path <- c(file.path(out_dir, acc, paste0(acc, ".sra")), file.path(out_dir, paste0(acc, ".sra")))
          sra_file <- sra_path[file.exists(sra_path)][1]
          
          if (!is.na(sra_file)) {
             cmd_dump <- paste("fasterq-dump --split-files -e", workers, "-O", out_dir, sra_file)
             system(cmd_dump, ignore.stdout = TRUE, ignore.stderr = TRUE)
          } else {
             append_log(paste("ERROR: SRA file not found for", acc))
          }
        }
        
        # Identify FASTQ files
        all_fastqs <- list.files(out_dir, pattern = "\\.(fq|fastq)(\\.gz)?$", full.names = TRUE)
        if (length(all_fastqs) == 0) {
           stop("No FASTQ files generated. Pipeline aborted.")
        }
        
        # Sort into Read 1 and Read 2
        r1_files <- sort(all_fastqs[grepl("_1\\.|\\.fastq|\\.fq", all_fastqs) & !grepl("_2\\.", all_fastqs)])
        r2_files <- if(is_paired) sort(all_fastqs[grepl("_2\\.", all_fastqs)]) else character(0)
        
        append_log(paste("Found", length(r1_files), "samples."))
        
        # ---------------------------------------------------------
        # 2. Quality Control (Optional)
        # ---------------------------------------------------------
        incProgress(0.2, detail = "Phase 2: Quality Control")
        if (input$qc_tool == "fastqc") {
           append_log("Phase 2: Running FastQC")
           qc_dir <- file.path(out_dir, "Auto_QC")
           dir.create(qc_dir, showWarnings = FALSE)
           for (f in c(r1_files, r2_files)) {
              system2("fastqc", args = c(shQuote(f), "--outdir", shQuote(qc_dir)), stdout = FALSE, stderr = FALSE)
           }
           append_log("FastQC completed.")
        } else {
           append_log("Phase 2: Quality Control skipped.")
        }
        
        # ---------------------------------------------------------
        # 3. Trimming (Optional)
        # ---------------------------------------------------------
        incProgress(0.3, detail = "Phase 3: Trimming")
        trim_dir <- file.path(out_dir, "Auto_Trimmed")
        if (input$trimming_tool != "none") {
           append_log(paste("Phase 3: Running Trimming using", input$trimming_tool))
           dir.create(trim_dir, showWarnings = FALSE)
           
           trim_r1 <- character(length(r1_files))
           trim_r2 <- character(length(r2_files))
           
           for (i in seq_along(r1_files)) {
              in1 <- r1_files[i]
              out1 <- file.path(trim_dir, paste0(tools::file_path_sans_ext(basename(in1)), "_trimmed.fastq.gz"))
              trim_r1[i] <- out1
              
              if (input$trimming_tool == "fastp") {
                 args <- c("-i", shQuote(in1), "-o", shQuote(out1))
                 if (is_paired) {
                    in2 <- r2_files[i]
                    out2 <- file.path(trim_dir, paste0(tools::file_path_sans_ext(basename(in2)), "_trimmed.fastq.gz"))
                    trim_r2[i] <- out2
                    args <- c(args, "-I", shQuote(in2), "-O", shQuote(out2))
                 }
                 if (!is.na(as.numeric(input$trim_start))) args <- c(args, "--trim_front1", input$trim_start)
                 if (!is.na(as.numeric(input$trim_minlen))) args <- c(args, "--length_required", input$trim_minlen)
                 
                 system2("fastp", args = args, stdout = FALSE, stderr = FALSE)
              } else if (input$trimming_tool == "trimmomatic") {
                 trim_cmd <- if(is_paired) "TrimmomaticPE" else "TrimmomaticSE"
                 if (is_paired) {
                    in2 <- r2_files[i]
                    out2 <- file.path(trim_dir, paste0(tools::file_path_sans_ext(basename(in2)), "_trimmed.fastq.gz"))
                    trim_r2[i] <- out2
                    u1 <- file.path(trim_dir, paste0(tools::file_path_sans_ext(basename(in1)), "_unpaired.fastq.gz"))
                    u2 <- file.path(trim_dir, paste0(tools::file_path_sans_ext(basename(in2)), "_unpaired.fastq.gz"))
                    args <- c("-threads", workers, shQuote(in1), shQuote(in2), shQuote(out1), shQuote(u1), shQuote(out2), shQuote(u2))
                 } else {
                    args <- c("-threads", workers, shQuote(in1), shQuote(out1))
                 }
                 if (!is.na(as.numeric(input$trim_start))) args <- c(args, paste0("HEADCROP:", input$trim_start))
                 if (!is.na(as.numeric(input$trim_minlen))) args <- c(args, paste0("MINLEN:", input$trim_minlen))
                 system2(trim_cmd, args = args, stdout = FALSE, stderr = FALSE)
              }
           }
           # Update pointers to trimmed files
           r1_files <- trim_r1
           r2_files <- trim_r2
           append_log("Trimming completed.")
        } else {
           append_log("Phase 3: Trimming skipped.")
        }
        
        # ---------------------------------------------------------
        # 4. Alignment
        # ---------------------------------------------------------
        incProgress(0.5, detail = "Phase 4: Alignment")
        append_log("Phase 4: Alignment")
        req(input$reference_fasta)
        ref_fa <- input$reference_fasta$datapath
        align_dir <- file.path(out_dir, "Auto_Alignment")
        dir.create(align_dir, showWarnings = FALSE)
        
        bam_files <- character(length(r1_files))
        
        if (input$aligner_tool == "rsubread") {
           append_log("Building Rsubread index...")
           Rsubread::buildindex(basename = file.path(align_dir, "ref_index"), reference = ref_fa)
           
           for (i in seq_along(r1_files)) {
              out_bam <- file.path(align_dir, paste0(tools::file_path_sans_ext(basename(r1_files[i])), ".bam"))
              bam_files[i] <- out_bam
              append_log(paste("Aligning:", basename(r1_files[i])))
              
              align_args <- list(
                 index = file.path(align_dir, "ref_index"),
                 readfile1 = r1_files[i],
                 output_format = "BAM",
                 output_file = out_bam,
                 nthreads = workers
              )
              if (is_paired) align_args$readfile2 <- r2_files[i]
              do.call(Rsubread::align, align_args)
           }
        } else if (input$aligner_tool == "hisat2") {
           append_log("Building HISAT2 index...")
           idx_base <- file.path(align_dir, "hisat2_index", "ref")
           dir.create(dirname(idx_base), showWarnings=FALSE)
           system2("hisat2-build", args = c(shQuote(ref_fa), shQuote(idx_base)), stdout=FALSE, stderr=FALSE)
           
           for(i in seq_along(r1_files)) {
              out_sam <- file.path(align_dir, paste0(tools::file_path_sans_ext(basename(r1_files[i])), ".sam"))
              out_bam <- file.path(align_dir, paste0(tools::file_path_sans_ext(basename(r1_files[i])), ".bam"))
              bam_files[i] <- out_bam
              
              args <- c("-p", workers, "-x", shQuote(idx_base), "-S", shQuote(out_sam))
              if (is_paired) args <- c(args, "-1", shQuote(r1_files[i]), "-2", shQuote(r2_files[i]))
              else args <- c(args, "-U", shQuote(r1_files[i]))
              
              append_log(paste("Aligning:", basename(r1_files[i])))
              system2("hisat2", args = args, stdout=FALSE, stderr=FALSE)
              
              append_log(paste("Converting SAM to BAM:", basename(r1_files[i])))
              system2("samtools", args = c("view", "-bS", shQuote(out_sam), "-o", shQuote(out_bam)))
              file.remove(out_sam)
           }
        } else {
           stop(paste("Automated run for", input$aligner_tool, "is not fully implemented. Use Rsubread or HISAT2."))
        }
        append_log("Alignment completed.")
        
        # ---------------------------------------------------------
        # 5. Quantification
        # ---------------------------------------------------------
        incProgress(0.7, detail = "Phase 5: Quantification")
        append_log("Phase 5: Quantification using featureCounts")
        anno_gtf <- input$annotation_gtf$datapath
        
        fc_res <- Rsubread::featureCounts(
           bam_files,
           isGTFAnnotationFile = TRUE,
           annot.ext = anno_gtf,
           isPairedEnd = is_paired,
           nthreads = workers
        )
        count_mat <- as.data.frame(fc_res$counts)
        # Clean column names to match SRA IDs
        colnames(count_mat) <- sapply(colnames(count_mat), function(x) {
           # Extract SRR pattern
           m <- regmatches(x, regexpr("SRR[0-9]+", x))
           if(length(m) > 0) m else x
        })
        append_log("Quantification completed.")
        
        # ---------------------------------------------------------
        # 6. Differential Expression
        # ---------------------------------------------------------
        incProgress(0.85, detail = "Phase 6: Differential Expression")
        append_log("Phase 6: Differential Expression (DESeq2)")
        
        meta_df <- metadata()
        cond_col <- input$select_group_column
        samp_col <- input$select_sample_column
        
        meta_sub <- meta_df[, c(samp_col, cond_col)]
        colnames(meta_sub) <- c("sample", "condition")
        meta_sub$condition <- as.factor(meta_sub$condition)
        
        # Align matrix and metadata
        common <- intersect(colnames(count_mat), meta_sub$sample)
        if (length(common) < 2) {
           stop("Sample names in count matrix do not match metadata.")
        }
        
        count_mat <- count_mat[, common, drop=FALSE]
        meta_sub <- meta_sub[match(common, meta_sub$sample), , drop=FALSE]
        rownames(meta_sub) <- meta_sub$sample
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
           countData = as.matrix(count_mat),
           colData = meta_sub,
           design = ~ condition
        )
        
        min_c <- as.numeric(input$deseq_min_counts)
        alpha <- as.numeric(input$deseq_alpha)
        
        keep <- rowSums(DESeq2::counts(dds)) >= min_c
        dds <- dds[keep, ]
        if (ncol(count_mat) < 2) {
          stop("At least two samples are required for DESeq2.")
        }
        dds <- DESeq2::DESeq(dds)
        res <- DESeq2::results(dds, alpha = alpha)
        res_df <- as.data.frame(res)
        res_df <- res_df[complete.cases(res_df), ]
        
        append_log("Differential Expression completed successfully.")
        incProgress(1.0, detail = "Pipeline Complete!")
        
        # Store results for UI
        deseq_res(res_df)
        
        # Generate Volcano Plot
        v_plot <- EnhancedVolcano::EnhancedVolcano(
            res_df,
            lab = rownames(res_df),
            x = "log2FoldChange",
            y = "padj",
            title = "Automated Pipeline: Differential Expression",
            pCutoff = alpha,
            FCcutoff = 1.0,
            pointSize = 2.5,
            labSize = 4.0
        )
        deseq_plot(v_plot)
        }, error = function(e) {
          append_log(paste("ERROR:", conditionMessage(e)))
          showNotification(paste("Automated pipeline failed:", conditionMessage(e)), type = "error", duration = 10)
        })
      })
    })
    
    output$pipeline_logs <- renderText({ logs() })
    output$final_deseq_table <- DT::renderDataTable({
       if (is.null(deseq_res())) return(data.frame(Message="Run pipeline first."))
       DT::datatable(deseq_res(), options = list(scrollX=TRUE, pageLength=10))
    })
    output$final_volcano_plot <- renderPlot({
       if (is.null(deseq_plot())) {
          plot(1, type="n", main="Run pipeline to see Volcano Plot", axes=FALSE, xlab="", ylab="")
          return()
       }
       print(deseq_plot())
    })
    
  })
}
