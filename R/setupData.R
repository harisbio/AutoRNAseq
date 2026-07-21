# Data setup UI

importDataUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      title = span(icon("file-upload"), " Import & Prepare Data"),
      width = NULL,
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = FALSE,
      style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",

      sidebarLayout(

        sidebarPanel(
          style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",

          h4("🧪 Load Demo Data", style = "color: #0092AC;"),
          actionButton(ns("btn_demo_data"), "Use Demo Data", icon = icon("play"),
                       style = "color: #ffffff; background-color: #F39C12;
                         border-color: #E67E22; padding: 6px 12px; font-size: 14px;
                                   border-radius: 5px; width: auto;"),
          hr(),


          box( title = span(icon("database"), " Upload Custom Data"),
               width = NULL,
               status = "primary",
               solidHeader = TRUE,
               collapsible = TRUE,
               collapsed = TRUE,  # << this collapses it on load
               style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
               #📂 Upload Metadata File
               h4("📋 Upload Metadata", style = "color: #0092AC;"),
               fileInput(ns("input_metadata"), "Upload Study Metadata (.csv)"),

               # 🏷 Select Group Column
               selectInput(ns("select_group_column"), "Choose Group Column:", choices = NULL),

               # 🔍 Data Type Selection
               h4("⚙️ Configuration", style = "color: #0092AC;"),
               selectInput(ns("select_pair_single"), "Choose Type:", choices = c("paired-end", "single-end")),

               # 🔢 Numeric Inputs
               numericInput(ns("numberOfWorkers"), "Parallel Workers (CPU cores):", value = 2, min = 1),
               helpText("Higher values speed processing, but can make your computer slower while running."),

               # 📂 FASTQ source mode
               h4("☁️ Fetch NCBI SRA Data", style = "color: #0092AC;"),
               shiny::textAreaInput(ns("sra_accessions"), 
                            "Enter SRA Accessions (one per line):", 
                            placeholder = "SRR1234567\nSRR1234568",
                            rows = 4),
               helpText("Note: The server will rapidly download the designated SRA files and automatically convert them to .fastq locally."),
               
               # 📄 Display detected files
               verbatimTextOutput(ns("sra_count")),

               br(),
               # 🚀 Load Data Button (Small)
               div(
                 style = "display: flex; justify-content: left;",
                 actionButton(ns("btn_import_data"), "Fetch & Load SRA", icon = icon("cloud-download-alt"),
                              style = "color: #ffffff; background-color: #0092AC;
                                   border-color: #007B9E; padding: 6px 12px; font-size: 14px;
                                   border-radius: 5px; width: auto;")
               )
          )


        ),

        mainPanel(
          tabBox(
            width = 12,

            # 📊 Metadata Table Tab
            tabPanel(
              title = span(icon("table"), " Experimental Design"),
              value = "tab-experiment",
              box(
                title = span(icon("database"), " Experimental Design Data"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("metadata_description"))
              )
            ),

            # 📂 Uploaded FASTQ Files Tab
            tabPanel(
              title = span(icon("file"), " Uploaded FASTQ Data"),
              value = "tab-uploadedData",
              box(
                title = span(icon("file"), " Uploaded Data Details"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("data_description"))
              )
            ),

            # 📌 Data Summary & Insights Panel
            tabPanel(
              title = span(icon("chart-bar"), " Data Summary & Insights"),
              value = "tab-data-summary",
              style = "background-color: #F8F9FA; border-radius: 10px; padding: 15px;",

              # 📂 File Selection
              h4("📂 Select Data File", style = "color: #0092AC;"),
              selectInput(ns("select_file"), "Choose a File:", choices = c("None")),

              # 📊 Select Data Representation
              h4("📊 Choose Visualisation", style = "color: #0092AC;"),
              selectInput(ns("select_representation"), "Select Data Representation:", choices = c(
                "Most Frequent Reads",
                "Read Width Distribution",
                "Unique/Duplicated Reads",
                "Mean Quality Distribution",
                "Cycle-Specific Quality",
                "Cycle-Specific Base Call"
              )),
              box(
                title = span(icon("chart-line"), " Data Overview"),
                width = NULL,
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                style = "background-color: #ffffff; border-radius: 10px; padding: 15px;",
                DT::dataTableOutput(ns("data_overview"))
              )
            )),

          br(),

          # ❓ Help Button (Small)
          div(
            style = "display: flex; justify-content: left;",
            actionButton(ns("btn_tutorial"),
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




#Data setup server

importDataServer <- function(input, output, session) {
    ns <- session$ns
    max_uploaded_files <- 50L



    # Show help modal
    observeEvent(input$btn_tutorial, {
      show_step_help(
        "Data Setup Help",
        c(
          "Metadata must be a CSV file with one row per sample.",
          "Recommended columns: sample ID, condition/group, and optional run/accession columns.",
          "For SRA import, enter one accession per line (SRR, ERR, or DRR).",
          "Demo data is built in if you want to explore the app without uploading files.",
          "Practical size guidance: small pilot studies are easiest to test first; large cohorts may take much longer to download and convert."
        ),
        "Tip: use short, consistent sample names and keep the sample column aligned with the FASTQ file names you expect to see later."
      )
    })

    # Reactive values to store metadata and QC result
    metadata <- reactiveVal(NULL)
    result <- reactiveVal(NULL)

    # Load metadata from file
    observeEvent(input$input_metadata, {
      req(input$input_metadata)
      df <- tryCatch(
        read.csv(input$input_metadata$datapath, header = TRUE, check.names = FALSE),
        error = function(e) {
          showNotification(paste("Metadata file could not be read:", conditionMessage(e)), type = "error")
          NULL
        }
      )
      if (!is.null(df)) {
        metadata(df)
        if (ncol(df) > 0) {
          updateSelectInput(session, "select_group_column", choices = names(df), selected = names(df)[1])
        }
      }
    })

    # Update selectInput for group column once metadata is available
    observe({
      meta <- metadata()
      if (!is.null(meta)) {
        updateSelectInput(session, "select_group_column",
                          choices = names(meta),
                          selected = names(meta)[1])
      }
    })

    # Display metadata
    output$metadata_description <- DT::renderDataTable({
      meta <- metadata()
      if (is.null(meta)) {
        DT::datatable(data.frame(Message = "No metadata loaded yet."))
      } else {
        DT::datatable(meta, options = list(scrollX = TRUE))
      }
    })

    # Demo FASTQ file detection
    demo_fastq_files <- reactive({
      path <- system.file("extdata", "Data", package = "AutoRNAseq")
      list.files(path, pattern = "\\.(fq|fastq)(\\.gz)?$", full.names = TRUE)
    })

    # Get SRA Accessions dynamically from TextArea OR Metadata
    custom_sra_accessions <- reactive({
      txt <- input$sra_accessions
      accs <- character(0)
      if (!is.null(txt) && nzchar(txt)) {
        accs <- unlist(strsplit(txt, "[\r\n, ]+"))
        accs <- accs[nzchar(accs)]
      }
      
      if (length(accs) == 0) {
        meta <- metadata()
        if (!is.null(meta)) {
          for (col in names(meta)) {
            if (any(grepl("^SRR|^ERR|^DRR", meta[[col]], ignore.case=TRUE))) {
              accs <- as.character(meta[[col]])
              accs <- accs[nzchar(accs)]
              break
            }
          }
        }
      }
      return(unique(accs))
    })

    # Show number of SRA Accessions
    output$sra_count <- renderText({
      accs <- custom_sra_accessions()
      count <- length(accs)
      if (count == 0) return("No SRA Accessions detected.")
      paste("SRA Accessions detected:", count)
    })

    # Process SRA fetch and conversion when "Fetch & Load SRA" is clicked
    observeEvent(input$btn_import_data, {
      req(metadata())
      accs <- custom_sra_accessions()
      
      if (length(accs) == 0) {
        showNotification("No SRA accessions provided. Enter SRR IDs or upload Metadata containing a Run column.", type = "error")
        return()
      }
      
      if (length(accs) > max_uploaded_files) {
        showNotification(paste0("Too many SRR files requested. Max allowed is ", max_uploaded_files), type = "error")
        return()
      }
      group <- factor(metadata()[[input$select_group_column]])
      if (length(group) != length(accs)) {
        showNotification("Metadata rows must match the number of detected SRA accessions.", type = "error")
        return()
      }
      workers <- input$numberOfWorkers

      out_dir <- Sys.getenv("FASTQ_MOUNT_PATH", "/data/fastq")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      total <- length(accs)
      showNotification(
        paste0("⏳ Starting SRA fetch for ", total, " accession(s)..."),
        type = "message", duration = 5
      )

      withProgress(message = "Fetching & Converting SRA...", value = 0, {
        for (i in seq_along(accs)) {
          acc <- accs[i]
          incProgress(1 / (total * 2), detail = paste0("[", i, "/", total, "] Processing: ", acc))

          # 1. Skip if FASTQ results already exist
          fastq_files <- list.files(out_dir, pattern = paste0("^", acc, "(_[1-2])?\\.(fq|fastq)(\\.gz)?$"))
          if (length(fastq_files) > 0) {
            showNotification(
              paste0("⏩ ", acc, " already exists. Skipping."),
              type = "message", duration = 3, id = paste0("skip_", acc)
            )
            next
          }

          # 2. ENA Direct Download Attempt
          showNotification(paste0("🌐 [", i, "/", total, "] Checking EBI ENA for ", acc, "..."), 
                           type = "message", duration = 4, id = paste0("dl_", acc))
          
          # Use ENA API to find FASTQ FTP links
          ena_url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", acc, 
                            "&result=read_run&fields=fastq_ftp&format=tsv")
          ena_success <- FALSE
          
          # Try reading ENA API response
          res_ena <- tryCatch({
            # Using read.delim for better TSV handling
            read.delim(url(ena_url), header = TRUE, stringsAsFactors = FALSE)
          }, error = function(e) NULL)
          
          if (!is.null(res_ena) && nrow(res_ena) > 0 && 
              "fastq_ftp" %in% names(res_ena) && 
              !is.na(res_ena$fastq_ftp[1]) && 
              res_ena$fastq_ftp[1] != "") {
             
             links <- unlist(strsplit(as.character(res_ena$fastq_ftp[1]), ";"))
             ena_success <- TRUE
             
             old_timeout <- getOption("timeout")
             options(timeout = max(7200, old_timeout)) # 2 hours for massive files
             
             for (link in links) {
                # Standardize link to http if protocol is missing
                full_link <- if (!grepl("^http|^ftp", link)) paste0("http://", link) else link
                dest_path <- file.path(out_dir, basename(full_link))
                
                showNotification(paste0("⬇️ ENA: Downloading ", basename(full_link)), 
                                 type = "message", duration = 10, id = paste0("dl_", acc))
                
                res_dl <- tryCatch({
                   download.file(full_link, destfile = dest_path, mode = "wb", quiet = TRUE)
                   0 # Success
                }, error = function(e) {
                   message(paste("Download error:", e$message))
                   1 # Failure
                })
                
                if (res_dl != 0) {
                   ena_success <- FALSE
                   if (file.exists(dest_path)) file.remove(dest_path)
                   break
                }
             }
             options(timeout = old_timeout)
          }
          
          # 3. Handle Fallback to NCBI if ENA failed
          if (ena_success) {
             showNotification(paste0("✅ ENA download successful for ", acc), 
                              type = "message", duration = 5, id = paste0("dl_", acc))
          } else {
            showNotification(paste0("⚠️ ", acc, " not found on ENA (or mirror lag). Falling back to NCBI SRA..."), 
                             type = "warning", duration = 6)

            showNotification(
              paste0("📥 [", i, "/", total, "] Prefetching from NCBI: ", acc, "..."),
              type = "message", duration = 8, id = paste0("dl_", acc)
            )

            # Check for existing .sra file
            sra_file <- file.path(out_dir, acc, paste0(acc, ".sra"))
            if (!file.exists(sra_file)) sra_file <- file.path(out_dir, paste0(acc, ".sra"))

            if (!file.exists(sra_file)) {
              # Use prefetch (SRA Toolkit)
              cmd_prefetch <- paste("prefetch", acc, "--output-directory", out_dir)
              exit_prefetch <- system(cmd_prefetch, ignore.stdout = TRUE, ignore.stderr = TRUE)

              if (exit_prefetch != 0) {
                showNotification(
                  paste0("❌ NCBI prefetch failed for ", acc, ". Data might be restricted or ID is invalid."),
                  type = "error", duration = 10
                )
                next
              }
            }

            # Convert SRA to FASTQ
            showNotification(
              paste0("🔄 [", i, "/", total, "] Converting SRA → FASTQ: ", acc, "..."),
              type = "message", duration = 10, id = paste0("conv_", acc)
            )

            incProgress(1 / (total * 2), detail = paste0("[", i, "/", total, "] Converting: ", acc))

            # Locate the downloaded SRA file
            sra_path_lookup <- c(
                file.path(out_dir, acc, paste0(acc, ".sra")),
                file.path(out_dir, paste0(acc, ".sra"))
            )
            final_sra_file <- sra_path_lookup[file.exists(sra_path_lookup)][1]

            if (!is.na(final_sra_file)) {
              # fasterq-dump with threading
              cmd_dump <- paste("fasterq-dump --split-files -e", workers, "-O", out_dir, final_sra_file)
              exit_dump <- system(cmd_dump, ignore.stdout = TRUE, ignore.stderr = TRUE)

              if (exit_dump != 0) {
                showNotification(paste0("❌ Conversion failed for ", acc), type = "error", duration = 10)
              } else {
                showNotification(paste0("✅ [", i, "/", total, "] ", acc, " ready!"), 
                                 type = "message", duration = 5, id = paste0("conv_", acc))
              }
            } else {
              showNotification(paste0("⚠️ .sra file missing for ", acc, " after prefetch."), type = "warning")
            }
          }
          # Final update for this accession
          incProgress(1 / (total * 2))
        }
      })

      files <- list.files(out_dir, pattern = "\\.(fq|fastq)(\\.gz)?$", full.names = TRUE)

      if (length(files) == 0) {
        showNotification("❌ No FASTQ files found after SRA conversion. Please check accession IDs.", type = "error", duration = 15)
        return()
      }

      showNotification(
        paste0("🧬 All SRA files ready! Running QC on ", length(files), " FASTQ file(s)..."),
        type = "message", duration = 8
      )

      workers <- input$numberOfWorkers

      res <- tryCatch(
        Rqc::rqcQA(files, group = group, workers = workers),
        error = function(e) {
          showNotification(paste("Data setup failed:", conditionMessage(e)), type = "error")
          NULL
        }
      )

      if (!is.null(res)) {
        showNotification("🎉 Data setup complete! Proceed to Quality Control.", type = "message", duration = 10)
      }
      result(res)
    })

    # Process demo metadata and FASTQ files when "Load Demo Data" is clicked
    observeEvent(input$btn_demo_data, {
      demo_path <- system.file("extdata", "SraRunTable.csv", package = "AutoRNAseq")
      df <- tryCatch(read.csv(demo_path, header = TRUE, check.names = FALSE), error = function(e) NULL)
      if (is.null(df)) {
        showNotification("Demo metadata could not be loaded.", type = "error")
        return()
      }
      metadata(df)

      files <- demo_fastq_files()
      if (length(files) == 0) {
        showNotification("Demo FASTQ files were not found in package extdata.", type = "error")
        return()
      }
      meta <- metadata()
      req(meta)
      group_col <- if ("Developmental_Stage" %in% names(meta)) "Developmental_Stage" else names(meta)[1]
      group <- factor(meta[[group_col]])
      if (length(group) != length(files)) {
        showNotification("Demo metadata does not match demo FASTQ file count.", type = "error")
        return()
      }
      res <- tryCatch(
        Rqc::rqcQA(files, group = group, workers = 2),
        error = function(e) {
          showNotification(paste("Demo data setup failed:", conditionMessage(e)), type = "error")
          NULL
        }
      )
      result(res)
    })

    # Show QC summary of loaded FASTQ files
    output$data_description <- DT::renderDataTable({
      res <- result()
      if (is.null(res)) {
        DT::datatable(data.frame(Message = "No FASTQ files processed yet."))
      } else {
        DT::datatable(as.data.frame(perFileInformation(res)), options = list(scrollX = TRUE))
      }
    })



    select_files(session, "select_file", result, "Select Files")

    output$data_overview <- DT::renderDataTable({
      res_obj <- result()
      if (is.null(res_obj)) {
        return(DT::datatable(data.frame(Message = " Setup your data for the analysis")))
      }

      # Check if selected_file is not NULL or empty
      if (is.null(input$select_file) || input$select_file == "None") {
        return(DT::datatable(data.frame(Message = " Setup your data for the analysis")))
      }

      selectedPresentation <- input$select_representation
      selectedFile <- input$select_file

      presentation_functions <- list(
        "Most Frequent Reads" = perFileTopReads,
        "Read Width Distribution" = perReadWidth,
        "Unique/Duplicated Reads" = perReadFrequency,
        "Mean Quality Distribution" = perReadQuality,
        "Cycle-Specific Quality" =  perCycleQuality,
        "Cycle-Specific Base Call" = perCycleBasecall
      )

      if (!(selectedPresentation %in% names(presentation_functions))) {
        return(DT::datatable(data.frame(Message = " Setup your data for the analysis")))
      }

      if (!(selectedFile %in% names(res_obj))) {
        return(DT::datatable(data.frame(Message = "Please choose a valid FASTQ file for preview.")))
      }

      viz_data <- tryCatch(
        presentation_functions[[selectedPresentation]](res_obj[selectedFile]),
        error = function(e) {
          data.frame(Message = paste("Unable to generate preview:", conditionMessage(e)))
        }
      )
      DT::datatable(as.data.frame(viz_data), options = list(scrollX = TRUE, pageLength = 5))
    })

    # Return final result for downstream use
    return(result)
  }



