# ============================================================
# Machine Learning Biomarker Discovery Module
# AutoRNAseq v2
# ============================================================

mlBiomarkerUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shinydashboard::box(
            width = NULL,
            title = shiny::span(shiny::icon("robot"), " Machine Learning Biomarker Discovery"),
            status = "primary",
            solidHeader = TRUE,
            style = "background-color: #EAF2F8; border-radius: 10px; padding: 15px;",
            shinydashboard::tabBox(
                width = 12,

                # ── Tab 1: Data Upload & Config ──────────────────────────────────
                shiny::tabPanel(
                    title = shiny::span(shiny::icon("upload"), "  Data Upload"),
                    value = "ml-upload",
                    shiny::fluidRow(
                        shiny::column(
                            4,
                            shiny::wellPanel(
                                shiny::h4("Data Input", style = "color: #0092AC;"),
                                shiny::helpText("Upload a normalised expression matrix (genes \u00d7 samples, CSV)."),
                                shiny::fileInput(ns("expr_matrix"), "Expression Matrix (.csv)"),
                                shiny::helpText("Upload sample labels: columns 'sample' and 'condition' (0/1 or factor)."),
                                shiny::fileInput(ns("sample_labels"), "Sample Labels (.csv)"),
                                shiny::hr(),
                                shiny::h4("Model Configuration", style = "color: #0092AC;"),
                                shiny::selectInput(ns("ml_method"), "ML Method:",
                                    choices = c(
                                        "Random Forest" = "rf",
                                        "Support Vector Machine" = "svmRadial",
                                        "Logistic Regression" = "glm",
                                        "XGBoost" = "xgbTree"
                                    ),
                                    selected = "rf"
                                ),
                                shiny::numericInput(ns("cv_folds"), "CV Folds:", value = 10, min = 2, max = 20),
                                shiny::numericInput(ns("top_features"), "Top N features to report:", value = 20, min = 5, max = 100),
                                shiny::checkboxInput(ns("use_lasso"), "Pre-select features via LASSO", TRUE),
                                shiny::hr(),
                                shiny::actionButton(ns("btn_run_ml"), "Train Model",
                                    icon  = shiny::icon("play"),
                                    style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px; width:100%;"
                                )
                            )
                        ),
                        shiny::column(
                            8,
                            shinydashboard::box(
                                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                title = "Expression Matrix Preview",
                                DT::dataTableOutput(ns("expr_preview"))
                            ),
                            shinydashboard::box(
                                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                title = "Sample Labels Preview",
                                DT::dataTableOutput(ns("labels_preview"))
                            )
                        )
                    )
                ),

                # ── Tab 2: Feature Importance ────────────────────────────────────
                shiny::tabPanel(
                    title = shiny::span(shiny::icon("star"), "  Feature Importance"),
                    value = "ml-features",
                    shinydashboard::box(
                        width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        title = "Top Feature Table",
                        DT::dataTableOutput(ns("feature_table"))
                    ),
                    shiny::plotOutput(ns("feature_plot"), height = "500px"),
                    shiny::br(),
                    shiny::div(
                        style = "display:flex; gap:10px;",
                        shiny::downloadButton(ns("dl_features"), "Download Table",
                            style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                        ),
                        shiny::downloadButton(ns("dl_feature_plot"), "Download Plot",
                            style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                        )
                    )
                ),

                # ── Tab 3: Model Performance ─────────────────────────────────────
                shiny::tabPanel(
                    title = shiny::span(shiny::icon("chart-line"), "  Model Performance"),
                    value = "ml-perf",
                    shiny::fluidRow(
                        shiny::column(
                            6,
                            shinydashboard::box(
                                width = NULL, status = "info", solidHeader = TRUE,
                                title = "Performance Metrics",
                                DT::dataTableOutput(ns("perf_table"))
                            )
                        ),
                        shiny::column(
                            6,
                            shinydashboard::box(
                                width = NULL, status = "info", solidHeader = TRUE,
                                title = "Confusion Matrix",
                                DT::dataTableOutput(ns("conf_matrix"))
                            )
                        )
                    ),
                    shiny::plotOutput(ns("roc_plot"), height = "450px"),
                    shiny::br(),
                    shiny::downloadButton(ns("dl_roc"), "Download ROC Plot",
                        style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                    )
                ),

                # ── Tab 4: Biomarker Ranking ─────────────────────────────────────
                shiny::tabPanel(
                    title = shiny::span(shiny::icon("trophy"), "  Biomarker Ranking"),
                    value = "ml-ranking",
                    shiny::helpText("Biomarker score = normalised ML importance + normalised DEG significance. Upload DEG CSV to enable full scoring."),
                    shiny::fileInput(ns("deg_for_ranking"), "Upload DEG CSV (optional)"),
                    shinydashboard::box(
                        width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        title = "Ranked Biomarker Candidates",
                        DT::dataTableOutput(ns("biomarker_rank"))
                    ),
                    shiny::br(),
                    shiny::downloadButton(ns("dl_biomarkers"), "Download Table",
                        style = "color:#fff; background-color:#0092AC; border-color:#007B9E; padding:6px 12px; font-size:14px; border-radius:5px;"
                    )
                )
            ),
            shiny::br(),
            shiny::actionButton(ns("help_ml"), "Help",
                icon  = shiny::icon("info-circle"),
                style = "color:#fff; background-color:#F39C12; border-color:#E67E22; padding:6px 12px; font-size:14px; border-radius:5px;"
            )
        )
    )
}


mlBiomarkerServer <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        shiny::observeEvent(input$help_ml, {
            show_step_help(
              "ML Biomarker Help",
              c(
                "Expression matrix must be a CSV with genes as rows and samples as columns.",
                "Sample labels must have sample and condition columns, with sample names matching the expression matrix columns.",
                "You need at least two classes and enough samples per class for cross-validation.",
                "Random Forest is a good first choice if you are not sure which model to start with.",
                "Optional DEG files used for ranking should include gene names, log2FoldChange, and padj."
              ),
              "If your dataset is very large, start with a smaller subset to confirm the file structure before training a full model."
            )
        })

        # ── Load data ────────────────────────────────────────────────────────────
        expr_data <- shiny::reactive({
            shiny::req(input$expr_matrix)
            utils::read.csv(input$expr_matrix$datapath, header = TRUE, check.names = FALSE, row.names = 1)
        })

        label_data <- shiny::reactive({
            shiny::req(input$sample_labels)
            utils::read.csv(input$sample_labels$datapath, header = TRUE, check.names = FALSE)
        })

        validate_ml_inputs <- shiny::reactive({
            expr <- expr_data()
            labels <- label_data()

            if (!all(c("sample", "condition") %in% names(labels))) {
                stop("Sample labels must contain the columns 'sample' and 'condition'.")
            }
            if (anyDuplicated(colnames(expr)) > 0) {
                stop("Expression matrix contains duplicate sample names. Please make columns unique.")
            }
            if (anyDuplicated(labels$sample) > 0) {
                stop("Sample labels contain duplicate sample names. Please keep one row per sample.")
            }
            common <- intersect(colnames(expr), labels$sample)
            if (length(common) < 2) {
                stop("At least two matching samples are required between the expression matrix and sample labels.")
            }
            labels_sub <- labels[match(common, labels$sample), , drop = FALSE]
            if (anyNA(labels_sub$condition)) {
                stop("Some matching samples are missing a condition value.")
            }
            if (length(unique(labels_sub$condition)) < 2) {
                stop("Sample labels must contain at least two condition groups for classification.")
            }
            counts <- table(labels_sub$condition)
            min_class <- min(counts)
            if (min_class < 2) {
                stop("Each condition group must contain at least two samples for cross-validation.")
            }
            list(expr = expr, labels = labels, common = common, labels_sub = labels_sub, min_class = min_class)
        })

        output$expr_preview <- DT::renderDataTable({
            shiny::req(expr_data())
            DT::datatable(head(expr_data(), 100), options = list(scrollX = TRUE, pageLength = 5))
        })

        output$labels_preview <- DT::renderDataTable({
            shiny::req(label_data())
            DT::datatable(label_data(), options = list(scrollX = TRUE, pageLength = 5))
        })

        # ── LASSO feature selection ───────────────────────────────────────────────
        lasso_features <- shiny::reactive({
            shiny::req(validate_ml_inputs())
            if (!isTRUE(input$use_lasso)) {
                return(rownames(validate_ml_inputs()$expr))
            }

            vi <- validate_ml_inputs()
            labels <- vi$labels
            expr <- vi$expr
            common <- vi$common
            expr <- expr[, common, drop = FALSE]
            y <- factor(labels$condition[match(common, labels$sample)])
            x <- t(as.matrix(expr))

            tryCatch(
                {
                    cv_fit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)
                    coef_m <- stats::coef(cv_fit, s = "lambda.1se")
                    selected <- rownames(coef_m)[coef_m[, 1] != 0]
                    selected <- selected[selected != "(Intercept)"]
                    if (length(selected) < 5) {
                        return(rownames(expr))
                    }
                    selected
                },
                error = function(e) rownames(expr)
            )
        })

        # ── Train caret model ─────────────────────────────────────────────────────
        ml_results <- shiny::reactiveVal(NULL)

        shiny::observeEvent(input$btn_run_ml, {
            shiny::req(validate_ml_inputs())
            shiny::withProgress(message = "Training model...", value = 0.2, {
                tryCatch(
                    {
                        vi <- validate_ml_inputs()
                        labels <- vi$labels
                        expr <- vi$expr
                        common <- vi$common
                        expr <- expr[, common, drop = FALSE]
                        y <- factor(labels$condition[match(common, labels$sample)])

                        feats <- lasso_features()
                        feats <- feats[feats %in% rownames(expr)]
                        if (length(feats) < 2) {
                            stop("Not enough usable features after filtering. Please provide a richer expression matrix.")
                        }
                        x_df <- as.data.frame(t(expr[feats, , drop = FALSE]))

                        shiny::setProgress(0.4, message = "Cross-validating...")
                        folds <- min(as.integer(input$cv_folds), vi$min_class)
                        folds <- max(2L, folds)
                        ctrl <- caret::trainControl(
                            method          = "cv",
                            number          = folds,
                            classProbs      = TRUE,
                            summaryFunction = caret::twoClassSummary,
                            savePredictions = "final"
                        )
                        y <- factor(make.names(y))
                        model <- caret::train(
                            x = x_df, y = y, method = input$ml_method,
                            metric = "ROC", trControl = ctrl
                        )
                        shiny::setProgress(0.9, message = "Computing metrics...")
                        ml_results(model)
                        shiny::showNotification("Model trained successfully!", type = "message")
                    },
                    error = function(e) {
                        shiny::showNotification(paste("ML Error:", conditionMessage(e)), type = "error")
                    }
                )
            })
        })

        # ── Feature importance ────────────────────────────────────────────────────
        feat_df <- shiny::reactive({
            shiny::req(ml_results())
            imp <- caret::varImp(ml_results())$importance
            if (is.null(imp) || nrow(imp) == 0) {
                return(data.frame(Gene = character(0), Importance = numeric(0)))
            }
            imp$Gene <- rownames(imp)
            imp$Importance <- round(imp[, 1], 4)
            imp <- imp[order(-imp$Importance), c("Gene", "Importance")]
            head(imp, input$top_features)
        })

        output$feature_table <- DT::renderDataTable({
            shiny::req(feat_df())
            DT::datatable(feat_df(), options = list(scrollX = TRUE, pageLength = 15))
        })

        output$feature_plot <- shiny::renderPlot({
            shiny::req(feat_df())
            df <- feat_df()
            ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Gene, Importance), y = Importance, fill = Importance)) +
                ggplot2::geom_col() +
                ggplot2::coord_flip() +
                ggplot2::scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
                ggplot2::labs(
                    title = paste("Top", nrow(df), "Features —", input$ml_method),
                    x = "Gene", y = "Importance"
                ) +
                ggplot2::theme_minimal(base_size = 13) +
                ggplot2::theme(legend.position = "none")
        })

        # ── Performance ────────────────────────────────────────────────────────────
        output$perf_table <- DT::renderDataTable({
            shiny::req(ml_results())
            preds <- ml_results()$pred
            cm <- caret::confusionMatrix(preds$pred, preds$obs)
            byClass <- cm$byClass
            metrics <- data.frame(
                Metric = c("Accuracy", "Kappa", "Sensitivity", "Specificity", "Precision", "Recall", "F1"),
                Value = round(c(
                    cm$overall["Accuracy"], cm$overall["Kappa"],
                    byClass["Sensitivity"], byClass["Specificity"],
                    byClass["Precision"], byClass["Recall"], byClass["F1"]
                ), 4)
            )
            DT::datatable(metrics, options = list(dom = "t", paging = FALSE))
        })

        output$conf_matrix <- DT::renderDataTable({
            shiny::req(ml_results())
            preds <- ml_results()$pred
            cm <- caret::confusionMatrix(preds$pred, preds$obs)
            DT::datatable(as.data.frame.matrix(cm$table), options = list(dom = "t", paging = FALSE))
        })

        output$roc_plot <- shiny::renderPlot({
            shiny::req(ml_results())
            preds <- ml_results()$pred
            lvls <- levels(preds$obs)
            prob_col <- lvls[2]
            roc_obj <- pROC::roc(preds$obs, preds[[prob_col]], levels = rev(lvls))
            auc_val <- round(pROC::auc(roc_obj), 4)
            graphics::plot(roc_obj,
                col = "#1A5276", lwd = 2,
                main = paste0("ROC Curve — AUC = ", auc_val),
                xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)"
            )
            graphics::abline(a = 0, b = 1, col = "gray", lty = 2)
        })

        # ── Biomarker ranking ──────────────────────────────────────────────────────
        biomarker_table <- shiny::reactive({
            shiny::req(feat_df())
            df <- feat_df()
            colnames(df) <- c("Gene", "ML_Importance")
            df$ML_Norm <- (df$ML_Importance - min(df$ML_Importance)) /
                (max(df$ML_Importance) - min(df$ML_Importance) + 1e-9)

            if (!is.null(input$deg_for_ranking) && !is.null(input$deg_for_ranking$datapath)) {
                deg <- utils::read.csv(input$deg_for_ranking$datapath, header = TRUE, check.names = FALSE)
                if (ncol(deg) < 3 || !all(c("log2FoldChange", "padj") %in% names(deg))) {
                    stop("DEG file must include gene names plus 'log2FoldChange' and 'padj' columns.")
                }
                colnames(deg)[1] <- "Gene"
                deg <- deg[, c("Gene", "log2FoldChange", "padj")]
                df <- merge(df, deg, by = "Gene", all.x = TRUE)
                max_padj <- suppressWarnings(max(df$padj, na.rm = TRUE))
                if (!is.finite(max_padj) || max_padj <= 0) {
                    df$DEG_Norm <- 0
                } else {
                    df$DEG_Norm <- 1 - (df$padj / max_padj)
                }
                df$DEG_Norm[is.na(df$DEG_Norm)] <- 0
                df$Biomarker_Score <- round((df$ML_Norm + df$DEG_Norm) / 2, 4)
                df <- df[, c("Gene", "ML_Importance", "log2FoldChange", "padj", "Biomarker_Score")]
            } else {
                df$Biomarker_Score <- round(df$ML_Norm, 4)
                df <- df[, c("Gene", "ML_Importance", "Biomarker_Score")]
            }
            df[order(-df$Biomarker_Score), ]
        })

        output$biomarker_rank <- DT::renderDataTable({
            shiny::req(biomarker_table())
            DT::datatable(biomarker_table(), options = list(scrollX = TRUE, pageLength = 15))
        })

        # ── Downloads ─────────────────────────────────────────────────────────────
        output$dl_features <- shiny::downloadHandler(
            filename = function() paste0("ML_features_", utils::Sys.Date(), ".csv"),
            content  = function(file) utils::write.csv(feat_df(), file, row.names = FALSE)
        )
        output$dl_biomarkers <- shiny::downloadHandler(
            filename = function() paste0("Biomarker_ranking_", utils::Sys.Date(), ".csv"),
            content  = function(file) utils::write.csv(biomarker_table(), file, row.names = FALSE)
        )
        output$dl_feature_plot <- shiny::downloadHandler(
            filename = function() paste0("Feature_importance_", utils::Sys.Date(), ".png"),
            content = function(file) {
                df <- feat_df()
                p <- ggplot2::ggplot(df, ggplot2::aes(
                    x = stats::reorder(Gene, Importance),
                    y = Importance, fill = Importance
                )) +
                    ggplot2::geom_col() +
                    ggplot2::coord_flip() +
                    ggplot2::scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
                    ggplot2::labs(title = "Feature Importance", x = "Gene", y = "Importance") +
                    ggplot2::theme_minimal(base_size = 13)
                grDevices::png(file, width = 1200, height = 800, res = 120)
                print(p)
                grDevices::dev.off()
            }
        )
        output$dl_roc <- shiny::downloadHandler(
            filename = function() paste0("ROC_curve_", utils::Sys.Date(), ".png"),
            content = function(file) {
                preds <- ml_results()$pred
                lvls <- levels(preds$obs)
                prob_col <- lvls[2]
                roc_obj <- pROC::roc(preds$obs, preds[[prob_col]], levels = rev(lvls))
                auc_val <- round(pROC::auc(roc_obj), 4)
                grDevices::png(file, width = 1000, height = 800, res = 120)
                graphics::plot(roc_obj,
                    col = "#1A5276", lwd = 2,
                    main = paste0("ROC Curve — AUC = ", auc_val)
                )
                graphics::abline(a = 0, b = 1, col = "gray", lty = 2)
                grDevices::dev.off()
            }
        )

        # ── Server call for mlBiomarkerServer ─────────────────────────────────────
    })
}

