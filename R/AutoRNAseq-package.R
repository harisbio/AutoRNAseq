#' @keywords internal
"_PACKAGE"

## usethis namespace: start

# Data manipulation
#' @importFrom dplyr filter mutate select rename
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom utils read.csv write.csv sessionInfo

# Data tables and rendering
#' @importFrom DT datatable renderDataTable dataTableOutput

# Shiny framework (split by purpose for clarity)
# Basic shiny
#' @importFrom shiny NS icon fluidPage tags tagList
# Inputs
#' @importFrom shiny fileInput selectInput textInput numericInput radioButtons checkboxInput actionButton downloadButton
# Outputs
#' @importFrom shiny uiOutput helpText textOutput plotOutput verbatimTextOutput renderText renderPrint renderUI
# Layout
#' @importFrom shiny sidebarLayout sidebarPanel mainPanel tabPanel hr h4 h5 br span div fluidRow column wellPanel
# Logic
#' @importFrom shiny reactive observeEvent observe reactiveVal req showModal showNotification callModule moduleServer modalDialog modalButton withProgress setProgress incProgress runApp HTML updateSelectInput includeMarkdown img downloadHandler

#' @importFrom shinydashboard menuItem menuSubItem sidebarMenu tabItems tabItem box tabBox dashboardBody dashboardSidebar dashboardHeader dashboardPage

# Widgets (shinyWidgets, shinycssloaders, alerts)
#' @importFrom shinyWidgets textInputIcon
#' @import shinyalert
#' @import shinycssloaders

# Tutorials and help
#' @import rintrojs
#' @importFrom rintrojs introjsUI

# RNA-seq related packages
#' @import Rqc
#' @import QuasR
#' @import Rsubread
#' @import DESeq2
#' @import apeglm

# Plotting
#' @import ggplot2
#' @import ggplotify
#' @import EnhancedVolcano
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom AnnotationDbi mapIds
#' @import visNetwork

# Other
#' @importFrom tools file_path_sans_ext

## usethis namespace: end

NULL
