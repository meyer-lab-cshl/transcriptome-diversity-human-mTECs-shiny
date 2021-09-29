## Libraries ####
library(shinydashboard)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
***REMOVED***
library(devtools)
library(remotes)

## Load Browser functions and annotations ####
source("browser.R")

## Import data ####
ht_tpm <- readRDS("data/tpm_heatmap.rds")
ht_lfc <- readRDS("data/fc_heatmap.rds")

gene_df <- data.frame(
  gene=names(genes(txdb)),
  location=as.character(genes(txdb)),
  row.names=NULL, stringsAsFactors=FALSE)

ui <- dashboardPage(
  dashboardHeader(title = "Epitope Diversity"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("TRS, 5Pseq, and RNAseq Data",
               tabName = "genes", icon = icon("th")),
      menuItem("Heatmaps", tabName = "hts", icon = icon("th")),
      menuItem("About", tabName = "About", icon = icon("th")))),
  dashboardBody(
    tabItems(
      tabItem("genes",
              fluidRow(
                titlePanel("TSR, 5Pseq, and RNAseq Data"),
                DT::dataTableOutput("genes"),
                p("This app displays TSR, 5Pseq, and RNAseq data from Carter et al."),
                p("The Genes track diplays the gene strucure for the given region."),
                p("The mTEC hi/lo TSR tracks display transcription start regions."),
                p("The mTEC hi/lo 5Pseq tracks display the alignments of the reads from the 5Pseq experiments."),
                p("The mTEC hi/lo RNAseq tracks display the alignments of the reads from the RNAseq experiments."),
                browser_ui("browser"))
      ),
      tabItem("hts",
              h2("Heatmaps"),
              h3("TPM Zscores"),
              InteractiveComplexHeatmapOutput("ht_tpm"),
              hr(),
              h3("LogFoldChanges"),
              InteractiveComplexHeatmapOutput("ht_lfc"),
              hr()),
      tabItem("About",
              h2("About"),
              p("About text here"),
              hr())
    )
  )
)

server <- function(input,output,session) {
  callModule(browser_server, "browser")

  output$Sample1 <- DT::renderDataTable(
    server = TRUE,
    selection = "single",
    options = list(pageLength=10), {
      gene_df
    })

  observeEvent(input$genes_rows_selected, {
    loc <- gene_df$location[input$genes_rows_selected]
    updateTextInput(session, "browser-location_str", value=loc)
  })

  makeInteractiveComplexHeatmap(input, output, session, ht_tpm, "ht_tpm")
  makeInteractiveComplexHeatmap(input, output, session, ht_lfc, "ht_lfc")
}

shinyApp(ui, server)
