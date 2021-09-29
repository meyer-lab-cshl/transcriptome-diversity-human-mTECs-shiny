## Libraries ####
library(shinydashboard)
library(InteractiveComplexHeatmap)

## Load Browser functions and annotations ####
source("browser.R")

## Create dashboard ####
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
                titlePanel("Genome Browser"),
                p("Genome-wide transcription start regions and gene expression in two maturaion stages of human medullary thymic epithelial cells (mTECs) - immature, MHCII low expressing mTECs (mTEC lo) and mature, MHCII high expressing mTECs (mTEC hi)."),
                DT::dataTableOutput("genes"),
                p("The Genes track diplays the gene structure for the given region."),
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
    loc <- gene_annotation$location[input$genes_rows_selected]
    updateTextInput(session, "browser-location_str", value=loc)
  })

    makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    #makeInteractiveComplexHeatmap(input, output, session, ht_tpm, "ht_tpm")
    #makeInteractiveComplexHeatmap(input, output, session, ht_lfc, "ht_lfc")
}

shinyApp(ui, server)
