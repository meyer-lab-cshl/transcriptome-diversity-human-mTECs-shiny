## Libraries ####
library(shinydashboard)
library(InteractiveComplexHeatmap)

## Load Browser functions and annotations ####
source("browser.R")

## Create dashboard ####
ui <- dashboardPage(
    skin = "green",
    dashboardHeader(title = "Epitope Diversity in human T cell education"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("About", tabName = "About", icon = icon("project-diagram")),
            menuItem("Epitope Browser",
                     tabName = "genes", icon = icon("wpexplorer")),
            menuItem("Expression Heatmaps", tabName = "hts", icon = icon("th"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem("About",
                    h2("About"),
                    p("The induction of central T-cell tolerance in the thymus depends on the  presentation of peripheral self-epitopes by medullary thymic epithelial cells (mTECs), enabled by a process known as promiscuous gene expression (pGE). Epitope diversity generated during pGE has many contributors, including non-canonical transcription initiation, alternative splicing and expression of endogenous retroelements. Here, we mapped the expression of genome-wide epitopes in immature and mature human mTECs using high-throughput 5'Cap and RNA sequencing."),
                    p("This application allows for the exploration of this rich dataset. In the",  em("Epitope Browser"), "RNA sequencing and 5'Cap reads, as well as the transcription start regions derived from the latter can be explored. In the", em("Expression Heatmap"), ", gene expression and their fold changes between mTEC hi and lo can be explored."),
                    p("We provide this comprehensive epitope map of the human thymus as a resource to the community, fto facilitate the identification of epitopes implicated in auto-immune responses against healthy tissue or immune responses against cancer cells."),
                    p("If you use this resource, please cite: Carter JA, Strömich L, Peacey M, Chapin S, Velten L, Steinmetz LM, Brors B, Pinto S, and Meyer HV (2021) Epitope diversity in human medullary thymusepithelial cells"),
                    p("Analysis code can be found at:", a(href="https://github.com/meyer-lab-cshl/thymus-epitope-mapping", "thymus-epitope-mapping")),
                    p("Data is available at EGA:", a(href="EGAXXX", "EGAXXX")),
                    p("For questions regarding the work, please contact: hmeyer [at] cshl.edu"),
                    hr()),
            tabItem("genes",
                    fluidRow(
                        titlePanel("Epitope Browser"),
                        p("Genome-wide transcription start regions and gene expression in two maturation stages of human medullary thymic epithelial cells (mTECs) - immature, MHCII low expressing mTECs (mTEC lo) and mature, MHCII high expressing mTECs (mTEC hi)."),
                        DT::dataTableOutput("genes"),
                        shiny::tags$li("The mTEC hi/lo TSR tracks display transcription start regions."),
                        shiny::tags$li("The mTEC hi/lo 5Pseq tracks display the alignments of the reads from the 5Pseq experiments."),
                        shiny::tags$li("The mTEC hi/lo RNAseq tracks display the alignments of the reads from the RNAseq experiments."),
                        browser_ui("browser"))
            ),
            tabItem("hts",
                    h2("Expression Heatmap"),
                    InteractiveComplexHeatmapOutput("ht"),
                    hr())
        )
    )
)

server <- function(input, output, session) {
    callModule(browser_server, "browser")

    output$Sample1 <- DT::renderDataTable(
        server = TRUE,
        selection = "single",
        options = list(pageLength=10), {
            gene_annotation
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
