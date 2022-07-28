## Libraries ####
library(shiny)
library(shinydashboard)
library(InteractiveComplexHeatmap)

## Load Browser functions and annotations ####
source("browser.R")

## Create dashboard ####
ui <- dashboardPage(
    skin = "green",
    dashboardHeader(title = "Transcriptomic Diversity in human T cell education",
                    titleWidth = 450),
    dashboardSidebar(
        # width = 450,
        sidebarMenu(
            menuItem(
                "About",
                tabName = "About",
                icon = icon("project-diagram")
            ),
            menuItem(
                "Transcriptome Browser",
                tabName = "genes",
                icon = icon("wpexplorer")
            ),
            menuItem(
                "Expression Heatmaps",
                tabName = "hts",
                icon = icon("th")
            ),
            menuItem(
              "Methods",
              tabName = "mtds",
              icon = icon("th")
            )
        )),
    dashboardBody(tabItems(
        tabItem("About",
                div(
                    h2("About"),
                    p(
                        "The induction of central T-cell tolerance in the thymus depends on the  presentation of peripheral self-epitopes by medullary thymic epithelial cells (mTECs), enabled by a process known as promiscuous gene expression (pGE). Transcriptomic diversity generated during pGE has many contributors, including non-canonical transcription initiation, alternative splicing and expression of endogenous retroelements. Here, we mapped the expression of genome-wide epitopes in immature and mature human mTECs (MHCII expression is low and high, respectively; refered to as mTEC lo and mTEC hi)  using high-throughput 5'Cap and RNA sequencing."
                    ),
                    p(
                        "This application allows for the exploration of this rich dataset. In the",
                        em("Transcriptome Browser"),
                        "RNA sequencing and 5'Cap reads, as well as the transcription start regions derived from the latter can be explored. In the",
                        em("Expression Heatmap"),
                        ", gene expression and their fold changes between mTEC hi and lo can be explored."
                    ),
                    p(
                        "We provide this comprehensive transcriptome map of the human thymus as a resource to the community, to facilitate the identification of epitopes implicated in auto-immune responses against healthy tissue or immune responses against cancer cells."
                    ),
                    p(
                        "If you use this resource, please cite: Carter JA, Strömich L, Peacey M, Chapin S, Velten L, Steinmetz LM, Brors B, Pinto S, and Meyer HV (2021)", em("Transcriptomic diversity in human medullary thymus epithelial cells.")
                    ),
                    p(
                        "Analysis code can be found at:",
                        a(href = "https://github.com/meyer-lab-cshl/thymus-epitope-mapping", "transcriptomic-diversity-human-mTECs.")
                    ),
                    p("Data is available at:", a(href = "https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE201720", "GEO Accession.")),
                    p(
                        "For questions regarding the work, please contact: hmeyer [at] cshl.edu"
                    ),
                    hr()
                )),
        tabItem(
            "genes",
            
            div(
                titlePanel("Transcriptome Browser"),
                p(
                    "Genome-wide transcription start regions and gene expression in two maturation stages of human medullary thymic epithelial cells (mTECs) - immature, MHCII low expressing mTECs (mTEC lo) and mature, MHCII high expressing mTECs (mTEC hi)."
                ),
                DT::dataTableOutput("genes"),
                shiny::tags$li("The mTEC hi/lo TSR tracks display transcription start regions."),
                shiny::tags$li(
                    "The mTEC hi/lo 5Pseq tracks display the alignments of the reads from the 5Pseq experiments."
                ),
                shiny::tags$li(
                    "The mTEC hi/lo RNAseq tracks display the alignments of the reads from the RNAseq experiments."
                ),
                p(
                    ""
                ),
                p(
                    "Please note that the transcriptome browser may intially take a few minutes to load below. You may also experience short delays when searching or navigating through the browser, during which the browser will be grayed out."
                ),
                shiny::tags$script(HTML(
                    "var socket_timeout_interval;
var n = 0;

$(document).on('shiny:connected', function(event) {
  socket_timeout_interval = setInterval(function() {
    Shiny.onInputChange('alive_count', n++)
  }, 10000);
});

$(document).on('shiny:disconnected', function(event) {
  clearInterval(socket_timeout_interval)
});"
                ))
            ),
            browser_ui("browser")
        ),
        tabItem(
            "hts",
            div(
                titlePanel("Expression Heatmap"),
                p(
                    "Gene expression (in transcripts per million) and log fold changes in expression of mTEC hi versus mTEC lo samples."
                ),
                p(
                    "Gene expression is displayed and searchable on transcript level (multiple Ensembl Transcript ID's corresponding to the same gene). In order to search for all instances of a particular gene in these cases, select the", em("Regular expression"), "option in the search tab, and then enter the gene name as the ", em("Keyword"),". All letters in the gene names and transcript ID's must be capitalized."
                )
            ),
            InteractiveComplexHeatmapOutput("ht"),
            hr()
        ),
        tabItem(
          "mtds",
          div(
            titlePanel("Methods"),
            p(
              "A detailed description of the methods can be found at: [link to publication]."
            ),
            p(
              strong("mTEC isolation:"), "Thymi were digested with collagenase/dispase followed by trypsin digestion. MTECs were enriched by magnetic cell sorting followed by cell staining and FACS. Magnetic cell depletion of thymocytes was performed using anti-CD45 Microbeads. The enriched stromal cell fraction (CD45-) was stained with biotinylated anti-epithelial cell adhesion molecule, anti-cortical dendritic reticulum antigen 2, anti-HLA-DR and anti-CD45. MTECs were sorted as CD45-, CDR2-, EpCAM+ cells and MHCII (HLA-DR) was used to separate immature mTEClo and mature mTEChi cell populations. Dead cells were excluded with propidium iodide."
            ),
            p(
              strong("Transcription start site analyses:"), "Transcription start sites were sequenced following the detailed protocol described in [Pelechano et al, 2014]. This method pre-selects only those RNA molecules from the RNA lysate that carry a 5'Cap, i.e. captures true start sites of transcripts, ignoring the 5' ends of degraded and fragmented RNA strands. Library preparation included the ligation of unique molecular identifiers to ensure PCR amplifications can be identified. Transcription start sites (TSS) were defined as the 5' position of the uniquely mapped, forward reads. The expression levels of the TSSs for each mTEC sample were normalized to tags per million using a power law normalization implemented in CAGEr (v1.32). Normalized TSSs in each sample were combined into transcription start regions (TSRs) using paraclu (v9) with a minimum tag cluster expression of 2 tags per million and a maximum cluster length of 20bp. TSRs across samples were combined with bedops (v2.4.38) and consensus, strand-specific TSRs in human mTEC samples called by merging TSRs derived from different samples within 20bp proximity using bedtools (v2.29.2). mTEChi-specific TSRs were called by finding the TSRs that were expressed in at least two mTEChi and not detected in any mTEClo samples; mTEClo-specific TSRs were called accordingly. TSRs were annotated using HOMER (v4.11.1) including mapping to the closest gene, calculation of TSR CpG/GC content and TATA motif search with the provided motif file."
            ),
            p(
              strong("Gene expression analyses:"), "Reads were filtered using fastp (v0.11.8) with for minimum phred quality > 25 of at most 10% unqualified bases, a minimum length of 50bp, at least 30% complexity and polyX tail trimming. Filtered reads were pseudo-aligned and quantified using Kallisto (v0.4.6) with 100 bootstrap samples. Sleuth (v0.30.0) was used for differential expression analysis between mTEChi and mTEClo samples, adjusting for patient effect in the reduced and full model. Transcripts with q-value < 0.05 were considered differentially expressed. Gene annotation was added to transcript IDs using biomaRt (v2.46.3)."
            ),
            p(
              strong("Endogenous retroelement (ERE) analysis:"), "Filtered reads (as in Gene expression analysis) were aligned to GRCh38 using STAR (v2.7.2b) with the following options:  sjdbOverhang 100, winAnchorMultimapNmax 200, outFilterMultimapNmax 100. To obtain transcript per million (TPM) counts at a subfamily level, SalmonTE (v0.4) was run on the filtered reads in 'quant' mode with the following options: reference=hs,  exprtype=TPM. ERE counts were extracted (RepeatMasker class 'LTR', 'LINE', 'SINE' or 'Retroposon') and filtered for entries with a minimum of 2 reads across samples. Differential expression analysis on this set of counts was performed with DESeq2 (v1.30.1). EREs with a Benjamini Hochberg -adjusted p-value < 0.05 were considered differentially expressed."
            ),
            p(
              "Pelechano, V., Wei, W., Jakob, P. & Steinmetz, L. Genome-wide identification of transcript start and end sites by transcript isoform sequencing. Nature Protocols 9, 1740–1759 (2014)."
            ),
          ),
        )
    ))
)

server <- function(input, output, session) {
    callModule(browser_server, "browser")
    
    output$Sample1 <- DT::renderDataTable(
        server = TRUE,
        selection = "single",
        options = list(pageLength = 10),
        {
            gene_annotation
        }
    )
    
    observeEvent(input$genes_rows_selected, {
        loc <- gene_annotation$location[input$genes_rows_selected]
        updateTextInput(session, "browser-location_str", value = loc)
    })
    
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    #makeInteractiveComplexHeatmap(input, output, session, ht_tpm, "ht_tpm")
    #makeInteractiveComplexHeatmap(input, output, session, ht_lfc, "ht_lfc")
    
    output$keep_alive <- renderText({
        req(input$alive_count)
        input$alive_count
    })
}

shinyApp(ui, server)
