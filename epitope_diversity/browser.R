## Libraries ####
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(shinydashboard)
library(InteractiveComplexHeatmap)

## Import data objects ####
options(ucscChromosomeNames=FALSE)
txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                            release=NA,
                            circ_seqs=NULL,
                            server="ensembldb.ensembl.org",
                            username="anonymous", password=NULL, port=0L,
                            tx_attrib=NULL)
seqlevelsStyle(txdb) <- "UCSC"
#txdb <- loadDb("data/txdb.db")

bgFile_hi <- rtracklayer::import.bedGraph("data/tsr_tpm_hi_sum.bedGraph")
bgFile_lo <- rtracklayer::import.bedGraph("data/tsr_tpm_lo_sum.bedGraph")

ht_tpm <- readRDS("data/tpm_heatmap.rds")
ht_lfc <- readRDS("data/fc_heatmap.rds")

ht <- ht_tpm + ht_lfc

## Create annotations ####
options(ucscChromosomeNames=FALSE)
genome <- BSgenome.Hsapiens.NCBI.GRCh38
seqlevelsStyle(genome) <- "UCSC"

gene_annotation <- data.frame(
    gene=names(genes(txdb)),
    location=as.character(genes(txdb)),
    row.names=NULL, stringsAsFactors=FALSE)

## Browser Tracks: Setting up ####
## Annotations
axis_track <- GenomeAxisTrack()
seq_track <- SequenceTrack(genome)
gene_track <- GeneRegionTrack(
 txdb, genome = "hg38", name="Genes", showId=TRUE, shape="arrow",
 transcriptAnnotation = "symbol")

## RNAseq data
alTrack_hi <- AlignmentsTrack(
    "data/mtec_hi_Aligned.sortedByCoord.out.bam",
    isPaired = TRUE,
    name = "mTEC hi RNAseq",
    background.title = "#4c72b0ff",
    fill="#4c72b0ff",
    type = c("coverage", "sashimi"))

alTrack_lo <- AlignmentsTrack(
    "data/mtec_lo_Aligned.sortedByCoord.out.bam",
    isPaired = TRUE,
    name = "mTEC lo RNAseq",
    background.title = "#dd8452ff",
    fill="#dd8452ff",
    type = c("coverage", "sashimi"))

## 5'Cap data
alTrack_5cap_hi <- AlignmentsTrack(
    "data/5cap_mtec_hi_Aligned.sortedByCoord.out.bam",
    isPaired = TRUE,
    name = "mTEC hi 5Pseq",
    background.title = "#4c72b0ff",
    fill="#4c72b0ff",
    type = "coverage")

alTrack_5cap_lo <- AlignmentsTrack(
    "data/5cap_mtec_lo_Aligned.sortedByCoord.out.bam",
    isPaired = TRUE,
    name = "mTEC lo 5Pseq",
    background.title = "#dd8452ff",
    fill="#dd8452ff",
    type = "coverage")

## processed TSR data
data_track_hi <- DataTrack(range = bgFile_hi,
                           genome = "hg38",
                           name = "mTEC hi TSR",
                           col ="#4c72b0ff",
                           background.title = "#4c72b0ff",
                           type = "h")

data_track_lo <- DataTrack(range = bgFile_lo,
                           genome = "hg38",
                           name = "mTEC lo TSR",
                           col ="#dd8452ff",
                           background.title = "#dd8452ff",
                           type = "h")


## Plot Browser Tracks ####
plot_genome <- function(location) {
  plotTracks(
    list(axis_track,
         seq_track,
         gene_track,
         data_track_hi, data_track_lo,
         alTrack_5cap_hi, alTrack_5cap_lo,
         alTrack_hi, alTrack_lo
         ),
    chromosome=as.character(seqnames(location)),
    sizes = c(0.5,1,2,1,1,1,1,2,2),
    from=start(location), to=end(location),
    transcriptAnnotation="symbol")
}

browser_ui <- function(id) {
  ns <- NS(id)
  div(
    fluidRow(
      column(3,
             textInput(ns("location_str"),
                       "Genomic location",
                       "chr1:1225000-1325000")),
      column(9,
             actionButton(ns("go_left"), "<<"),
             actionButton(ns("go_right"), ">>"),
             actionButton(ns("zoom_in"), "+"),
             actionButton(ns("zoom_out"), "-"))
      ),

    plotOutput(ns("genome_plot"), height = "1200px"),

    )
}

browser_server <- function(input, output, session) {
  #To observe input (eg exactly what brushing provides):
  #observe({ print(as.list(input)) })

  location <- reactive({
    GRanges(input$location_str, seqinfo=seqinfo(genome))
  })

  observeEvent(input$go_left, {
    amount <- max(1, width(location()) %/% 4)
    new_location <- shift(location(), -amount)
    new_location <- trim(new_location)
    updateTextInput(session, "location_str", value=as.character(new_location))
  })

  observeEvent(input$go_right, {
    amount <- max(1, width(location()) %/% 4)
    new_location <- shift(location(), amount)
    new_location <- trim(new_location)
    updateTextInput(session, "location_str", value=as.character(new_location))
  })

  observeEvent(input$zoom_in, {
    amount <- width(location()) %/% 4
    new_location <- location()
    start(new_location) <- start(new_location) + amount
    end(new_location) <- end(new_location) - amount
    updateTextInput(session, "location_str", value=as.character(new_location))
  })

  observeEvent(input$zoom_out, {
    amount <- width(location()) %/% 2
    new_location <- location()
    start(new_location) <- start(new_location) - amount
    end(new_location) <- end(new_location) + amount
    new_location <- trim(new_location)
    updateTextInput(session, "location_str", value=as.character(new_location))
  })

  output$genome_plot <- renderPlot({
    plot_genome( location() )
  })
}


