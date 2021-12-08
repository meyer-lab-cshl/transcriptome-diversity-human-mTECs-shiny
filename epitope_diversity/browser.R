## Libraries ####
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(shinydashboard)
library(InteractiveComplexHeatmap)

options(shiny.port = 6372)
options(shiny.host = "0.0.0.0")
options(shiny.sanitize.errors = TRUE)

data_path <- "../data"

## Import data objects ####
txdb <- loadDb(file.path(data_path, "txdb_ucsc.db"))

bgFile_hi <- rtracklayer::import.bedGraph(file.path(data_path,
                                                    "tsr_tpm_hi_sum.bedGraph"))
bgFile_lo <- rtracklayer::import.bedGraph(file.path(data_path,
                                                    "tsr_tpm_lo_sum.bedGraph"))

ht_tpm <- readRDS(file.path(data_path, "tpm_heatmap.rds"))
ht_lfc <- readRDS(file.path(data_path, "fc_heatmap.rds"))

ht <- ht_tpm + ht_lfc

##Load locations for transcript id's and gene names in txdb
tx_locs <- readRDS(file.path(data_path, "tx_locs.rds"))

find_loc <- function(term) {
  term <- if(startsWith(term, "chr")) {
    term
  } else {
    toupper(term)
  }

  index <- if(term %in% tx_locs$ensembl_transcript_id_version) {
    grep(term, tx_locs$ensembl_transcript_id_version)
  } else if(term %in% tx_locs$ensembl_transcript_id) {
    grep(term, tx_locs$ensembl_transcript_id)
  } else if(term %in% tx_locs$external_gene_name) {
    grep(term, tx_locs$external_gene_name)
  } else if(startsWith(term, "chr")) {
    term
  } else {
    "error"
  }

  if(startsWith(term, "chr")) {
    loc <- term
  } else if(index == "error") {
    loc <- "error"
  } else {loc <- paste0("chr",
                        tx_locs$chromosome_name[index],
                        ":",
                        tx_locs$start_position[index],
                        "-",
                        tx_locs$end_position[index])
  }
  
  return(unique(loc))
}


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

## RNAseq data
alTrack_hi <- AlignmentsTrack(
  #file.path(data_path, "5cap_mtec_hi_Aligned.sortedByCoord.out.bam"),
  file.path(data_path, "mtec_hi_Aligned.sortedByCoord.out.bam"),
  isPaired = TRUE,
  name = "mTEC hi RNAseq",
  fill="#4c72b0ff",
  type = c("coverage", "sashimi"))

alTrack_lo <- AlignmentsTrack(
  #file.path(data_path, "5cap_mtec_hi_Aligned.sortedByCoord.out.bam"),
  file.path(data_path, "mtec_lo_Aligned.sortedByCoord.out.bam"),
  isPaired = TRUE,
  name = "mTEC lo RNAseq",
  fill="#dd8452ff",
  type = c("coverage", "sashimi"))

## 5'Cap data
alTrack_5cap_hi <- AlignmentsTrack(
  file.path(data_path, "5cap_mtec_hi_Aligned.sortedByCoord.out.bam"),
  isPaired = TRUE,
  name = "mTEC hi 5Pseq",
  fill="#4c72b0ff",
  type = "coverage")

alTrack_5cap_lo <- AlignmentsTrack(
  file.path(data_path, "5cap_mtec_lo_Aligned.sortedByCoord.out.bam"),
  isPaired = TRUE,
  name = "mTEC lo 5Pseq",
  fill="#dd8452ff",
  type = "coverage")

## processed TSR data
data_track_hi <- DataTrack(range = bgFile_hi,
                           genome = "hg38",
                           name = "mTEC hi TSR",
                           col ="#4c72b0ff",
                           type = "h")

data_track_lo <- DataTrack(range = bgFile_lo,
                           genome = "hg38",
                           name = "mTEC lo TSR",
                           col ="#dd8452ff",
                           type = "h")

## Plot Browser Tracks ####
plot_genome <- function(location) {
  chr_name <- as.character(seqnames(location))
  ideogram_track <- IdeogramTrack(genome = "hg38", chromosome = chr_name)
  
  gene_track <- GeneRegionTrack(
    txdb, genome = "hg38", name="Genes", showId=TRUE, shape="arrow",
    transcriptAnnotation = "symbol", chromosome = chr_name,
    start=start(location), end=end(location))
  
  plotTracks(
    list(ideogram_track,
         axis_track,
         gene_track,
         data_track_hi, data_track_lo,
         alTrack_5cap_hi, alTrack_5cap_lo,
         alTrack_hi, alTrack_lo
    ),
    chromosome=chr_name,
    sizes = c(0.1,
              0.5,
              2.5,
              0.8,
              0.8,
              1,
              1,
              2,
              2),
    from=start(location), to=end(location),
    background.title = "white",
    col.title = "black",
    col.axis = "black",
    transcriptAnnotation="symbol")
}

browser_ui <- function(id) {
  ns <- NS(id)
  div(
    div(textInput(
      ns("location_str"),
      "Genomic Location, Gene Name, or Ensembl Transcript ID",
      "chr1:1250000-1300000"
    ),
    style = "margin-top:25px"),
    div(
      actionButton(ns("go_left"), "<<"),
      actionButton(ns("go_right"), ">>")
    ),
    plotOutput(ns("genome_plot"), height = "1200px"),
  )
}

browser_server <- function(input, output, session) {
  
  location <- reactive({
    GRanges(find_loc(input$location_str), seqinfo=seqinfo(genome))
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

  output$genome_plot <- renderPlot({
    if (find_loc(input$location_str) == "error") {
      stop(safeError(
        "Input unknown: please check for typos"
      ))
    } else {
      plot_genome( location() )
    }
  })
}


