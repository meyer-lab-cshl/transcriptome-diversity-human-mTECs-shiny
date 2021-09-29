library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(shinydashboard)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)

#Import data
genome <- BSgenome.Hsapiens.NCBI.GRCh38
seqlevelsStyle(genome) <- "UCSC"
bgFile_hi <- rtracklayer::import.bedGraph("tsr_tpm_hi_avg.bedGraph")

bgFile_lo <- rtracklayer::import.bedGraph("tsr_tpm_lo_avg.bedGraph")

options(ucscChromosomeNames=FALSE)



txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                                    release=NA,
                                    circ_seqs=NULL,
                                    server="ensembldb.ensembl.org",
                                    username="anonymous", password=NULL, port=0L,
                                    tx_attrib=NULL)

seqlevelsStyle(txdb) <- "UCSC"

axis_track <- GenomeAxisTrack()
seq_track <- SequenceTrack(genome)

gene_track <- GeneRegionTrack(
 txdb, genome = "hg38", name="Genes", showId=TRUE, shape="arrow",
 transcriptAnnotation = "symbol")


alTrack_hi <- AlignmentsTrack("mtec_hi_Aligned.sortedByCoord.out.bam",
                           isPaired = TRUE,
                           name = "mTEC hi RNAseq",
                           fill="#4c72b0ff")

alTrack_lo <- AlignmentsTrack("mtec_lo_Aligned.sortedByCoord.out.bam",
                              isPaired = TRUE,
                              name = "mTEC lo RNAseq",
                              fill="#dd8452ff")

alTrack_5cap_hi <- AlignmentsTrack(
  "5cap_mtec_hi_Aligned.sortedByCoord.out.bam",
  isPaired = TRUE,
  name = "mTEC hi 5Pseq",
  fill="#4c72b0ff")

alTrack_5cap_lo <- AlignmentsTrack(
  "5cap_mtec_lo_Aligned.sortedByCoord.out.bam",
  isPaired = TRUE,
  name = "mTEC lo 5Pseq",
  fill="#dd8452ff")


plot_genome <- function(location) {
  # Load data, at a reasonable level of detail for width(location)
  data_track_hi <- DataTrack(range = bgFile_hi, genome = "hg38", type = "h", 
                          name = "mTEC hi TSR",
                          col ="#4c72b0ff",
                          cex.sampleNames = 10)

  
  data_track_lo <- DataTrack(range = bgFile_lo, genome = "hg38", type = "h", 
                             name = "mTEC lo TSR",
                             col ="#dd8452ff",
                             cex.sampleNames = 10)
  
  plotTracks(
    list(axis_track, 
         seq_track, 
         gene_track,
         data_track_hi, 
         data_track_lo),
    chromosome=as.character(seqnames(location)),
    #sizes = c(0.5,1,2,3,3),
    from=start(location), to=end(location),
    transcriptAnnotation="symbol")

}

plot_cap <- function(location) {
  # Load data, at a reasonable level of detail for width(location)

  plotTracks(
    list(alTrack_5cap_hi, alTrack_5cap_lo),
    chromosome=as.character(seqnames(location)),
    from=start(location), to=end(location),
    type = c("coverage"))
}

plot_sashimi <- function(location) {
  # Load data, at a reasonable level of detail for width(location)
  
  plotTracks(
    list(alTrack_hi, alTrack_lo),
    chromosome=as.character(seqnames(location)),
    from=start(location), to=end(location),
    type = c("coverage", "sashimi"))
}


browser_ui <- function(id) {
  ns <- NS(id)
  
  div(
    fluidRow(
      column(3,
             textInput(ns("location_str"), "Genomic location", "chr11:1300000-1400000")),
      column(9,
             actionButton(ns("zoom_in"), "+"),
             actionButton(ns("zoom_out"), "-"))),
    # plotOutput(ns("axis_plot")),
    plotOutput(ns("genome_plot")),
    plotOutput(ns("cap_plot")),
    plotOutput(ns("sashimi_plot")))
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
  
  output$cap_plot <- renderPlot({
    plot_cap( location() )
  })
  
  output$sashimi_plot <- renderPlot({
    plot_sashimi( location() )
  })
}


