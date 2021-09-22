library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(AnnotationDbi)


# Ideally these globals would be wrapped in some sort of object...

keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

#genome <- BSgenome.Hsapiens.UCSC.hg38
genome <- BSgenome.Hsapiens.NCBI.GRCh38
seqlevelsStyle(genome) <- "UCSC"
# bgFile <- rtracklayer::import.bedGraph("pt214_hi_Aligned.sortedByCoord.out.bedGraph")
# 
# sequences_to_keep <- levels(seqnames(bgFile))
# genome <- keepBSgenomeSequences(genome, sequences_to_keep)


#txdb.human <- TxDb.Hsapiens.UCSC.hg38.knownGene

#genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
options(ucscChromosomeNames=FALSE)



txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                                    release=NA,
                                    circ_seqs=NULL,
                                    server="ensembldb.ensembl.org",
                                    username="anonymous", password=NULL, port=0L,
                                    tx_attrib=NULL)

seqlevelsStyle(txdb) <- "UCSC"

# txdb <- loadDb(file = "txdb.human.db")
axis_track <- GenomeAxisTrack()
seq_track <- SequenceTrack(genome)

gene_track <- GeneRegionTrack(
 txdb, genome = "hg38", name="Genes", showId=TRUE, shape="arrow")

alTrack <- AlignmentsTrack("pt214_hi_Aligned.sortedByCoord.out.bam",
                           isPaired = TRUE)

#dTrack2 <- DataTrack(range = bgFile, genome = "mm10", type = "a", 
                    # name = "bedGraph")


#gene_track <- DataTrack(range = bgFile, genome = "mm10", type = "a", 
                     #name = "bedGraph")

  
plot_genome <- function(location) {
  # Load data, at a reasonable level of detail for width(location)
  data_track <- DataTrack(range = bgFile, genome = "hg38", type = "h", 
                          name = "bedGraph")
  
  plotTracks(
    list(axis_track, seq_track, gene_track, data_track),
    chromosome=as.character(seqnames(location)),
    from=start(location), to=end(location))
  
  # plotTracks(
  #   list(gene_track, alTrack),
  #   chromosome=as.character(seqnames(location)),
  #   from=start(location), to=end(location),
  #   type = c("coverage", "sashimi"))
}

plot_sashimi <- function(location) {
  # Load data, at a reasonable level of detail for width(location)
  
  plotTracks(
    list(alTrack),
    chromosome=as.character(seqnames(location)),
    from=start(location), to=end(location),
    type = c("coverage", "sashimi"))
}

browser_ui <- function(id) {
  ns <- NS(id)
  
  div(
    fluidRow(
      column(3,
             textInput(ns("location_str"), "Genomic location", "chr11:140000-180000")),
      column(9,
             actionButton(ns("go_left"), "<<"),
             actionButton(ns("go_right"), ">>"),
             actionButton(ns("zoom_in"), "+"),
             actionButton(ns("zoom_out"), "-"))),
    plotOutput(ns("genome_plot")),
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
  
  output$sashimi_plot <- renderPlot({
    plot_sashimi( location() )
  })
}


