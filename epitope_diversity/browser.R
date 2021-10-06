## Libraries ####
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(shinydashboard)
library(InteractiveComplexHeatmap)

##Test comment

data_path <- "../data"

## Import data objects ####
if (TRUE) {
    options(ucscChromosomeNames=FALSE)
    txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                                release=NA,
                                circ_seqs=NULL,
                                server="ensembldb.ensembl.org",
                                username="anonymous", password=NULL, port=0L,
                                tx_attrib=NULL)
    seqlevelsStyle(txdb) <- "UCSC"
}
#txdb <- loadDb(file.path(data_path, "txdb.db"))

bgFile_hi <- rtracklayer::import.bedGraph(file.path(data_path, "tsr_tpm_hi_sum.bedGraph"))
bgFile_lo <- rtracklayer::import.bedGraph(file.path(data_path, "tsr_tpm_lo_sum.bedGraph"))

ht_tpm <- readRDS(file.path(data_path, "tpm_heatmap.rds"))
ht_lfc <- readRDS(file.path(data_path, "fc_heatmap.rds"))

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
gene_track <- GeneRegionTrack(
    txdb, genome = "hg38", name="Genes", showId=TRUE, shape="arrow",
    transcriptAnnotation = "symbol")

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
            "Genomic location",
            "chr1:1250000-1300000"
        ),
        style = "margin-top:25px"),
        div(
            actionButton(ns("go_left"), "<<"),
            actionButton(ns("go_right"), ">>"),
            actionButton(ns("zoom_in"), "+"),
            actionButton(ns("zoom_out"), "-")
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


