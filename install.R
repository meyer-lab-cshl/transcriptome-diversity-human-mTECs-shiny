## CRAN ####
install.packages(c('devtools', 'magick', 'BiocManager', 'data.table',
                   'shinydashboard', 'RMariaDB', 'shiny', 'DT', 'shinyWidgets'), 
                 repos = "http://cran.us.r-project.org")

## Bioconductor

###Install specific BiocManager version (3.14)
BiocManager::install(version = "3.14")
###Install packages
BiocManager::install(c('GenomicRanges', 'Gviz', 'rtracklayer',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'BSgenome.Hsapiens.NCBI.GRCh38',
                       'ComplexHeatmap', 'GenomeInfoDb')
)
## github ####
# currently from github due to: https://www.biostars.org/p/9490625/
devtools::install_github("jokergoo/InteractiveComplexHeatmap")