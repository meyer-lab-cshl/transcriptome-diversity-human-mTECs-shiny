## CRAN ####
install.packages(c('devtools', 'magick', 'BiocManager', 'data.table',
                   'shinydashboard', 'RMariaDB'))

## Bioconductor
BiocManager::install(c('GenomicRanges', 'Gviz', 'rtracklayer',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'BSgenome.Hsapiens.NCBI.GRCh38',
                       'ComplexHeatmap')
)
## github ####
# currently from github due to: https://www.biostars.org/p/9490625/
devtools::install_github("jokergoo/InteractiveComplexHeatmap")