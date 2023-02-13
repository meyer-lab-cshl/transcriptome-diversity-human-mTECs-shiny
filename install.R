#.libPaths("/opt/R/4.2.1/lib/R/library")
#.libPaths("/bin/R")

getRversion()


## CRAN ####
install.packages(c('devtools', 'magick', 'BiocManager', 'data.table',
                   'shinydashboard', 'RMariaDB', 'shiny', 'DT', 'berryFunctions'),
                 repos = "http://cran.us.r-project.org")

## Bioconductor

###Install specific BiocManager version (3.14)
#BiocManager::install()
#BiocManager::install(version = "3.16")
###Install packages
###Install packages
BiocManager::install(c('GenomicRanges', 'Gviz', 'rtracklayer',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'BSgenome.Hsapiens.NCBI.GRCh38',
		       'ComplexHeatmap', 'InteractiveComplexHeatmap'), force = TRUE
)

BiocManager::install("Bioconductor/GenomeInfoDb")

## github ####
# currently from github due to: https://www.biostars.org/p/9490625/
#devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/InteractiveComplexHeatmap")
