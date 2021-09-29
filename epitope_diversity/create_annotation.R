library(TxDb.Hsapiens.UCSC.hg38.knownGene)
options(ucscChromosomeNames=FALSE)

txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                            release=NA,
                            circ_seqs=NULL,
                            server="ensembldb.ensembl.org",
                            username="anonymous", password=NULL, port=0L,
                            tx_attrib=NULL)
seqlevelsStyle(txdb) <- "UCSC"
saveDb(txdb, "epitope_diversity/data/txdb.db")
