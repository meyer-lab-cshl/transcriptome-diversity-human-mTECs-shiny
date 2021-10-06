library(TxDb.Hsapiens.UCSC.hg38.knownGene)
options(ucscChromosomeNames=FALSE)

txdb_ensembl <- makeTxDbFromEnsembl(organism="Homo sapiens",
                            release=NA,
                            circ_seqs=NULL,
                            server="ensembldb.ensembl.org",
                            username="anonymous", password=NULL, port=0L,
                            tx_attrib=NULL)
saveDb(txdb_ensembl, "data/txdb_ensembl.db")

txdb_ucsc <- makeTxDbFromUCSC("hg38", tablename="knownGene")


saveDb(txdb_ucsc, "data/txdb_ucsc.db")
