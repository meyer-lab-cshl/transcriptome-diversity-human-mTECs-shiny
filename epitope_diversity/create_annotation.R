#.libPaths("/opt/R/4.2.1/lib/R/library")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)

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


##Get locations for transcript id's in txdb
tx_names <- mcols(transcripts(txdb_ucsc, "tx_name"))

ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")

tx_locs <- getBM(attributes=c('ensembl_transcript_id_version',
                              'ensembl_transcript_id',
                              'external_gene_name',
                              'chromosome_name',
                              'start_position',
                              'end_position'),
                 values=tx_names$tx_name,
                 filters = "ensembl_transcript_id_version",
                 mart=ensembl)

saveRDS(tx_locs, file = "data/tx_locs.rds")


