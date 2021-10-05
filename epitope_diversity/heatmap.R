
## libraries ####
library("DESeq2")
library("cowplot")
library("limma")
library("ggrepel")
library("tidyverse")
library("scales")
library("ggnewscale")
library("ComplexHeatmap")
library("circlize")
library("EnhancedVolcano")
library("readr")
library("janitor")
library("biomaRt")

## functions ####
theme_set(theme_half_open(font_family = 'sans' ))

title_size <- 10
text_size <- 8
point_size <- 0.8

## settings ####
map_color_range <- function(matrix, color_vec) {
  matrix %>%
    range(., na.rm=TRUE) %>%
    abs %>%
    max %>%
    seq(-., ., length.out = length(color_vec)) %>%
    colorRamp2(., color_vec)
}

## data ####
diffexp <- read.csv(file = "epitope_data/mTEC_RNAseq_results_all.csv")

#keep only significant genes

diffexp <- diffexp[!(is.na(diffexp$pval)),] 
diffexp <- diffexp[diffexp$qval < 0.05,]
rownames(diffexp) <- NULL

#remove rows with duplicated ens_gene ids
diffexp_no_dup <- diffexp[!duplicated(diffexp$ens_gene),]

diffexp_no_dup <- diffexp_no_dup %>%
  drop_na(ens_gene)

#Switch from ensembl id to gene symbol
ensembl_human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genemap_human <- getBM(attributes=c("ensembl_gene_id",
                                               "external_gene_name"),
                                  filters = "ensembl_gene_id",
                                  values = diffexp_no_dup$ens_gene,
                                  mart = ensembl_human) 

#Remove genemap entries with duplicates in external_gene_name
genemap_human <- genemap_human[!duplicated(genemap_human$external_gene_name),]

rownames(genemap_human) <- NULL

#Get rid of genes that don't have entries in genemap_human
diffexp_no_dup <- diffexp_no_dup[(diffexp_no_dup$ens_gene 
                                                    %in% genemap_human$ensembl_gene_id), ]
rownames(diffexp_no_dup) <- NULL

diffexp_no_dup <- diffexp_no_dup %>%
  column_to_rownames(var = "ens_gene")

genemap_human <- genemap_human %>%
  column_to_rownames(var = "ensembl_gene_id") 

diffexp_no_dup_symbols <- cbind(diffexp_no_dup, genemap_human) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  column_to_rownames(var = "external_gene_name")


# obtain tpms
counts <- diffexp_no_dup_symbols %>%
  dplyr::select(c("pt214_hi_tpm", 
                  "pt214_lo_tpm",
                  "pt221_hi_tpm",
                  "pt221_lo_tpm",
                  "pt226_hi_tpm",
                  "pt226_lo_tpm")) 

#caclualte tpm z score
counts_zscore <- t(counts) %>%
  scale %>%
  t %>%
  as.data.frame()


fc <- diffexp_no_dup_symbols %>%
  dplyr::select("b") 

#Combine tpms and fc

combined <- cbind(counts_zscore, fc)


#counts matrix
counts_matrix <- combined %>%
  dplyr::select(-c(b)) %>%
  dplyr::select(contains("hi"), contains("lo")) %>%
  as.matrix

fc_matrix <- combined %>%
  dplyr::select(b) %>%
  as.matrix


#Reformat TPM's and create coldata
labels <- gsub('pt\\d*_([hilo]{2})_tpm', '\\1', colnames(counts_matrix))
sample <- colnames(counts_matrix)
coldata <- data.frame(sample, labels)






# Set heatmap colors, themes, and sizes
cond_type_vals <- c('#4c72b0ff','#dd8452ff')
names(cond_type_vals) <- unique(coldata$labels)
color <- list(cond_type = cond_type_vals)
counts_color <- c('#a6611a','#dfc27d','#f5f5f5','#80cdc1','#018571')
#counts_color <- c('#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494')
fc_color <- c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020')


# Create column annotations for rnaseq data
coldata <- coldata %>%
  as_tibble

columnAnno <- columnAnnotation(
  cond_type = coldata$labels,
  col = list(cond_type = color$cond_type),
  annotation_legend_param = list(cond_type = list(title = "MHCII")),
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = text_size)
)

# Create the heatmap


hm_counts_sig <- Heatmap(counts_matrix,
                         name="TPM zscore",
                         col=map_color_range(counts_matrix, counts_color),
                         cluster_rows=TRUE,
                         cluster_columns=FALSE,
                         show_row_dend = FALSE,
                         show_row_names=FALSE,
                         show_column_names=TRUE,
                         top_annotation = columnAnno,
                         row_names_gp = gpar(fontsize = text_size),
                         column_names_gp = gpar(fontsize = text_size),
                         heatmap_legend_param =
                           list(direction = "horizontal",
                                title_gp = gpar(fontsize = text_size,
                                                fontface = "bold"),
                                labels_gp = gpar(fontsize = text_size)),
                         use_raster=TRUE,
                         row_title=NULL,
                         column_title=NULL)

hm_fc <- Heatmap(fc_matrix,
                 col=map_color_range(fc_matrix, fc_color),
                 use_raster=TRUE,
                 show_row_names=FALSE,
                 show_row_dend = FALSE,
                 show_column_names = FALSE,
                 row_title = NULL,
                 column_title = NULL,
                 cluster_columns=FALSE,
                 cluster_rows=TRUE,
                 row_names_gp = gpar(fontsize = text_size),
                 column_names_gp = gpar(fontsize = text_size),
                 heatmap_legend_param =
                   list(direction = "horizontal",
                        title_gp = gpar(fontsize = text_size,
                                        fontface = "bold"),
                        labels_gp = gpar(fontsize = text_size)),
                 name = "lFC",
                 width = unit(3, "cm"))


#Save counts heatmap
saveRDS(hm_counts_sig, file = "tpm_heatmap.rds")
saveRDS(hm_fc, "fc_heatmap.rds")

