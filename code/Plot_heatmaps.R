library(CALDER)
library(tidyr)
library(jgplot2)
library(ggpubr)
library(pheatmap)
library(dendextend)
library(NMF)
if(any(grepl("package:plyr", search()))) detach("package:plyr")
source("tcga_palettes.R")

#-----------------
# Unsupervised hierarchical clustering based on HiChIP sub-compartments (Figure 1f)
#-----------------

# Call subcompartments with CALDER using hic files (on publication page)
hic_files <- list.files("TCGA_HiChIP_hic", pattern = "*.hic$", full.names = TRUE)
chr <- 1:22
res <- 100E3

for (hic in hic_files) {
    CALDER(contact_file_hic=hic, 
           chrs=chr, 
           bin_size=res,
           genome='hg38',
           save_dir=paste0("CALDER/",gsub(".allValidPairs.hic$", "", basename(hic))),
           save_intermediate_data=TRUE,
           n_cores=20,
           sub_domains=TRUE)
}

samples <- list.dirs("CALDER", recursive = FALSE)
samples <- samples[grepl("T1_H3K27ac", samples)]
subcompartments_all <- data.frame()
for (i in samples) {
    subcompartments <- read.delim(paste0(i, "/sub_compartments/all_sub_compartments.tsv"))
    subcompartments$bin <- paste(subcompartments$chr, subcompartments$pos_start, subcompartments$pos_end
                                 , sep = "_")
    subcompartments$sample <- basename(i)
    subcompartments_all <- rbind(subcompartments_all, subcompartments)
}

subcompartments_all_spread <- spread(subcompartments_all[,c('comp_rank', 'bin', 'sample')]
                                     , key = sample, value = comp_rank)
rownames(subcompartments_all_spread) <- subcompartments_all_spread$bin
subcompartments_all_spread <- subcompartments_all_spread[,2:ncol(subcompartments_all_spread)]
colnames(subcompartments_all_spread) <- gsub("-", "_", colnames(subcompartments_all_spread))

# Compute pairwise pearson correlation using subcompartment rank
subcompartments_all_cor <- cor(subcompartments_all_spread, use = "pairwise.complete.obs")

# Plot heatmap
heatmap_subcompartment_cor <- pheatmap(subcompartments_all_cor, silent = TRUE)
colData <- data.frame(submitter_ID = colnames(subcompartments_all_cor))
colData$project_ID <- gsub("x", "", gsub("_.*", "", colData$submitter_ID))
rownames(colData) <- colData$submitter_ID
colData <- colData[,'project_ID',drop = FALSE]
colors <- tcga_discrete_palettes$cancer_types[
    names(tcga_discrete_palettes$cancer_types) %in% colData$project_ID]
tree_row <- heatmap_subcompartment_cor$tree_row %>% sort()
tree_col <- heatmap_subcompartment_cor$tree_col %>% sort()
pheatmap(subcompartments_all_cor, annotation_legend = FALSE,cluster_rows = tree_row
         , cluster_cols = tree_col, annotation_col = colData, show_rownames = FALSE
         , fontsize = 12, border_color = NA
         , annotation_row = colData
         , treeheight_row = 0
         , annotation_names_col = FALSE
         , annotation_names_row = FALSE
         , annotation_colors = list(project_ID = colors)
         , color = colorRampPalette(create_pal_c())(100)
         , breaks = seq(0,0.8, length.out = 100)
         , show_colnames = FALSE)

# Calculate clustering purity and entropy
NMF::purity(factor(colData$project_ID)
            , factor(cutree(tree_row , k = length(unique(colData$project_ID)))))
NMF::entropy(factor(colData$project_ID)
            , factor(cutree(tree_row , k = length(unique(colData$project_ID)))))