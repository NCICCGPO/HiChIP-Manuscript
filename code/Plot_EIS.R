library(rtracklayer)
library(ArchR)
library(plyranges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Organism.dplyr)
library(Homo.sapiens)
library(dplyr)

source("plotting_functions.R")

#-----------------
# Plot EIS and H3K27ac signal tracks (MYC locus, COAD & LIHC, Figure 1e)
#-----------------

#https://reftss.riken.jp/datafiles/3.3/human/refTSS_v3.3_human_coordinate.hg38.bed.gz
tss <- 'refTSS_v3.3_human_coordinate.hg38.bed' 

# Plotting MYC locus
region <- GRanges(seqnames = "chr8", ranges = IRanges(126200000, 129000000)) 

# Liver H3K27ac ChIP-seq from ENCODE https://www.encodeproject.org/files/ENCFF905FLR/@@download/ENCFF905FLR.bigWig
bigwig1 <- import('ENCFF905FLR.bigWig', format = "BigWig")
bigwig1$RG <- "Liver H3K27ac ChIP-seq"
# Colon H3K27ac ChIP-seq from ENCODE https://www.encodeproject.org/files/ENCFF873MWG/@@download/ENCFF873MWG.bigWig
bigwig2 <- import('ENCFF873MWG.bigWig', format = "BigWig")
bigwig2$RG <- "Colon H3K27ac ChIP-seq"
bigwig <- c(bigwig1, bigwig2)

# Plot H3K27ac ChIP-seq signal
p1 <- bigWig_track(bigwig_gr = bigwig, region = region, tss = tss)

# COAD H3K27ac HiChIP 1D signal bed from publication page
fragment1 <- import("COAD_H3K27ac_1D-signal.sort.bed.gz", format = "bed")
fragment1$RG <- "COAD"
# LIHC H3K27ac HiChIP 1D signal bed from publication page
fragment2 <- import("LIHC_H3K27ac_1D-signal.sort.bed.gz", format = "bed")
fragment2$RG <- "LIHC"
fragment <- c(fragment1, fragment2)

# Plot H3K27ac HiChIP 1D signal
p2 <- H3K27Ac_1D_track(fragment, region = region, tss = tss)

# H3K27ac HiChIP hic files from publication page
hic_files_merged <- c("COAD_H3K27ac.allValidPairs.hic", "LIHC_H3K27ac.allValidPairs.hic")
names(hic_files_merged) <- gsub("_.*", "",hic_files_merged)

# number of valid pairs
validPairs_merged <- setNames(c(221818488, 128424571), c("COAD", "LIHC"))

# Plot EIS and loops
p3 <- virtual4C_plot(hic_files_merged[c('COAD')]
                     , names = names(hic_files_merged[c('COAD')])
                     , res = 10000
                     , color = paletteDiscrete(values = c("COAD", "LIHC"))['COAD']
                     , norm = validPairs_merged[c('COAD')]
                     , bedpe = "COAD_H3K27ac.FitHiChIP.interactions_Q0.1_MergeNearContacts.bedpe", 
                     , region = region
                     , anchor_gene = "MYC")

p4 <- virtual4C_plot(hic_files_merged[c('LIHC')]
               , pal = paletteDiscrete(values = c("COAD", "LIHC"))
               , names = names(hic_files_merged[c('LIHC')])
               , res = 10000
               , color = paletteDiscrete(values = c("COAD", "LIHC"))['LIHC']
               , norm = validPairs_merged[c('COAD')]
               , bedpe = "LIHC_H3K27ac.FitHiChIP.interactions_Q0.1_MergeNearContacts.bedpe", 
               , region = region
               , anchor_gene = "MYC")

# Plot genes
p5 <- plot_genes(region = region, gene_list = "MYC")

# Combine plots
patchwork::wrap_plots(p1,p2,p3,p4,p5, ncol = 1, heights = c(4,4,4,4,1))