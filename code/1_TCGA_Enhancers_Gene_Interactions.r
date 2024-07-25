library(GenomicRanges)
library(InteractionSet)
library(MultiAssayExperiment)
library(RaggedExperiment)
library(dplyr)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(RColorBrewer)
library(ggrepel)
library(ggExtra)
library(tidyr)
library(ggridges)

MAEobject <- readRDS('intermediates/TCGA_MAEobject.rds')
MAEobject

#granges with all genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
seqlevelsStyle(genes_gr) <- 'UCSC'
genes_gr$Gene <- mapIds(org.Hs.eg.db, genes_gr$gene_id, "SYMBOL","ENTREZID")
genes_gr

gene_TSS <- data.frame(genes_gr) %>% 
mutate(TSS_start=ifelse(strand=='-', end,start)) 

gene_TSS_gr <- makeGRangesFromDataFrame(
    gene_TSS[,c('seqnames','TSS_start','strand','gene_id','Gene')],
    start.field = 'TSS_start', end.field = 'TSS_start',
    keep.extra.columns = TRUE)
names(gene_TSS_gr) <- gene_TSS_gr$gene_id
gene_TSS_gr

saveRDS(genes_gr, 'intermediates/genes_gr.rds')
saveRDS(gene_TSS_gr, 'intermediates/gene_TSS_gr.rds')

#get all peaks within 1Mb from each gene
overlaps <- findOverlaps(gene_TSS_gr, experiments(MAEobject)$normcounts_hichip, 
                         maxgap=1000000, ignore.strand=TRUE) 
#distance between enhancer-gene
overlaps

counts_expanded <- assays(MAEobject)$normcounts_hichip[subjectHits(overlaps),]
head(counts_expanded)
nrow(counts_expanded)

hichip_genes_gi <- GInteractions(gene_TSS_gr[queryHits(overlaps)], 
                                 rowRanges(experiments(MAEobject)$normcounts_hichip)[subjectHits(overlaps)])
hichip_genes_gi

merged_loops <- read.delim("intermediates/TCGA_HiChIP/loops/Merged.FitHiChIP.interactions_Q0.1_MergeNearContacts.bed")
merged_loops %>% head()

merged_loops_left <- makeGRangesFromDataFrame(merged_loops, 
                                              seqnames.field = 'chr1', start.field = 's1', end.field = 'e1')
merged_loops_right <- makeGRangesFromDataFrame(merged_loops, 
                                              seqnames.field = 'chr2', start.field = 's2', end.field = 'e2')
merged_loops_gi <- GInteractions(merged_loops_left, merged_loops_right)
merged_loops_gi

saveRDS(merged_loops_gi, 'intermediates/merged_loops_gi.rds')

#identify hichip peak-gene pairs supported by hichip loops
hichip_genes_loops <- findOverlaps(hichip_genes_gi, merged_loops_gi, ignore.strand=TRUE)
hichip_genes_loops

unique_gene_peak_hits <- queryHits(hichip_genes_loops) %>% unique() 
counts_gene_peaks <- counts_expanded[unique_gene_peak_hits,]
rownames(counts_gene_peaks) <- NULL
head(counts_gene_peaks)
nrow(counts_gene_peaks)
hichip_gene_peaks_gi <- hichip_genes_gi[unique_gene_peak_hits,]
hichip_gene_peaks_gi

colData <- colnames(counts_gene_peaks) %>% DataFrame()
colnames(colData) <- 'primarynames_hichip_gene_peaks'
colData
hichip_gene_peaks_iset <- InteractionSet(list(counts=counts_gene_peaks), hichip_gene_peaks_gi, colData=colData)
hichip_gene_peaks_iset

#append gene-peak interactionSet to MAEobject
MAEobject2 <- c(MAEobject, hichip_gene_peaks_iset=hichip_gene_peaks_iset)
MAEobject2
sampleMap(MAEobject2)

saveRDS(MAEobject2, 'intermediates/TCGA_MAEobject2.rds')
