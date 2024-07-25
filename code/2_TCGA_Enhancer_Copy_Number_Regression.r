library(dplyr)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(RaggedExperiment)
library(InteractionSet)
library(gridExtra)
library(umap)
library(ggplot2)

source('TCGA_King_misc_functions.r')
source('tcga_palettes.R')

MAEobject2 <- readRDS('intermediates/TCGA_MAEobject2.rds')
MAEobject2

MAEobject2M_hichip_CN <- intersectColumns(MAEobject2[ , , c('segmentCN_WGS_ploidy_corrected',
                                                         'normcounts_hichip',
                                                         'hichip_gene_peaks_iset')])
MAEobject2M_hichip_CN

#create CN matrix for 1D peaks
timecheck_1 <- Sys.time()
samples <- colData(MAEobject2M_hichip_CN)$submitter_ID
enhancerCN <- lapply(samples, function(sample_name){
    segmentCN_df <- assays(MAEobject2M_hichip_CN)$segmentCN_WGS_ploidy_corrected %>% data.frame() 
    idx <- grep(TRUE, !is.na(segmentCN_df[,sample_name]))
    overlaps <- findOverlaps(rowRanges(experiments(MAEobject2M_hichip_CN)$segmentCN_WGS_ploidy_corrected[idx]), 
                             rowRanges(experiments(MAEobject2M_hichip_CN)$normcounts_hichip))
        
    enhancerCN_sample <- assays(MAEobject2M_hichip_CN)$segmentCN_WGS_ploidy_corrected[idx,sample_name][queryHits(overlaps)] %>% 
        data.frame(CN=.)
    enhancerCN_sample$PeakID <- 
        as.character(rowData(experiments(MAEobject2M_hichip_CN)$normcounts_hichip)$PeakID)[subjectHits(overlaps)]
    enhancerCN_sample <- enhancerCN_sample %>%
        group_by(PeakID) %>%
        dplyr::summarize(CN=mean(CN))
    colnames(enhancerCN_sample)[colnames(enhancerCN_sample)=='CN'] <- sample_name
    return(enhancerCN_sample)
})
enhancerCN <- Reduce(function(x,y){merge(x,y,by='PeakID')}, enhancerCN) %>% data.frame()
peak_info <- data.frame(rowRanges(experiments(MAEobject2M_hichip_CN)$normcounts_hichip))
enhancerCN <- merge(peak_info, enhancerCN, by='PeakID')
rownames(enhancerCN) <- enhancerCN$PeakID
timecheck_2 <- Sys.time()
timecheck_2 - timecheck_1

head(enhancerCN)

primarynames_enhancerCN <- colnames(enhancerCN) %>% subset(grepl('^TCGA',.))
counts <- as.matrix(enhancerCN[ , primarynames_enhancerCN]) 
head(counts)
rowRanges <- makeGRangesFromDataFrame(enhancerCN[,1:6], keep.extra.columns = TRUE)
rowRanges
colData <- DataFrame(primarynames_enhancerCN)
colData

enhancerCN_summarized <- SummarizedExperiment(assays=list(counts=counts), rowRanges=rowRanges, colData=colData)
enhancerCN_summarized

cn_matrix <- assay(enhancerCN_summarized)
#old array CN data: enhancerCN values are log2(cn/2), where diploid mean is zero (log2(2/2)). 
#cn_matrix <- (2^cn_matrix)*2
hichip_matrix <- assays(MAEobject2)$normcounts_hichip[rownames(cn_matrix) , colnames(cn_matrix)]
hichip_matrix_cn_regressed <- hichip_matrix / (cn_matrix*2 + 1)
head(hichip_matrix_cn_regressed)

primarynames_hichip_cn_regressed <- colnames(hichip_matrix_cn_regressed) %>% subset(grepl('^TCGA',.))
counts <- as.matrix(hichip_matrix_cn_regressed[ , primarynames_hichip_cn_regressed]) 
head(counts)
rowRanges <- rowRanges( enhancerCN_summarized )
rowRanges
colData <- DataFrame(primarynames_hichip_cn_regressed)
colData

hichip_cn_regressed_summarized <- SummarizedExperiment(assays=list(counts=counts), 
                                                       rowRanges=rowRanges, 
                                                       colData=colData)
hichip_cn_regressed_summarized

iset_gi <- interactions(experiments(MAEobject2)$hichip_gene_peaks_iset) %>%
    subset(anchor2.PeakID %in% rownames(hichip_matrix_cn_regressed))
iset_gi

iset_PeakIDs <- iset_gi$anchor2.PeakID %>% as.character()

counts <- hichip_matrix_cn_regressed[iset_PeakIDs , ]
rownames(counts) <- NULL
head(counts)

colData <- colnames(counts) %>% DataFrame()
colnames(colData) <- 'primarynames_hichip_gene_peaks_cn_regressed'
hichip_gene_peaks_iset_cn_regressed <- InteractionSet(list(counts=counts), iset_gi, colData=colData)
hichip_gene_peaks_iset_cn_regressed

normcounts_fithichip_loops <- read.delim(
    'intermediates/Merged.FitHiChIP.interactions_Q0.1_MergeNearContacts_normcounts.txt')
colnames(normcounts_fithichip_loops) <- colnames(normcounts_fithichip_loops) %>% 
    gsub('_','-',.) %>% 
    LibraryNames2Submitters(sampleid)
head(normcounts_fithichip_loops)

leftbin <- makeGRangesFromDataFrameInColumnOrder(normcounts_fithichip_loops[,1:3])
rightbin <- makeGRangesFromDataFrameInColumnOrder(normcounts_fithichip_loops[,4:6])
gi <- GInteractions(leftbin, rightbin)
gi

counts <- normcounts_fithichip_loops[,8:ncol(normcounts_fithichip_loops)] %>% as.matrix()
colData <- DataFrame(primaryname=colnames(counts))

normcounts_fithichip_loops_iset <- InteractionSet(counts, gi, colData=colData)
normcounts_fithichip_loops_iset

MAEobject2M_loopcounts <- c(MAEobject2, 
                            normcounts_fithichip_loops = normcounts_fithichip_loops_iset)

MAEobject2M_loopcounts <- intersectColumns(
    MAEobject2M_loopcounts[ , , c('segmentCN_WGS_ploidy_corrected','normcounts_fithichip_loops')]
)
MAEobject2M_loopcounts

samples <- colData(MAEobject2M_loopcounts)$submitter_ID
loopcounts_iset <- experiments(MAEobject2M_loopcounts)$normcounts_fithichip_loops
segmentCN <- experiments(MAEobject2M_loopcounts)$segmentCN_WGS_ploidy_corrected

timecheck_1 <- Sys.time()

left_anchors <- anchors(interactions(loopcounts_iset))$first
right_anchors <- anchors(interactions(loopcounts_iset))$second                        
loopcountsCN <- lapply(samples, function(sample_name){
    segmentCN_sample <- experiments(MAEobject2M_loopcounts)$segmentCN_WGS_ploidy_corrected[,sample_name]
    segmentCN_df <- assay(segmentCN_sample) %>% data.frame() 
    idx <- grep(TRUE, !is.na(segmentCN_df))
    segmentCN_sample <- segmentCN_sample[idx , ]
    overlaps_left <- findOverlaps(rowRanges(segmentCN_sample), left_anchors)
    overlaps_right <- findOverlaps(rowRanges(segmentCN_sample), right_anchors)    
    
    overlaps_merged <- merge(data.frame(overlaps_left), data.frame(overlaps_right), by = 'subjectHits') %>%
        dplyr::rename(loopcounts_idx = subjectHits, 
                      left_segmentCN_idx = queryHits.x, 
                      right_segmentCN_idx = queryHits.y)

    segmentCN_anchors <- data.frame(loopcounts_idx = overlaps_merged$loopcounts_idx,
               segmentCN_left = assay( segmentCN_sample[overlaps_merged$left_segmentCN_idx,] ) %>% as.vector(), 
               segmentCN_right = assay( segmentCN_sample[overlaps_merged$right_segmentCN_idx,] ) %>% as.vector()) 
    segmentCN_anchors <-  segmentCN_anchors %>%
        mutate(segmentCN_left = segmentCN_left * 2,
               segmentCN_right = segmentCN_right * 2) %>% 
        #mutate(segmentCN_left = (2^segmentCN_left)*2,
               #segmentCN_right = (2^segmentCN_right)*2) %>%
        group_by(loopcounts_idx) %>%
        dplyr::summarize(segmentCN_left = mean(segmentCN_left),
                         segmentCN_right = mean(segmentCN_right))
    segmentCN_anchors$product <- segmentCN_anchors$segmentCN_left * segmentCN_anchors$segmentCN_right

    segmentCN_products <- segmentCN_anchors %>% dplyr::select(-c(segmentCN_left, segmentCN_right)) 
    colnames(segmentCN_products)[colnames(segmentCN_products) == 'product'] <- sample_name
    return(segmentCN_products)
})
loopcountsCN <- Reduce(function(x,y){merge(x,y,by='loopcounts_idx')}, loopcountsCN) %>% data.frame()
loopcountsCN <- loopcountsCN[ , c('loopcounts_idx', colnames(loopcounts_iset))]

timecheck_2 <- Sys.time()
timecheck_2 - timecheck_1

head(loopcountsCN)

#divide loopcount by (CN product + 1)
stopifnot( colnames(loopcounts_iset) == colnames(loopcountsCN[,2:ncol(loopcountsCN)]) )

loopcounts_iset_cn_regressed <- loopcounts_iset[loopcountsCN$loopcounts_idx, ]
names(assays(loopcounts_iset_cn_regressed)) <- 'counts'
assays(loopcounts_iset_cn_regressed)$cn_products <- as.matrix(loopcountsCN[,2:ncol(loopcountsCN)])
assays(loopcounts_iset_cn_regressed)$counts_cn_regressed <-
    assays(loopcounts_iset_cn_regressed)$counts / ( assays(loopcounts_iset_cn_regressed)$cn_products + 1 )
assays(loopcounts_iset_cn_regressed) <- assays(loopcounts_iset_cn_regressed)['counts_cn_regressed']
loopcounts_iset_cn_regressed

MAEobject3 <- c(MAEobject2, 
                enhancerCN = enhancerCN_summarized,
                hichip_cn_regressed = hichip_cn_regressed_summarized,
                hichip_gene_peaks_iset_cn_regressed = hichip_gene_peaks_iset_cn_regressed,
                normcounts_fithichip_loops = normcounts_fithichip_loops_iset,
                normcounts_fithichip_loops_cn_regressed = loopcounts_iset_cn_regressed)
MAEobject3
sampleMap(MAEobject3)

saveRDS(MAEobject3, 'intermediates/TCGA_MAEobject3.rds')

#QC CN regressed loop counts

plotCancerTypeUMAP <- function(data){
    t_data <- t(data)
    data_umap <- umap(t_data)
    
    to_plot <- data.frame(data_umap$layout)
    to_plot$project_ID <- Submitters2ProjectIDs(submitter_IDs = rownames(to_plot), sampleid = sampleid)
    mycolors <- tcga_discrete_palettes[['cancer_types']]
    names(mycolors) <- paste('TCGA', names(mycolors), sep = '-')
    p <- ggplot(to_plot , aes(x = X1, y = X2, fill = project_ID) ) +
    geom_point(shape = 21, size = 4, color = '#444444', stroke = 0.1) +
    theme_classic() +
    scale_fill_manual(values = mycolors)
    return(p)
}

iset_regressed <- experiments(MAEobject3)[['normcounts_fithichip_loops_cn_regressed']] 
data_regressed <- assays(iset_regressed)[['counts_cn_regressed']]

iset_not_regressed <- experiments(MAEobject3)[['normcounts_fithichip_loops']] 
data_not_regressed <- assays(iset_not_regressed)[[1]]

p1 <- plotCancerTypeUMAP(data_not_regressed) + ggtitle('loop counts')
p2 <- plotCancerTypeUMAP(data_regressed) + ggtitle('CN-regressed loop counts')

pdf('figures/UMAPs_cn_regressed_loopcounts2.pdf', useDingbats = FALSE, height = 4.5, width = 10)
grid.arrange(p1, p2, ncol = 2)
dev.off()
