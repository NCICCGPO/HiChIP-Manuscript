library(dplyr)
library(VCA)
library(grid)
library(ggplot2)
library(viridis)
library(ggpubr)
library(lme4)
library(sjstats)
library(lmerTest)
library(mixOmics)
library(tidyr)
library(relaimpo)
library(ggrepel)
library(plyr)
library(MultiAssayExperiment)
library(RaggedExperiment)
library(RColorBrewer)
library(ggExtra)
library(gridExtra)
library(ppcor)
library(ComplexHeatmap)

source('King_R_functions.r')
source('TCGA_King_RNA_modeling_functions.r')

MAEobject3 <- readRDS('intermediates/TCGA_MAEobject3.rds')
MAEobject3

intersect(colnames(assays(MAEobject3)$normcounts_hichip), 
          colnames(assays(MAEobject3)$enhancerCN)) %>% 
intersect(colnames(assays(MAEobject3)$normcounts_RNAseq)) %>% 
length

#add interaction_ID to rowData of iset. add this earlier
rowData(experiments(MAEobject3)$hichip_gene_peaks_iset_cn_regressed) <- 
    rowData(experiments(MAEobject3)$hichip_gene_peaks_iset_cn_regressed) %>%
    data.frame() %>%
    mutate(interaction_ID = paste(anchor2.PeakID,anchor1.Gene,sep='_')) %>%
    DataFrame()

MAEobject3_matchedsamples <- intersectColumns(MAEobject3[,,names(MAEobject3)!='AA_structures'])
MAEobject3_matchedsamples_matchedgenes <- 
    intersectRows(MAEobject3_matchedsamples[,,c('normcounts_RNAseq','TPM_RNAseq','geneCN_WGS_ploidy_corrected')])

#filter out genes with low expression: fewer than 3 samples with >10 TPM
TPM_cutoff <- 10
n_samples_above_cutoff <- 3
expressed_genes <- apply(assays(MAEobject3_matchedsamples)$TPM_RNAseq, 1, 
                         function(x){subset(x, x > TPM_cutoff)}) %>% 
    lapply(length) %>% 
    subset(. > n_samples_above_cutoff) %>% 
    names()
head(expressed_genes)
length(expressed_genes)

dir = 'intermediates/RNA_modeling_cn_regressed_WGS_ploidy_corrected/'
dir.create(dir)
modeled_genes <- 
    model_rna_expression(MAE_matchedsamples=MAEobject3_matchedsamples,
                         MAE_matchedsamples_matchedgenes=MAEobject3_matchedsamples_matchedgenes,
                         iset_name='hichip_gene_peaks_iset_cn_regressed',
                         rna_experiment_name='normcounts_RNAseq',
                         cn_experiment_name='geneCN_WGS_ploidy_corrected',
                         nPC=5,
                         genes_to_retain=expressed_genes,
                         rna_cn_iset_genes=NULL,
                         rna_cn_only_genes=NULL,
                         model_cn_only=FALSE,
                         dir=dir)

#modeling rna expression with non-cn-regressed hichip iset
dir = 'intermediates/RNA_modeling_not_cn_regressed_WGS_ploidy_corrected/'
dir.create(dir)
modeled_genes_not_regressed <- 
    model_rna_expression(MAE_matchedsamples=MAEobject3_matchedsamples,
                         MAE_matchedsamples_matchedgenes=MAEobject3_matchedsamples_matchedgenes,
                         iset_name='hichip_gene_peaks_iset',
                         rna_experiment_name='normcounts_RNAseq',
                         cn_experiment_name='geneCN_WGS_ploidy_corrected',
                         nPC=5,
                         genes_to_retain=expressed_genes,
                         rna_cn_iset_genes=NULL,
                         rna_cn_only_genes=NULL,
                         model_cn_only=FALSE,
                         dir=dir)

#modeling rna expression with cn only
dir = 'intermediates/RNA_modeling_cn_only_WGS_ploidy_corrected/'
dir.create(dir)
modeled_genes_cn_only <- 
    model_rna_expression(MAE_matchedsamples=MAEobject3_matchedsamples,
                         MAE_matchedsamples_matchedgenes=MAEobject3_matchedsamples_matchedgenes,
                         iset_name='hichip_gene_peaks_iset_cn_regressed',
                         rna_experiment_name='normcounts_RNAseq',
                         cn_experiment_name='geneCN_WGS_ploidy_corrected',
                         nPC=5,
                         genes_to_retain=expressed_genes,
                         rna_cn_iset_genes=NULL,
                         rna_cn_only_genes=NULL,
                         model_cn_only=TRUE,
                         dir=dir)

modeled_genes %>% lapply(length)
modeled_genes_not_regressed %>% lapply(length)
modeled_genes_cn_only %>% lapply(length)

model_rna_expression_loopcounts <- function(MAE_matchedsamples, #MAEobject with matched sample columns,
                                 MAE_matchedsamples_matchedgenes, #MAEobject with matched genes between RNA and CN
                                 iset_name, #interactionSet experiment name 
                                 rna_experiment_name, #RNA experiment name
                                 cn_experiment_name, #gene-level CN experiment name
                                 nPC, #number of principal components for E-P correlations
                                 genes_to_retain, #genes to use for the modeling
                                 gene_TSS_gr, #GRanges with gene TSS coordinates. column "Gene" contains gene names
                                 rna_cn_iset_genes=NULL, #or specify using a vector
                                 rna_cn_only_genes=NULL, #or specify using a vector
                                 model_cn_only=FALSE, #if true, use only CN to model RNA expression
                                 dir #directory to save to
                                ){
    stopifnot(dir.exists(dir))
    stopifnot(iset_name %in% names(experiments(MAE_matchedsamples)))
    stopifnot(rna_experiment_name %in% names(experiments(MAE_matchedsamples_matchedgenes)))
    stopifnot(cn_experiment_name %in% names(experiments(MAE_matchedsamples_matchedgenes)))
    set.seed(69)
    require(mixOmics)

    rna_experiment <- experiments(MAE_matchedsamples_matchedgenes[genes_to_retain,,])[[rna_experiment_name]]
    rna_data <- get_rna_data(rna_experiment)
    rna_matrix <- rna_data[['matrix']]
    rna_rowdata <- rna_data[['rowdata']]

    cn_experiment <- experiments(MAE_matchedsamples_matchedgenes[genes_to_retain,,])[[cn_experiment_name]]
    cn_data <- get_cn_data(cn_experiment)
    cn_matrix <- cn_data[['matrix']]
    cn_matrix <- transform_cn_matrix_colnames(cn_matrix=cn_matrix, 
                                              MAE=MAE_matchedsamples_matchedgenes, 
                                              cn_experiment_name=cn_experiment_name)
    cn_rowdata <- cn_data[['rowdata']]

    rna_cn_genes <- rna_rowdata$gene_name %>% subset(. %in% cn_rowdata$gene_name) %>% unique()
    
    if ( model_cn_only ) {
        iset_genes <- NULL
        rna_cn_iset_genes <- NULL
    } else {
        iset <- experiments(MAE_matchedsamples)[[iset_name]] 
        gi <- interactions(iset)
        retained_gene_TSS_gr <- gene_TSS_gr %>% subset(Gene %in% genes_to_retain)
        iset_genes <- subsetByOverlaps(retained_gene_TSS_gr, gi, ignore.strand = TRUE) %>% .$Gene
    } 
    if ( is.null(rna_cn_iset_genes) & is.null(rna_cn_only_genes) ) {
        rna_cn_iset_genes <- rna_cn_genes %>% subset(. %in% iset_genes)
        rna_cn_only_genes <- rna_cn_genes %>% subset(!. %in% iset_genes)
    }
    
    pca_result <- list()
    cumvar_pca <- list()
    dat_PC <- list()
    dat_PC_null <- list()
    percent_var_grouped <- list()
    percent_var_grouped_null <- list()
    percent_var <- list()
    percent_var_null <- list()
    model <- list()
    null_model <- list()
    enhancer_groups <- list()

    timecheck_1 <- Sys.time()
    logfile <- paste0(dir,'rna_cn_iset_genes_log.txt')
    cat('', file=logfile, append=FALSE)

    for (gene in rna_cn_iset_genes){
        skip_boolean <- FALSE
        tryCatch({    
            
        target_gene_TSS_gr <- gene_TSS_gr %>% subset(Gene == gene)
        iset_idx <- findOverlaps(gi, target_gene_TSS_gr, ignore.strand = TRUE) %>% queryHits()
        iset_matrix_subset <- assay( iset[iset_idx,] ) %>% data.frame() %>% t()
        loopIDs <- paste0('loop', iset_idx)
        colnames(iset_matrix_subset) <- loopIDs

        if (length(loopIDs)>nPC){
            #PCA for putative enhancers to reduce dimensionality
            pca_result[[gene]] <- pca(iset_matrix_subset, ncomp = nPC, center = TRUE, scale = TRUE)
            cumvar_pca[[gene]] <- pca_result[[gene]]$cum.var
            enhancer_groups[[gene]] <- paste0('PC',1:nPC)
            hichip_columns <- pca_result[[gene]]$x
            hichip_formula <- paste(enhancer_groups[[gene]], collapse ='+')
        } else if (length(loopIDs)>1) { 
            enhancer_groups[[gene]] <- loopIDs
            hichip_columns <- iset_matrix_subset 
            hichip_formula <- paste(loopIDs, collapse ='+')
        } else {
            enhancer_groups[[gene]] <- NULL
            hichip_columns <- iset_matrix_subset 
            hichip_formula <- paste(loopIDs, collapse ='+')        
        }

        #merge PCs with RNA and CN
        gene_rna_cn_matrix <- get_gene_rna_cn(rna_matrix, cn_matrix, gene)
        dat_PC[[gene]] <- gene_rna_cn_matrix %>% 
        merge(hichip_columns, by='row.names', all.x = TRUE)

        #scramble RNA-1D/CN pairs to create null model
        dat_PC_null[[gene]] <- dat_PC[[gene]]
        randomrows1 <- sample(nrow(dat_PC_null[[gene]]))
        randomrows2 <- sample(nrow(dat_PC_null[[gene]]))
        dat_PC_null[[gene]]$RNA <- dat_PC_null[[gene]]$RNA[randomrows1]
        dat_PC_null[[gene]]$CN <- dat_PC_null[[gene]]$CN[randomrows2]

        #linear regression model
        formula <- as.formula(paste0('RNA ~ CN + ', hichip_formula))
        model[[gene]] <- eval(bquote(lm(formula, data=dat_PC[[gene]])))
        null_model[[gene]] <- eval(bquote(lm(formula, data=dat_PC_null[[gene]])))

        #percent explained variance of each predictor of model; enhancer group vs CN
        percent_var_grouped[[gene]] <- calc.relimp(model[[gene]],type="lmg",rela=FALSE, 
                                                   groups = enhancer_groups[[gene]], groupnames = 'E_P')
        percent_var_grouped_null[[gene]] <- calc.relimp(null_model[[gene]],type="lmg",rela=FALSE, 
                                                        groups = enhancer_groups[[gene]], groupnames = 'E_P')    

        #percent explained variance of each predictor of model
        percent_var[[gene]] <- calc.relimp(model[[gene]],type="lmg",rela=FALSE)
        percent_var_null[[gene]] <- calc.relimp(null_model[[gene]],type="lmg",rela=FALSE)    

        }, error = function(cond) {
            skip_boolean <<- TRUE
            cat(paste(gene, as.character(cond), sep='\n'), file=logfile, sep="\n", append=TRUE)
        })
        if (skip_boolean) { next } 
    }
    timecheck_2 <- Sys.time()
    mins_spent <- difftime(timecheck_2,timecheck_1, units = 'mins')
    cat(paste0('Time spent: ', mins_spent, 'mins'), file=logfile, sep="\n", append=TRUE)

    timecheck_1 <- Sys.time()
    logfile <- paste0(dir,'rna_cn_only_genes_log.txt')
    cat('', file=logfile, append=FALSE)

    for (gene in rna_cn_only_genes){
        skip_boolean <- FALSE
        tryCatch({    

        gene_rna_cn_matrix <- get_gene_rna_cn(rna_matrix, cn_matrix, gene)
        dat_PC[[gene]] <- gene_rna_cn_matrix

        #scramble RNA-1D/CN pairs to create null model
        dat_PC_null[[gene]] <- dat_PC[[gene]]
        randomrows1 <- sample(nrow(dat_PC_null[[gene]]))
        randomrows2 <- sample(nrow(dat_PC_null[[gene]]))
        dat_PC_null[[gene]]$RNA <- dat_PC_null[[gene]]$RNA[randomrows1]
        dat_PC_null[[gene]]$CN <- dat_PC_null[[gene]]$CN[randomrows2]

        #linear regression model
        formula <- as.formula('RNA ~ CN')
        model[[gene]] <- eval(bquote(lm(formula, data=dat_PC[[gene]])))
        null_model[[gene]] <- eval(bquote(lm(formula, data=dat_PC_null[[gene]])))

        }, error = function(cond) {
            skip_boolean <<- TRUE
            cat(paste(gene, as.character(cond), sep='\n'), file=logfile, sep="\n", append=TRUE)
        })
        if (skip_boolean) { next } 
    }
    timecheck_2 <- Sys.time()
    mins_spent <- difftime(timecheck_2,timecheck_1, units = 'mins')
    cat(paste0('Time spent: ', mins_spent, 'mins'), file=logfile, sep="\n", append=TRUE)
    
    modeled_genes <- list(rna_cn_iset_genes=rna_cn_iset_genes, rna_cn_only_genes=rna_cn_only_genes)

    saveRDS(pca_result, paste0(dir, 'pca_result.rds'))
    saveRDS(cumvar_pca, paste0(dir, 'cumvar_pca.rds'))
    saveRDS(dat_PC, paste0(dir, 'dat_PC.rds'))
    saveRDS(dat_PC_null, paste0(dir, 'dat_PC_null.rds'))
    saveRDS(percent_var_grouped, paste0(dir, 'percent_var_grouped.rds'))
    saveRDS(percent_var_grouped_null, paste0(dir, 'percent_var_grouped_null.rds'))
    saveRDS(percent_var, paste0(dir, 'percent_var.rds'))
    saveRDS(percent_var_null, paste0(dir, 'percent_var_null.rds'))
    saveRDS(model, paste0(dir, 'model.rds'))
    saveRDS(null_model, paste0(dir, 'null_model.rds'))
    saveRDS(modeled_genes, paste0(dir, 'modeled_genes.rds'))
    return(modeled_genes)
}

gene_TSS_gr <- readRDS('intermediates/gene_TSS_gr.rds')
gene_TSS_gr

dir = 'intermediates/RNA_modeling_loopcounts_cn_regressed_WGS_ploidy_corrected/'
dir.create(dir)
modeled_genes_loopcounts_cn_regressed <- 
    model_rna_expression_loopcounts(MAE_matchedsamples=MAEobject3_matchedsamples,
                         MAE_matchedsamples_matchedgenes=MAEobject3_matchedsamples_matchedgenes,
                         iset_name='normcounts_fithichip_loops_cn_regressed',
                         rna_experiment_name='normcounts_RNAseq',
                         cn_experiment_name='geneCN_WGS_ploidy_corrected',
                         nPC=5,
                         genes_to_retain=expressed_genes,
                         gene_TSS_gr = gene_TSS_gr,
                         rna_cn_iset_genes=NULL,
                         rna_cn_only_genes=NULL,
                         model_cn_only=FALSE,
                         dir=dir)

dir = 'intermediates/RNA_modeling_loopcounts_not_cn_regressed_WGS_ploidy_corrected/'
dir.create(dir)
modeled_genes_loopcounts_not_cn_regressed <- 
    model_rna_expression_loopcounts(MAE_matchedsamples=MAEobject3_matchedsamples,
                         MAE_matchedsamples_matchedgenes=MAEobject3_matchedsamples_matchedgenes,
                         iset_name='normcounts_fithichip_loops',
                         rna_experiment_name='normcounts_RNAseq',
                         cn_experiment_name='geneCN_WGS_ploidy_corrected',
                         nPC=5,
                         genes_to_retain=expressed_genes,
                         gene_TSS_gr = gene_TSS_gr,
                         rna_cn_iset_genes=NULL,
                         rna_cn_only_genes=NULL,
                         model_cn_only=FALSE,
                         dir=dir)



#model cancer type contribution

assignments <- read.delim('intermediates/HiChIP_sample_assignments.txt')
assignments <- assignments %>% 
mutate(submitter_id = submitter_id %>% gsub('-','\\.', .))
head(assignments)
dim(assignments)

idx <- colData(MAEobject3_matchedsamples)$submitter_ID %>% 
gsub('\\.01.*','',.) %>% 
match(., assignments$submitter_id)

colData(MAEobject3_matchedsamples) <- colData(MAEobject3_matchedsamples) %>% 
data.frame %>% 
mutate(Cohort = assignments[idx , 'Cohort']) %>% 
DataFrame

colData(MAEobject3_matchedsamples)

MAEobject3_matchedsamples_matchedgenes <- 
    intersectRows(MAEobject3_matchedsamples[,,c('normcounts_RNAseq','TPM_RNAseq','geneCN_WGS_ploidy_corrected')])

#multiple linear regression modeling of rna expression;
#genes that only appear in rna and cn matrices but not the iset matrix are modeled using only cn as input;
#rna_cn_iset_genes and rna_cn_only_genes can be user-defined using vectors of gene names
#function returns these two gene vectors as outputs, and saves modeling results as .rds files.
#.rds files can be found in the user-specified directory dir
model_rna_expression_with_cancer_type <- function(
    MAE_matchedsamples, #MAEobject with matched sample columns,
    MAE_matchedsamples_matchedgenes, #MAEobject with matched genes between RNA and CN
    iset_name, #interactionSet experiment name 
    rna_experiment_name, #RNA experiment name
    cn_experiment_name, #gene-level CN experiment name
    nPC, #number of principal components for E-P correlations
    genes_to_retain, #genes to use for the modeling
    rna_cn_iset_genes=NULL, #or specify using a vector
    rna_cn_only_genes=NULL, #or specify using a vector
    model_cn_only=FALSE, #if true, use only CN to model RNA expression
    dir #directory to save to
){
    stopifnot(dir.exists(dir))
    stopifnot(iset_name %in% names(experiments(MAE_matchedsamples)))
    stopifnot(rna_experiment_name %in% names(experiments(MAE_matchedsamples_matchedgenes)))
    stopifnot(cn_experiment_name %in% names(experiments(MAE_matchedsamples_matchedgenes)))
    set.seed(69)
    require(mixOmics)

    rna_experiment <- experiments(MAE_matchedsamples_matchedgenes[genes_to_retain,,])[[rna_experiment_name]]
    rna_data <- get_rna_data(rna_experiment)
    rna_matrix <- rna_data[['matrix']]
    rna_rowdata <- rna_data[['rowdata']]

    cn_experiment <- experiments(MAE_matchedsamples_matchedgenes[genes_to_retain,,])[[cn_experiment_name]]
    cn_data <- get_cn_data(cn_experiment)
    cn_matrix <- cn_data[['matrix']]
    cn_matrix <- transform_cn_matrix_colnames(cn_matrix=cn_matrix, 
                                              MAE=MAE_matchedsamples_matchedgenes, 
                                              cn_experiment_name=cn_experiment_name)
    cn_rowdata <- cn_data[['rowdata']]

    rna_cn_genes <- rna_rowdata$gene_name %>% subset(. %in% cn_rowdata$gene_name) %>% unique()
    
    if ( model_cn_only ) {
        iset_genes <- NULL
        rna_cn_iset_genes <- NULL
    } else {
        iset <- experiments(MAE_matchedsamples)[[iset_name]] 
        iset_data <- get_iset_data(iset, genes_to_retain)
        iset_matrix <- iset_data[['matrix']]
        iset_rowdata <- iset_data[['rowdata']]
        stopifnot( 'anchor2.PeakID' %in% colnames(iset_rowdata) )
        iset_rowdata$anchor2.PeakID <- as.character(iset_rowdata$anchor2.PeakID)
        iset_genes <- unique(iset_rowdata$anchor1.Gene)
    } 
    if ( is.null(rna_cn_iset_genes) & is.null(rna_cn_only_genes) ) {
        rna_cn_iset_genes <- rna_cn_genes %>% subset(. %in% iset_genes)
        rna_cn_only_genes <- rna_cn_genes %>% subset(!. %in% iset_genes)
    }
    
    pca_result <- list()
    cumvar_pca <- list()
    dat_PC <- list()
    dat_PC_null <- list()
    percent_var_grouped <- list()
    percent_var_grouped_null <- list()
    percent_var <- list()
    percent_var_null <- list()
    model <- list()
    null_model <- list()
    enhancer_groups <- list()

    timecheck_1 <- Sys.time()
    logfile <- paste0(dir,'rna_cn_iset_genes_log.txt')
    cat('', file=logfile, append=FALSE)

    for (gene in rna_cn_iset_genes){
        skip_boolean <- FALSE
        tryCatch({    

        iset_idx <- which(iset_rowdata$anchor1.Gene==gene)
        gene_1D_matrix <- iset_matrix[iset_idx, , drop = FALSE] %>% t()
        peaks <- iset_rowdata$anchor2.PeakID[iset_idx]            
        colnames(gene_1D_matrix) <- peaks

        if (length(peaks)>nPC){
            #PCA for putative enhancers to reduce dimensionality
            pca_result[[gene]] <- pca(gene_1D_matrix, ncomp = nPC, center = TRUE, scale = TRUE)
            cumvar_pca[[gene]] <- pca_result[[gene]]$cum.var
            enhancer_groups[[gene]] <- paste0('PC',1:nPC)
            hichip_columns <- pca_result[[gene]]$x
            hichip_formula <- paste(enhancer_groups[[gene]], collapse ='+')
        } else if (length(peaks)>1) { 
            enhancer_groups[[gene]] <- peaks
            hichip_columns <- gene_1D_matrix 
            hichip_formula <- paste(peaks, collapse ='+')
        } else {
            enhancer_groups[[gene]] <- NULL
            hichip_columns <- gene_1D_matrix 
            hichip_formula <- paste(peaks, collapse ='+')        
        }

        #merge PCs with RNA and CN
        gene_rna_cn_matrix <- get_gene_rna_cn(rna_matrix, cn_matrix, gene)
        
        coldata_idx <- match(rownames(gene_rna_cn_matrix), colData(MAE_matchedsamples_matchedgenes)$submitter_ID)
        gene_rna_cn_matrix <- gene_rna_cn_matrix %>% 
            mutate(Cohort = colData(MAE_matchedsamples_matchedgenes)[coldata_idx , 'Cohort'])

        dat_PC[[gene]] <- gene_rna_cn_matrix %>% 
        merge(hichip_columns, by='row.names', all.x = TRUE)

        #scramble RNA-1D/CN pairs to create null model
        dat_PC_null[[gene]] <- dat_PC[[gene]]
        randomrows1 <- sample(nrow(dat_PC_null[[gene]]))
        randomrows2 <- sample(nrow(dat_PC_null[[gene]]))
        dat_PC_null[[gene]]$RNA <- dat_PC_null[[gene]]$RNA[randomrows1]
        dat_PC_null[[gene]]$CN <- dat_PC_null[[gene]]$CN[randomrows2]

        #linear regression model
        formula <- as.formula(paste0('RNA ~ CN + Cohort + ', hichip_formula))
        model[[gene]] <- eval(bquote(lm(formula, data=dat_PC[[gene]])))
        null_model[[gene]] <- eval(bquote(lm(formula, data=dat_PC_null[[gene]])))

        #percent explained variance of each predictor of model; enhancer group vs CN
        percent_var_grouped[[gene]] <- calc.relimp(model[[gene]],type="lmg",rela=FALSE, 
                                                   groups = enhancer_groups[[gene]], groupnames = 'E_P')
        percent_var_grouped_null[[gene]] <- calc.relimp(null_model[[gene]],type="lmg",rela=FALSE, 
                                                        groups = enhancer_groups[[gene]], groupnames = 'E_P')    

        #percent explained variance of each predictor of model
        percent_var[[gene]] <- calc.relimp(model[[gene]],type="lmg",rela=FALSE)
        percent_var_null[[gene]] <- calc.relimp(null_model[[gene]],type="lmg",rela=FALSE)    

        }, error = function(cond) {
            skip_boolean <<- TRUE
            cat(paste(gene, as.character(cond), sep='\n'), file=logfile, sep="\n", append=TRUE)
        })
        if (skip_boolean) { next } 
    }
    timecheck_2 <- Sys.time()
    mins_spent <- difftime(timecheck_2,timecheck_1, units = 'mins')
    cat(paste0('Time spent: ', mins_spent, 'mins'), file=logfile, sep="\n", append=TRUE)

    timecheck_1 <- Sys.time()
    logfile <- paste0(dir,'rna_cn_only_genes_log.txt')
    cat('', file=logfile, append=FALSE)

    for (gene in rna_cn_only_genes){
        skip_boolean <- FALSE
        tryCatch({    

        gene_rna_cn_matrix <- get_gene_rna_cn(rna_matrix, cn_matrix, gene)
        dat_PC[[gene]] <- gene_rna_cn_matrix

        #scramble RNA-1D/CN pairs to create null model
        dat_PC_null[[gene]] <- dat_PC[[gene]]
        randomrows1 <- sample(nrow(dat_PC_null[[gene]]))
        randomrows2 <- sample(nrow(dat_PC_null[[gene]]))
        dat_PC_null[[gene]]$RNA <- dat_PC_null[[gene]]$RNA[randomrows1]
        dat_PC_null[[gene]]$CN <- dat_PC_null[[gene]]$CN[randomrows2]

        #linear regression model
        formula <- as.formula('RNA ~ CN')
        model[[gene]] <- eval(bquote(lm(formula, data=dat_PC[[gene]])))
        null_model[[gene]] <- eval(bquote(lm(formula, data=dat_PC_null[[gene]])))

        }, error = function(cond) {
            skip_boolean <<- TRUE
            cat(paste(gene, as.character(cond), sep='\n'), file=logfile, sep="\n", append=TRUE)
        })
        if (skip_boolean) { next } 
    }
    timecheck_2 <- Sys.time()
    mins_spent <- difftime(timecheck_2,timecheck_1, units = 'mins')
    cat(paste0('Time spent: ', mins_spent, 'mins'), file=logfile, sep="\n", append=TRUE)
    
    modeled_genes <- list(rna_cn_iset_genes=rna_cn_iset_genes, rna_cn_only_genes=rna_cn_only_genes)

    saveRDS(pca_result, paste0(dir, 'pca_result.rds'))
    saveRDS(cumvar_pca, paste0(dir, 'cumvar_pca.rds'))
    saveRDS(dat_PC, paste0(dir, 'dat_PC.rds'))
    saveRDS(dat_PC_null, paste0(dir, 'dat_PC_null.rds'))
    saveRDS(percent_var_grouped, paste0(dir, 'percent_var_grouped.rds'))
    saveRDS(percent_var_grouped_null, paste0(dir, 'percent_var_grouped_null.rds'))
    saveRDS(percent_var, paste0(dir, 'percent_var.rds'))
    saveRDS(percent_var_null, paste0(dir, 'percent_var_null.rds'))
    saveRDS(model, paste0(dir, 'model.rds'))
    saveRDS(null_model, paste0(dir, 'null_model.rds'))
    saveRDS(modeled_genes, paste0(dir, 'modeled_genes.rds'))
    return(modeled_genes)
}

dir = 'intermediates/RNA_modeling_cn_regressed_WGS_ploidy_corrected_with_cancer_type/'
dir.create(dir)
modeled_genes_with_cancer_type <- model_rna_expression_with_cancer_type(
    MAE_matchedsamples=MAEobject3_matchedsamples,
    MAE_matchedsamples_matchedgenes=MAEobject3_matchedsamples_matchedgenes,
    iset_name='hichip_gene_peaks_iset_cn_regressed',
    rna_experiment_name='normcounts_RNAseq',
    cn_experiment_name='geneCN_WGS_ploidy_corrected',
    nPC=5,
    genes_to_retain=expressed_genes,
    rna_cn_iset_genes=NULL,
    rna_cn_only_genes=NULL,
    model_cn_only=FALSE,
    dir=dir)
