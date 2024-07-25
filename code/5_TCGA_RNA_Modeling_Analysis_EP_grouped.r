library(dplyr)
library(VCA)
library(grid)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(ggpubr)
library(sjstats)
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

source('TCGA_King_RNA_modeling_functions.r')
source('~/jupyter/resources/King_R_functions.r')
source('~/jupyter/resources/King_custom_colors.r')

dirs_to_exclude <- NULL

dirs <- list.dirs('intermediates', full.names = TRUE, recursive = FALSE)
dirs <- dirs[grep('RNA_modeling.*_WGS',dirs)]
names(dirs) <- dirs %>% gsub('.*/RNA_modeling_','',.)
dirs <- dirs %>% 
subset(! names(.) %in% dirs_to_exclude)
dirs

oncogenes_subset <- read.delim('intermediates/oncogenes_subset.txt')
head(oncogenes_subset)
oncogene_names <- oncogenes_subset %>%
    subset(grepl('oncogene', Type)) %>% 
    .$Gene %>% 
    unique() %>% 
    as.character()
head(oncogene_names)

oncogene_names
oncogene_names %>% length

cumvar_pca <- lapply(dirs, function(dir){
    readRDS(paste0(dir, '/cumvar_pca.rds'))
})

modeled_genes <- lapply(dirs, function(dir){
    readRDS(paste0(dir, '/modeled_genes.rds'))
})

percent_var <- lapply(dirs, function(dir){
    readRDS(paste0(dir, '/percent_var_grouped.rds'))
})

model <- lapply(dirs, function(dir){
    readRDS(paste0(dir, '/model.rds'))
})

dat_PC <- lapply(dirs, function(dir){
    readRDS(paste0(dir, '/dat_PC.rds'))
})

modeled_genes %>% lapply(function(x)lapply(x,length))

cumvar_pca_all <- lapply(cumvar_pca, function(cumvar_pca_model){
    cumvar_pca_model %>% 
    do.call(cbind,.) %>%
    data.frame() %>%
    mutate(principal_components=rownames(.)) %>%
    gather('gene','cum_var', 1:(ncol(.)-1)) 
})
cumvar_pca_all %>% lapply(head)

names(cumvar_pca_all) %>% subset(grepl('ploidy_corrected', .)) %>% subset(!grepl('^cn_only', .))

modeling_runs = names(cumvar_pca_all) %>% subset(grepl('ploidy_corrected', .)) %>% subset(!grepl('^cn_only', .))
for (modeling_run in modeling_runs){
    p <- ggplot(cumvar_pca_all[[modeling_run]], aes(x=principal_components,y=cum_var)) + 
    geom_boxplot(width = 0.5, fill='#477CBF') +
    ylim(0,1) +
    theme_classic()
    pdf(paste0('figures/gex_model_all_enhancer_PC_cum_var_', modeling_run, '.pdf'), height=3,width=3,useDingbats = FALSE)
    grid.draw(p)
    dev.off()
}

dat_PC %>% lapply(head)

genes_of_interest <- oncogene_names
for (modeling_run in names(dat_PC)){
    plist <- sapply(genes_of_interest, function(gene){
        if ('PC1' %in% colnames(dat_PC[[modeling_run]][[gene]])) {
            dat_PC[[modeling_run]][[gene]] %>% 
            ggplot(aes(x=PC1, y=PC2, color=RNA)) + 
            geom_point(size=3) + 
            scale_color_viridis(option = 'magma', direction = -1) +
            theme_classic() +
            ggtitle(gene) +
            labs(x='H3K27ac HiChIP PC1', y='H3K27ac HiChIP PC2', color='log2(RNA+1)')   
        } else { NULL }
    }, simplify=FALSE, USE.NAMES=TRUE)

    path = paste0('figures/PCA_enhancer_RNA_plots_examples_', modeling_run)
    dir.create(path = path)
    lapply(names(plist), function(gene){

        pdf(paste0(path,'/PCA_enhancer_RNA_',gene,'.pdf'), 
            width=4, height=3, useDingbats = FALSE)
        grid.draw(plist[[gene]])
        dev.off()
    })    
    
}
#PCA plots (PC1 vs PC2), color by RNA expression, color by CN, color by project.ID

#correlation between variables
mycolors1 <- getColorRampGradient(dark_color = 'firebrick', mid_color = 'white', color_dark_portion = 5/6)
mycolors2 <- getColorRampGradient()
genes_of_interest <- c(oncogene_names) #user-defined

for (modeling_run in names(dat_PC)){
    dir <- paste0('figures/variable_correlation_matrices_', modeling_run, '/')
    dir.create(path = dir)
    save_variable_correlation_heatmaps_for_genes(genes_of_interest = genes_of_interest, 
                                                 dat_PC_list = dat_PC[[modeling_run]],
                                                 matrix1_colors = mycolors1, 
                                                 matrix2_colors = mycolors2, 
                                                 dir = dir)
}

names(modeled_genes)

percent_var_summary <- sapply(names(modeled_genes), function(modeling_run){
    summarize_percent_var(rna_cn_iset_genes = modeled_genes[[modeling_run]][['rna_cn_iset_genes']],
                          rna_cn_only_genes = modeled_genes[[modeling_run]][['rna_cn_only_genes']],
                          percent_var_list = percent_var[[modeling_run]],
                          model_list = model[[modeling_run]])
}, simplify=FALSE, USE.NAMES=TRUE)
percent_var_summary[['loopcounts_cn_regressed_WGS']] %>% head()

percent_var_summary_ranked <- sapply(names(modeled_genes), function(modeling_run){
    get_percent_var_summary_ranked(percent_var_summary[[modeling_run]], 
                                   model[[modeling_run]])
}, simplify=FALSE, USE.NAMES=TRUE)
percent_var_summary_ranked[['loopcounts_cn_regressed_WGS']] %>% head()

MAEobject3 <- readRDS('intermediates/TCGA_MAEobject3.rds')
MAEobject3

rna_counts <- assays(MAEobject3)$TPM_RNAseq 
rna_counts <- log2(rna_counts + 1)
rna_summary <- data.frame(gene=rownames(rna_counts),
                          mean_rna=apply(rna_counts, 1, mean),
                          var_rna=apply(rna_counts, 1, var))
head(rna_summary)

percent_var_CN_EP <- sapply(names(modeled_genes), function(modeling_run){
    percent_var_CN_EP <- percent_var_summary_ranked[[modeling_run]] %>% 
    get_percent_var_CN_EP()
    percent_var_CN_EP <- add_p_percent_var_CN_EP(percent_var_summary_ranked[[modeling_run]],
                                                 percent_var_CN_EP)
    percent_var_CN_EP <- merge(percent_var_CN_EP, rna_summary, by='gene') %>% dplyr::arrange(desc(CN))
    return(percent_var_CN_EP)
}, simplify=FALSE, USE.NAMES=TRUE)
percent_var_CN_EP[['loopcounts_cn_regressed_WGS']] %>% head()

#bookmark, run following next

var_rna_cutoff <- 1
percent_var_CN_EP_subset <- lapply(percent_var_CN_EP, function(df){
    df %>% 
    subset(var_rna > var_rna_cutoff) %>% 
    mutate(gene_order=seq_len(nrow(.)),
           gene=factor(gene, levels = gene))
})
percent_var_CN_EP_subset[['loopcounts_cn_regressed_WGS']] %>% head()

#plot CN ranked gene plot for model with only CN as input vs full model: 
#how much does EP component impact CN component?
to_plot <- merge(percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']],
                 percent_var_CN_EP_subset[['cn_only_WGS_ploidy_corrected']][,c('gene','CN')],
                 by='gene') %>%
    dplyr::rename(CN_with_EP=CN.x,CN_without_EP=CN.y)
head(to_plot)

pdf('figures/ranked_gene_plot_CN_withorwithout_EP_WGS_ploidy_corrected.pdf', width=7, height=4)
ranked_gene_plot(data=to_plot, x='gene_order', y_above='CN_with_EP', y_below='CN_without_EP') +
scale_x_continuous(expand = c(0, 0), ) +
labs(y='fraction RNA variance', fill='explained by')
dev.off()

fdr_cutoff <- 0.05
percent_var_CN_EP_subset <- lapply(percent_var_CN_EP_subset, function(df){
    df %>%
    mutate(CN_p_adj = p.adjust(CN_p_value, method = 'fdr'),
           EP_p_adj = p.adjust(EP_p_value, method = 'fdr')) %>%
    mutate(CN_significance = ifelse(CN_p_adj<fdr_cutoff, 'CN_driven', 'NS'),
           EP_significance = ifelse(EP_p_adj<fdr_cutoff, 'EP_driven', 'NS'))
})
percent_var_CN_EP_subset %>% lapply(head)

saveRDS(percent_var_CN_EP_subset, 'intermediates/percent_var_CN_EP_subset.rds')

#ranked gene plot for percent variance explained by CN and EP 
for (modeling_run in c('cn_regressed_WGS_ploidy_corrected','loopcounts_cn_regressed_WGS_ploidy_corrected',
                       'loopcounts_not_cn_regressed_WGS_ploidy_corrected','not_cn_regressed_WGS_ploidy_corrected')){
    plist <- ranked_gene_plot_with_insets(percent_var_CN_EP_subset[[modeling_run]], 
                                          genes_of_interest=genes_of_interest, 
                                          n_top_genes=500)

    pdf(paste0('figures/ranked_gene_plot_CN_EP_',modeling_run,'.pdf'), width=7, height=4)
    print(plist[['p']])
#    print(plist[['CN_inset']],vp=plist[['CN_vp']])
#    print(plist[['EP_inset']],vp=plist[['EP_vp']])
    dev.off()    
}

#ranked gene plot for percent variance explained by CN and EP with fdr
to_plot <- percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']] 

mycolors <- c('CN_driven' = my_custom_palettes[['crookback_bog']][1],
              'EP_driven' = my_custom_palettes[['crookback_bog']][5],
              'NS' = '#99999933')

pdf(paste0('figures/ranked_gene_plot_CN_EP_fdr_cn_regressed_WGS_ploidy_corrected.pdf'), width=7, height=4)
ggplot(to_plot) + 
geom_bar(aes(x=gene_order, y=CN, fill=CN_significance), stat="identity", position="identity") +
geom_bar(aes(x=gene_order, y=-1*E_P, fill = EP_significance), stat="identity", position="identity") +
scale_fill_manual(values = mycolors) +
scale_y_continuous(labels = abs) + 
theme_classic() +
hide_labels_theme(hide_labels=TRUE, flip=FALSE) + 
scale_x_continuous(expand = c(0, 0)) +
labs(y='fraction RNA variance', fill = element_blank())
dev.off()

write.table(percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']], 
            'intermediates/percent_var_CN_EP_subset_cn_regressed_WGS_ploidy_corrected.txt', 
            quote = FALSE, sep = '\t', row.names = FALSE)

percent_var_CN_EP_subset[['loopcounts_cn_regressed_WGS_ploidy_corrected']] %>% 
subset(CN_significance == 'CN_driven') %>% nrow

#ranked gene plot for percent variance explained by CN and EP with fdr
to_plot <- percent_var_CN_EP_subset[['loopcounts_cn_regressed_WGS_ploidy_corrected']] 

mycolors <- c('CN_driven' = my_custom_palettes[['crookback_bog']][1],
              'EP_driven' = my_custom_palettes[['crookback_bog']][5],
              'NS' = '#99999933')

pdf(paste0('figures/ranked_gene_plot_CN_EP_fdr_loopcounts_cn_regressed_WGS_ploidy_corrected.pdf'), width=7, height=4)
ggplot(to_plot) + 
geom_bar(aes(x=gene_order, y=CN, fill=CN_significance), stat="identity", position="identity") +
geom_bar(aes(x=gene_order, y=-1*E_P, fill = EP_significance), stat="identity", position="identity") +
scale_fill_manual(values = mycolors) +
scale_y_continuous(labels = abs) + 
theme_classic() +
hide_labels_theme(hide_labels=TRUE, flip=FALSE) + 
scale_x_continuous(expand = c(0, 0)) +
labs(y='fraction RNA variance', fill = element_blank())
dev.off()

#p value heatmap for ranked gene plot
for (modeling_run in names(percent_var_CN_EP_subset)){
    to_plot <- percent_var_CN_EP_subset[[modeling_run]]
    dat <- to_plot[,c('CN_p_adj','EP_p_adj')] %>% as.matrix()
    dat <- t(dat)
    dat <- -log10( dat )
    dat <- scale_by_row(dat)
    mycolors <- getColorRampGradient(color_dark_portion = 0.9)
    plot <- Heatmap(dat, name='-log10(fdr) scaled', show_column_names = FALSE,
                    col = mycolors, cluster_rows = FALSE, cluster_columns = FALSE)
    
    pdf(paste0('figures/ranked_gene_plot_CN_EP_p_values_',modeling_run,'_WGS.pdf'), width=10, height=1)
    draw(plot)
    dev.off()    
}

df1 <- percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']] %>% mutate(enhancer_data = '1D_signal')
df2 <- percent_var_CN_EP_subset[['loopcounts_cn_regressed_WGS_ploidy_corrected']] %>% 
mutate(enhancer_data = 'loopcounts')

to_plot <- rbind(df1, df2) %>%
mutate(significance = ifelse(CN_significance == 'CN_driven' & EP_significance == 'EP_driven', 'CN_EP', 
                             ifelse(CN_significance == 'CN_driven', 'CN_driven',
                                    ifelse(EP_significance == 'EP_driven', 'EP_driven', 'NS'))))
head(to_plot)

freq_summary <- to_plot %>% 
group_by(enhancer_data, significance) %>%
dplyr::summarize(Freq = n())
freq_summary

mypalette <- my_custom_palettes[['crookback_bog']]
mycolors <- c('CN_driven' = mypalette[1], 
              'CN_EP' = mypalette[3], 
              'EP_driven' = mypalette[5], 
              'NS' = '#99999966')

pdf('figures/CN_EP_bar_chart_WGS_ploidy_corrected.pdf', width = 3.5, height = 2.5)
ggplot(freq_summary, aes(x=enhancer_data, y=Freq, fill=significance))+
geom_bar(stat = "identity") +
scale_fill_manual(values = mycolors) +
theme_classic() +
labs(y = 'number of genes', fill = element_blank()) 
dev.off()

for (modeling_run in names(percent_var_CN_EP_subset)){
    to_plot <- percent_var_CN_EP_subset[[modeling_run]]
    dat <- to_plot[,c('mean_rna','var_rna')] %>% as.matrix()
    dat <- t(dat)
    dat <- scale_by_row(dat)
    mycolors <- getColorRampGradient(color_dark_portion = 0.9)
    plot <- Heatmap(dat, name='scaled value', show_column_names = FALSE,
                    col = mycolors, cluster_rows = FALSE, cluster_columns = FALSE)
    
    pdf(paste0('figures/ranked_gene_plot_mean_var_rna_',modeling_run,'.pdf'), width=10, height=1)
    draw(plot)
    dev.off()    
}

percent_var_CN_EP_subset %>% names

mycolors <- c('CN_driven' = my_custom_palettes[['crookback_bog']][1],
              'EP_driven' = my_custom_palettes[['crookback_bog']][5],
              'NS' = '#99999933')

to_plot <- percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']] %>%
    subset(gene %in% oncogene_names) %>%
    dplyr::arrange(desc(CN)) %>%
    subset(CN_significance != 'NS') %>%
    head(50)
to_plot$gene <- factor(to_plot$gene, levels = to_plot$gene)
to_plot

pdf('figures/top_CN_genes_barplot_ploidy_corrected.pdf', width=4, height=3)
ggplot(to_plot) + 
geom_bar(aes(x=gene, y=CN, fill = CN_significance), stat="identity", position="identity") +
geom_bar(aes(x=gene, y=-1*E_P, fill = EP_significance), stat="identity", position="identity") +
scale_fill_manual(values = mycolors) +
scale_y_continuous(labels = abs, limits = c(-0.8,0.8)) + 
theme_classic() +
coord_flip() +
scale_x_discrete(limits = rev(levels(to_plot[,'gene']))) +
labs(y='fraction RNA variance', x = element_blank(), fill = element_blank()) +
ggtitle('CN-driven cancer driver genes')
dev.off()

to_plot <- percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']] %>%
    subset(gene %in% oncogene_names) %>%
    dplyr::arrange(desc(E_P)) %>%
    subset(EP_significance != 'NS') %>%
    head(50)
to_plot$gene <- factor(to_plot$gene, levels = to_plot$gene)
head(to_plot)

pdf('figures/top_EP_genes_barplot_ploidy_corrected.pdf', width=4, height=7)
ggplot(to_plot) + 
geom_bar(aes(x=gene, y=CN, fill = CN_significance), stat="identity", position="identity") +
geom_bar(aes(x=gene, y=-1*E_P, fill = EP_significance), stat="identity", position="identity") +
scale_fill_manual(values = mycolors) +
scale_y_continuous(labels = abs, limits = c(-0.8,0.8)) + 
theme_classic() +
coord_flip() +
scale_x_discrete(limits = rev(levels(to_plot[,'gene']))) +
labs(y='fraction RNA variance', x = element_blank(), fill = element_blank()) +
ggtitle('EP-driven cancer driver genes')
dev.off()

n_enhancer_inputs <- rowData(experiments(MAEobject3)$hichip_gene_peaks_iset_cn_regressed) %>% 
data.frame() %>% 
group_by(anchor1.Gene) %>%
dplyr::summarize(n_enhancer_inputs = length(unique(as.character(anchor2.PeakID)))) 
head(n_enhancer_inputs)

to_plot <- percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']] %>%
    subset(gene %in% oncogene_names) %>%
    merge(n_enhancer_inputs, by.x = 'gene', by.y = 'anchor1.Gene', all.x = TRUE)
head(to_plot)

mycolors <- rev(getColorRampGradient(dark_color = 'navy'))

pdf('figures/enhancer_input_RNA_variance_oncogenes_ploidy_corrected.pdf', width = 6, height = 5, useDingbats = FALSE)
ggplot(to_plot, aes(x = n_enhancer_inputs, y = -log10(EP_p_adj), fill = -log10(CN_p_adj), label = gene)) +
geom_point(data = to_plot %>% subset(CN_significance == 'NS'), size = 4, shape = 21) +
geom_point(data = to_plot %>% subset(CN_significance != 'NS'), size = 4, shape = 21) +
geom_text_repel() +
theme_classic() +
scale_fill_gradient(high = 'firebrick', low = '#999999') +
labs(x = 'number of enhancer inputs', 
     y = 'enhancer-RNA covariation \n [-log10(adjusted p)]', 
     fill = 'CN-RNA covariation \n [-log10(adjusted p)]')
dev.off()

to_plot <- percent_var_CN_EP_subset[['loopcounts_cn_regressed_WGS_ploidy_corrected']] %>%
    subset(gene %in% oncogene_names) %>%
    merge(n_enhancer_inputs, by.x = 'gene', by.y = 'anchor1.Gene', all.x = TRUE)
head(to_plot)

mycolors <- rev(getColorRampGradient(dark_color = 'navy'))

pdf('figures/enhancer_input_RNA_variance_loopcounts_oncogenes_ploidy_corrected.pdf', width = 6, height = 5, useDingbats = FALSE)
ggplot(to_plot, aes(x = n_enhancer_inputs, y = -log10(EP_p_adj), fill = -log10(CN_p_adj), label = gene)) +
geom_point(data = to_plot %>% subset(CN_significance == 'NS'), size = 4, shape = 21) +
geom_point(data = to_plot %>% subset(CN_significance != 'NS'), size = 4, shape = 21) +
geom_text_repel() +
theme_classic() +
scale_fill_gradient(high = 'firebrick', low = '#999999') +
labs(x = 'number of enhancer inputs', 
     y = 'enhancer-RNA covariation \n [-log10(adjusted p)]', 
     fill = 'CN-RNA covariation \n [-log10(adjusted p)]')
dev.off()

to_plot <- percent_var_CN_EP_subset[['cn_regressed_WGS_ploidy_corrected']] %>%
    merge(n_enhancer_inputs, by.x = 'gene', by.y = 'anchor1.Gene', all.x = TRUE)
head(to_plot)

mycolors <- rev(getColorRampGradient(dark_color = 'navy'))

pdf('figures/enhancer_input_RNA_variance_all_genes_ploidy_corrected.pdf', width = 6, height = 5, useDingbats = FALSE)
ggplot(to_plot, aes(x = n_enhancer_inputs, y = -log10(EP_p_adj), fill = -log10(CN_p_adj), label = gene)) +
geom_point(data = to_plot %>% subset(CN_significance == 'NS'), size = 4, shape = 21) +
geom_point(data = to_plot %>% subset(CN_significance != 'NS'), size = 4, shape = 21) +
geom_text_repel() +
theme_classic() +
scale_fill_gradient(high = 'firebrick', low = '#999999') +
labs(x = 'number of enhancer inputs', 
     y = 'enhancer-RNA covariation \n [-log10(adjusted p)]', 
     fill = 'CN-RNA covariation \n [-log10(adjusted p)]')
dev.off()

percent_var_EP_sum <- sapply(names(percent_var_summary), function(model_run){
    df <- percent_var_summary[[model_run]] %>%
    group_by(gene, category) %>%
    dplyr::summarize(percent_var = sum(percent_var)) %>% 
    spread(category, percent_var) %>%
    dplyr::select(-CN) %>%
    mutate(enhancer_data = model_run)
    if( ! 'E_P' %in% colnames(df) ){ df$E_P <- 0 }
    df$E_P[is.na(df$E_P)] <- 0
    return(df)
}, simplify = FALSE, USE.NAMES = TRUE)
percent_var_EP_sum %>% lapply(head, 20)

to_plot <- rbind(percent_var_EP_sum[['cn_regressed_WGS_ploidy_corrected']], 
                 percent_var_EP_sum[['loopcounts_cn_regressed_WGS_ploidy_corrected']]) %>%
    mutate(enhancer_data = sub('loopcounts_cn_regressed_WGS_ploidy_corrected', 'loopcounts', enhancer_data)) %>% 
    mutate(enhancer_data = sub('cn_regressed_WGS_ploidy_corrected', '1D_signal', enhancer_data))
head(to_plot)

pdf('figures/1D_vs_loopcounts_EP_ploidy_corrected.pdf', useDingbats = FALSE, width = 3.5, height = 3)
ggplot(to_plot, aes(x = enhancer_data, y = E_P)) +
geom_point(color = '#333333') +
geom_violin(trim = FALSE, scale = 'width', fill = '#999999') +
geom_boxplot(width = 0.1, fill = '#666666', color = 'black') +
theme_classic() +
stat_compare_means()
dev.off()

p_value_cutoff <- 0.05
percent_var_CN_EP_sig <- sapply(names(modeled_genes), function(modeling_run){
    percent_var_summary_ranked[[modeling_run]] %>% 
    group_by(gene, category) %>%
    dplyr::summarize( percent_var=sum(percent_var[ coeff_p_value < p_value_cutoff ]) ) %>%
    get_percent_var_CN_EP()
}, simplify=FALSE, USE.NAMES=TRUE)
percent_var_CN_EP_sig[['cn_regressed_WGS_ploidy_corrected']] %>% head()

percent_var_boxplots <- sapply(c('cn_regressed_WGS_ploidy_corrected', 'not_cn_regressed_WGS_ploidy_corrected'), 
                               function(modeling_run){
    dat_percent_var_categories <- get_dat_percent_var_categories(percent_var_summary_ranked[[modeling_run]])
    plot <- plot_percent_var_boxplot(dat_percent_var_categories) + ggtitle(modeling_run)
    return(plot)
}, simplify=FALSE, USE.NAMES=TRUE)

percent_var_significant_boxplots <- sapply(c('cn_regressed_WGS_ploidy_corrected', 
                                             'not_cn_regressed_WGS_ploidy_corrected'), 
                                           function(modeling_run){
    dat_percent_var_categories <- get_dat_percent_var_categories(percent_var_summary_ranked[[modeling_run]],
                                                                coeff_p_value_cutoff=0.05)
    plot <- plot_percent_var_boxplot(dat_percent_var_categories) + ggtitle(modeling_run)
    return(plot)
}, simplify=FALSE, USE.NAMES=TRUE)

pdf('figures/RNA_modeling_percent_var_boxplots_ploidy_corrected.pdf', width=6, height=5)
grid.arrange(grobs = percent_var_boxplots, ncol = 1)
dev.off()

pdf('figures/RNA_modeling_percent_var_significant_boxplots_ploidy_corrected.pdf', width=6, height=5)
grid.arrange(grobs = percent_var_significant_boxplots, ncol = 1)
dev.off()

#get strongest predictor component from E-P interactions
get_max_E_P_predictor_component <- function(data_summary, input_gene_name){
    max_E_P_predictor_component <- data_summary %>% 
        subset(gene==input_gene_name) %>% 
        subset(category=='E_P') %>%
        subset(coeff_p_value<p_value_cutoff) %>%
        subset(percent_var==max(.$percent_var)) %>%
        .$type
    return(max_E_P_predictor_component)
}

#plot scatter function
plot_top_predictors_scatter <- function(data_summary, dat_PC_df, input_gene_name){
    max_E_P_predictor_component <- get_max_E_P_predictor_component(data_summary = data_summary, input_gene_name)
    to_plot <- dat_PC_df[[input_gene_name]]
    
    #inverse PC signs if slope is negative
    if( cor(x = to_plot[,max_E_P_predictor_component], y = to_plot[,'RNA']) < 0){ 
        to_plot[,max_E_P_predictor_component] <- -1 * to_plot[,max_E_P_predictor_component]
    }
    
    to_plot %>%
    ggplot(aes_string(x=max_E_P_predictor_component, y='RNA', color='CN')) + 
    geom_point() +
    scale_color_viridis(option = 'magma', direction = -1) +
    theme_classic() +
    labs(color='CN (log scaled)', x='top E-P component (log scaled)', y='RNA (log scaled)') +
    ggtitle(input_gene_name)
}

plist <- lapply(oncogene_names, function(gene){ 
    tryCatch({
        p <- plot_top_predictors_scatter(percent_var_summary_ranked[['cn_regressed_WGS_ploidy_corrected']], 
                                         dat_PC[['cn_regressed_WGS_ploidy_corrected']], 
                                         gene) 
    }, error = function(cond){ p <<- NULL })

    pdf(paste0('figures/scatter_EP_CN_RNA_ploidy_corrected/scatter_EP_CN_RNA_',gene,'.pdf'), 
        height = 3, width = 4, useDingbats = FALSE)
    grid.draw(p)
    dev.off()
    
    return(p)
})







#extract enhancer PC loadings
#top loadings vs enhancer-gex correlation
input_gene <- 'UBR5'
sig_EP_components <- percent_var_summary_ranked %>% 
subset(gene==input_gene) %>% 
subset(coeff_p_value<p_value_cutoff) %>% 
subset(category=='E_P') %>%
dplyr::arrange(desc(percent_var))
sig_EP_components

enhancer_loadings <- pca_result[[input_gene]]$rotation %>% 
abs() %>% 
data.frame()
enhancer_loadings %>% head(20)

peak_coords <- data.frame(rowRanges(experiments(MAEobject2_matchedsamples)$normcounts_hichip))

enhancer_loadings_gex <- RNA_H3K27ac_ChIP_correlations %>% subset(Gene==input_gene) %>%
merge(enhancer_loadings, by.x='PeakID',by.y='row.names') %>%
merge(peak_coords, ., by = 'PeakID')
head(enhancer_loadings_gex)

ggplot(enhancer_loadings_gex, aes_string(x=sig_EP_components$type[1], y='correlation')) + 
geom_point()

#increase boot runs to 1000
#genes with statistically sign. increased %variance by CN (boot real vs null). same for E-P
#plot %explained variance by hichip vs RNA_var - variable genes explained by E-P?
#delete percent_var_grouped and percent_var_grouped_null

#k27ac 1D pca all peaks vs selected peaks from modeling -> better clustering?

#go back to scatters, genes that are strongly correlated with enhancer signals vs those that correlate with CN

#model improvements:
#increase stringency for enhancer selection for calculating mean enhancer signal - correlation cutoff
#better processing of RNA and CN data? PCA regression for each cancer type??
#repeat modeling with all genes
#add ATAC (score from gex-correlated elements) to gex predicting model
#incorporate looping (count matrix!)
