#filter interactionSet object using a vector of genes
filter_iset <- function(iset, genes_to_retain){
    iset_rowdata <- data.frame(rowData(iset))
    stopifnot( 'anchor1.Gene' %in% colnames(iset_rowdata) )
    iset_idx <- which(iset_rowdata$anchor1.Gene %in% genes_to_retain)
    filtered_iset <- iset[iset_idx]
    return(filtered_iset)
}

#get matrices and rowdata from iset, rna, and cn experiments
get_iset_data <- function(iset, genes_to_retain){
    iset <- filter_iset(iset=iset, genes_to_retain=genes_to_retain)
    iset_matrix <- assay(iset)
    iset_matrix <- log2( iset_matrix + 1)
    iset_matrix <- scale_by_row( iset_matrix )
    iset_rowdata <- data.frame(rowData(iset))
    return(list(matrix=iset_matrix, rowdata=iset_rowdata))
}

get_rna_data <- function(rna_experiment){
    rna_matrix <- assay(rna_experiment)
    rna_matrix <- log2( rna_matrix + 1 )
    rna_matrix <- scale_by_row( rna_matrix )
    rna_rowdata <- rowData(rna_experiment) %>% data.frame()  
    return(list(matrix=rna_matrix, rowdata=rna_rowdata))
}
    
get_cn_data <- function(cn_experiment){
    cn_matrix <- assay(cn_experiment)
    cn_matrix <- log2( cn_matrix + 1 )
    cn_matrix <- scale_by_row( cn_matrix )
    cn_rowdata <- rowData(cn_experiment) %>% data.frame()
    return(list(matrix=cn_matrix, rowdata=cn_rowdata))
}

transform_cn_matrix_colnames <- function(cn_matrix, MAE, cn_experiment_name){
    colnames(cn_matrix) <- colnames(cn_matrix) %>%
    sapply(function(name){
        sampleMap(MAE) %>%
        subset(assay==cn_experiment_name) %>%
        subset(colname==name) %>%
        .$primary
    })
    return(cn_matrix)
}

#get sample-paired rna and cn values for specified gene from matrices
get_gene_rna_cn <- function(rna_matrix, cn_matrix, gene){
    same_rownames_colnames <- all(rownames(rna_matrix)==rownames(cn_matrix)) & all(colnames(rna_matrix)==colnames(cn_matrix))
    stopifnot(same_rownames_colnames)
    df <- data.frame(RNA=rna_matrix[gene,],
                     CN=cn_matrix[gene,])
    return(df)
}

#multiple linear regression modeling of rna expression;
#genes that only appear in rna and cn matrices but not the iset matrix are modeled using only cn as input;
#rna_cn_iset_genes and rna_cn_only_genes can be user-defined using vectors of gene names
#function returns these two gene vectors as outputs, and saves modeling results as .rds files.
#.rds files can be found in the user-specified directory dir
model_rna_expression <- function(MAE_matchedsamples, #MAEobject with matched sample columns,
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

#get percent variances for each variable component from lmg
extract_lmg <- function(gene, percent_var){
    df <- percent_var@lmg %>% data.frame('percent_var'=.) 
    df$gene <- gene
    df$type <- rownames(df)
    df$category <- ifelse(grepl('^PC|^MergePeakID|^loop',df$type), 'E_P', df$type)
    return(df)
}

#get percent variances for CN from r.squared of model
extract_r_squared_cn_only <- function(gene, model){
    if(is.null(model)) { 
        return(NULL)
    } else {
        df <- summary(model)$r.squared %>% data.frame('percent_var'=.) 
        df$gene <- gene
        rownames(df) <- 'CN'
        df$type <- rownames(df)
        df$category <- 'CN'
        return(df)
    }
}

#combine percent variances to make a summary table
summarize_percent_var <- function(rna_cn_iset_genes, rna_cn_only_genes, percent_var_list, model_list){
    percent_var_lmg <- sapply(rna_cn_iset_genes, function(gene){
        if (gene %in% names(percent_var_list)){
            extract_lmg(gene, percent_var_list[[gene]])
        }
    }, simplify=FALSE, USE.NAMES=TRUE)
    
    percent_var_cn_only <- sapply(rna_cn_only_genes, function(gene){
        if (gene %in% names(model_list)){
            extract_r_squared_cn_only(gene,model_list[[gene]])
        }
    }, simplify=FALSE, USE.NAMES=TRUE)

    percent_var_summary <- c(percent_var_lmg, percent_var_cn_only) %>% do.call(rbind,.) 
    return(percent_var_summary)
}

#rank EP components by percent variances
rank_percent_var_EP <- function(percent_var_summary){
    percent_var_EP_ranked <- percent_var_summary %>%
    subset(category=='E_P') %>%
    group_by(gene) %>%
    dplyr::arrange(desc(percent_var), .by_group = TRUE) 
    
    if(nrow(percent_var_EP_ranked)>0){
        percent_var_EP_ranked <- percent_var_EP_ranked %>%
        dplyr::summarize(percent_var=percent_var,
                         type=type,
                         category=category,
                         ranked_category=paste0('enhancer_component_top_',1:length(gene)))
    } 
    percent_var_EP_ranked <- data.frame(percent_var_EP_ranked)
    return(percent_var_EP_ranked)
}
                                
get_percent_var_CN <- function(percent_var_summary){
    percent_var_CN <- percent_var_summary %>% subset(category=='CN')
    percent_var_CN$ranked_category <- percent_var_CN$category
    rownames(percent_var_CN) <- NULL
    return(percent_var_CN)
}

get_percent_var_p_values <- function(percent_var_summary, model_list){
    lapply(unique(percent_var_summary$gene), function(input_gene){
        df <- data.frame(coeff_p_value=summary(model_list[[input_gene]])$coefficients[,4])
        df$gene <- input_gene
        df$type <- rownames(df)
        return(df)
    }) %>% do.call(rbind,.)  
}

get_percent_var_summary_ranked <- function(percent_var_summary, model_list){
    percent_var_ranked_EP <- rank_percent_var_EP(percent_var_summary)
    percent_var_CN <- get_percent_var_CN(percent_var_summary)
    percent_var_p_values <- get_percent_var_p_values(percent_var_summary, model_list)
    percent_var_summary_ranked <- rbind(percent_var_ranked_EP, percent_var_CN) %>% arrange(gene)
    percent_var_summary_ranked <- merge(percent_var_summary_ranked, percent_var_p_values,
                                        by=c('gene','type'))
    return(percent_var_summary_ranked)
}

get_percent_var_CN_EP <- function(percent_var_summary_ranked){
    percent_var_CN_EP <- percent_var_summary_ranked %>%
    dplyr::select(gene,category,percent_var) %>%
    spread(category,percent_var)
    percent_var_CN_EP[is.na(percent_var_CN_EP)] <- 0
    
    if ('E_P' %in% colnames(percent_var_CN_EP)){
        percent_var_CN_EP <- percent_var_CN_EP %>%
        dplyr::arrange(desc(CN), E_P) 
    } else { 
        percent_var_CN_EP <- percent_var_CN_EP %>%
        dplyr::arrange(desc(CN))
    }
    percent_var_CN_EP <- percent_var_CN_EP %>%
    mutate(gene_order=seq_len(nrow(.))) %>%
    mutate(gene=factor(gene, levels=gene))
    return(percent_var_CN_EP)
}

add_p_percent_var_CN_EP <- function(percent_var_summary_ranked, 
                                    percent_var_CN_EP,
                                    ranked_category_CN='CN', 
                                    ranked_category_EP='enhancer_component_top_1'){
    CN_summary <- percent_var_summary_ranked %>%
    subset(ranked_category==ranked_category_CN) %>%
    dplyr::rename(CN_p_value=coeff_p_value)
    
    EP_summary <- percent_var_summary_ranked %>%
    subset(ranked_category==ranked_category_EP) %>%
    dplyr::rename(EP_p_value=coeff_p_value)

    percent_var_CN_EP <- percent_var_CN_EP %>%
    merge(CN_summary[,c('gene','CN_p_value')], by='gene', all.x=TRUE) %>%
    merge(EP_summary[,c('gene','EP_p_value')], by='gene', all.x=TRUE) %>%
    arrange(gene_order)
    
    percent_var_CN_EP$CN_p_value[is.na(percent_var_CN_EP$CN_p_value)] <- 1
    percent_var_CN_EP$EP_p_value[is.na(percent_var_CN_EP$EP_p_value)] <- 1
        
    return(percent_var_CN_EP)
}

hide_labels_theme <- function(hide_labels=TRUE, flip=FALSE){
    if (!hide_labels){ 
        my_theme <- NULL 
    } else if (!flip){
        my_theme <- theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank()) 
    } else {
        my_theme <- theme(axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())  
    }
    return(my_theme)
}

ranked_gene_plot <- function(data, x, y_above, y_below=NULL, flip=FALSE, hide_gene_names=TRUE){
    data[,y_below] <- -1*data[,y_below]
    if (!flip){
        flip_command <- NULL
        axis_order <- NULL
    } else {
        flip_command <- coord_flip() 
        axis_order <- scale_x_discrete(limits = levels(data[,x]))
    }

    p <- ggplot(data) + 
    geom_bar(aes_string(x=x, y=y_above, fill=shQuote(y_above)), stat="identity", position="identity") +
    scale_fill_brewer(palette='Set2') +
    scale_y_continuous(labels = abs) + 
    flip_command +
    axis_order +
    theme_classic() +
    hide_labels_theme(hide_labels=hide_gene_names, flip=flip)

    if (!is.null(y_below)){
        p <- p + geom_bar(aes_string(x=x, y=y_below, fill = shQuote(y_below)), stat="identity", position="identity")
    }
    return(p)
}

ranked_gene_plot_with_insets <- function(percent_var_CN_EP, genes_of_interest, n_top_genes){
    top_CN_genes <- percent_var_CN_EP %>%
    dplyr::arrange(desc(CN)) %>%
    head(n_top_genes)
    top_EP_genes <- percent_var_CN_EP %>%
    dplyr::arrange(desc(E_P)) %>%
    head(n_top_genes)

    p <- ranked_gene_plot(data=percent_var_CN_EP, x='gene_order', y_above='CN', y_below='E_P') + 
    scale_x_continuous(expand = c(0, 0)) +
    labs(y='fraction RNA variance', fill='explained by')

    CN_inset <- ranked_gene_plot(data=top_CN_genes, x='gene', y_above='CN', y_below='E_P') +
    geom_text_repel(aes(x = gene, y = CN, label = gene), 
                    size = 1.7, 
                    data = top_CN_genes %>% subset(gene %in% genes_of_interest)) +
    labs(y=element_blank()) +
    guides(fill=FALSE) +
    ggtitle('top CN') +
    theme(plot.title = element_text(size = 10))

    CN_vp <- viewport(width = 0.35, height = 0.3, x = 0.4, y = 0.8)

    EP_inset <- ranked_gene_plot(data=top_EP_genes, x='gene', y_above='CN', y_below='E_P') +
    geom_text_repel(aes(x = gene, y = -E_P, label = gene), 
                    size = 1.7, 
                    data = top_EP_genes %>% subset(gene %in% genes_of_interest)) +
    labs(y=element_blank()) +
    guides(fill=FALSE) +
    ggtitle('top E-P') +
    theme(plot.title = element_text(size = 10))

    EP_vp <- viewport(width = 0.35, height = 0.3, x = 0.8, y = 0.8)    
    
    return(list(p=p,CN_inset=CN_inset,CN_vp=CN_vp, EP_inset=EP_inset, EP_vp=EP_vp))
}

get_variable_correlations <- function(dat_variables){
    variables_matrix <- dplyr::select_if(dat_variables, is.numeric)
    variable_correlations <- ppcor::pcor(variables_matrix, method = "pearson")
    return(variable_correlations)
}

get_variable_correlation_heatmaps <- function(dat_variables, matrix1_colors, matrix2_colors){
    variable_correlation_matrices <- get_variable_correlations(dat_variables)
    variable_correlation_abs_matrix <- variable_correlation_matrices$estimate %>% abs()
    variable_correlation_p_values_matrix <- variable_correlation_matrices$p.value
    min_non_zero_p <- min(variable_correlation_p_values_matrix[variable_correlation_p_values_matrix>0])
    p_values_offset_for_diagonal <- min_non_zero_p*0.1
    variable_correlation_p_log_matrix <- 
        -log10( variable_correlation_p_values_matrix + p_values_offset_for_diagonal)

    p1 <- Heatmap(variable_correlation_abs_matrix, 
                 name = 'abs(correlation)',
                 col = matrix1_colors,
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 width = unit(4, 'cm'), height = unit(4, 'cm'))
    p2 <- Heatmap(variable_correlation_p_log_matrix,
                  name='-log10(p values)',
                  col = matrix2_colors,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  width = unit(4, 'cm'), height = unit(4, 'cm'))
    return(p1+p2)    
}

save_variable_correlation_heatmaps_for_genes <- function(genes_of_interest,
                                                         dat_PC_list, 
                                                         matrix1_colors, 
                                                         matrix2_colors,
                                                         dir){
    stopifnot(dir.exists(dir))
    set.seed(2)
    logfile <- paste0(dir,'variable_correlation_heatmaps_log.txt')
    cat('', file=logfile, append=FALSE)
    
    timecheck_1 <- Sys.time()
    for (gene in genes_of_interest){
        skip_boolean <- FALSE
        tryCatch({    
        
        heatmaps <- get_variable_correlation_heatmaps(dat_PC_list[[gene]], matrix1_colors, matrix2_colors)
            
        pdf(paste0(dir,'variable_correlations_',gene,'.pdf'), width=6, height=4)
        draw(heatmaps)
        dev.off()

        }, error = function(cond) {
            skip_boolean <<- TRUE
            cat(paste(gene, as.character(cond), sep='\n'), file=logfile, sep="\n", append=TRUE)
       })
        if (skip_boolean) { next } 
    }
    timecheck_2 <- Sys.time()
    mins_spent <- difftime(timecheck_2,timecheck_1, units = 'mins')
    cat(paste0('Time spent: ', mins_spent, 'mins'), file=logfile, sep="\n", append=TRUE)
}

get_dat_percent_var_categories <- function(percent_var_summary_ranked, coeff_p_value_cutoff=1){
    to_plot <- percent_var_summary_ranked %>%
    subset(coeff_p_value < coeff_p_value_cutoff) %>%
    dplyr::select(gene,percent_var,ranked_category)
    to_plot_total <- to_plot %>%
    group_by(gene) %>%
    dplyr::summarize(percent_var=sum(percent_var), ranked_category='total')
    to_plot2 <- rbind(to_plot, to_plot_total) %>% arrange(gene)
    return(to_plot2)
}

plot_percent_var_boxplot <- function(dat_pervent_var_categories, 
                                     percent_var_column='percent_var',
                                     category_column='ranked_category',
                                     var_comparisons=list(c('CN','enhancer_component_top_1'))){
    ggplot(dat_pervent_var_categories, aes_string(x=category_column, y=percent_var_column)) + 
    geom_boxplot(width = 0.5, aes_string(fill=category_column), coef=1e100) +
    scale_fill_brewer(palette = 'Set2') +
    stat_compare_means(comparisons = var_comparisons, size = 3) + 
    labs(x=element_blank()) +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) 
}