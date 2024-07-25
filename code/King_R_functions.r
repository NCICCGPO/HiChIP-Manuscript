#create quantile color breaks for heatmap
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#replace outlier values
ReplaceOutliers <- function(x, quantCut = c(.05, .95), removeNA = TRUE) {
     quantiles <- quantile(x, probs=quantCut, na.rm = removeNA)
     x[x < quantiles[1]] <- quantiles[1]
     x[x > quantiles[2]] <- quantiles[2]
     return(x)
}

#calculate concentrations of DNA samples from KAPA Quant, given light cycler output text file and sample_loading_order file
calculate_concentrations_KAPAQ <- function(dat, 
                                           sample_loading_order, 
                                           standard_fragment_size=452,
                                           standards_to_exclude=NULL){
    require(dplyr)
    replicate_order <- unique(sample_loading_order$replicate_order) %>% as.character()
    replicate_loading <- unique(sample_loading_order$replicate_order)
    if (replicate_loading=='consecutive'){
        sample_loading_order_repeated <- sample_loading_order[rep(seq_len(nrow(sample_loading_order)),
                                                                  times = sample_loading_order$n_replicates), 
                                                              c('sample_order', 'dilution_factor', 'mean_fragment_size')]
    }
    
    dat_summary <- cbind(dat, sample_loading_order_repeated) %>% 
        group_by(sample_order,dilution_factor) %>%
        dplyr::summarize(mean_Cp=mean(Cp),
                         mean_fragment_size=mean(as.numeric(mean_fragment_size)))

    standard_concentrations <- c(standard_1=20,
                                 standard_2=2,
                                 standard_3=0.2,
                                 standard_4=0.02,
                                 standard_5=0.002,
                                 standard_6=0.0002) %>% data.frame(concentration=.)
    standard_concentrations$log2_concentration <- log2(standard_concentrations$concentration)
    standard_concentrations <- merge(standard_concentrations, dat_summary, 
                                     by.x = 'row.names', by.y = 'sample_order')
    standard_concentrations <- standard_concentrations %>% subset(! Row.names %in% standards_to_exclude)
    
    linear_regression <- lm(log2_concentration ~ mean_Cp, data = standard_concentrations)
    intercept <- linear_regression$coefficients['(Intercept)']
    slope <- linear_regression$coefficients['mean_Cp']

    dat_summary$log2_concentration <- dat_summary$mean_Cp*slope + intercept
    dat_summary$concentration_pM <- 2^dat_summary$log2_concentration
    dat_summary$concentration_original_nM <- (dat_summary$concentration_pM*dat_summary$dilution_factor)/1000
    dat_summary2 <- dat_summary %>% 
        group_by(sample_order) %>%
        dplyr::summarize(concentration_original_nM=mean(concentration_original_nM),
                         mean_fragment_size=mean(mean_fragment_size)) %>%
        mutate(concentration_original_nM_size_corrected=concentration_original_nM*standard_fragment_size/mean_fragment_size) %>%
        mutate(concentration_original_nM_size_corrected=format(concentration_original_nM_size_corrected, scientific=FALSE))
    dat_summary2 <- dat_summary2[match(unique(sample_loading_order$sample_order), dat_summary2$sample_order),]

    return(list(all_concentrations=dat_summary, mean_concentrations=dat_summary2))
}

#calculate scaled, mean-centered values for a matrix by row
scale_by_row <- function(matrix){
    matrix_scaled <- t(matrix) %>% apply(2, function(column){ scale(column) }) 
    matrix_scaled <- t(matrix_scaled)
    colnames(matrix_scaled) <- colnames(matrix)
    rownames(matrix_scaled) <- rownames(matrix)
    return(matrix_scaled)
}

#reformat Mergedpeaks_counts
reformat_Mergedpeaks_counts <- function(Mergedpeaks_counts, sample_names){
    colnames(Mergedpeaks_counts) <- c(c("chr","start","end", "PeakID","length","strand"), sample_names)
    rownames(Mergedpeaks_counts) <- Mergedpeaks_counts$PeakID
    return(Mergedpeaks_counts)
}

#make dds DESeq object from countsTable and colData
make_dds <- function(countsTable, colData, formula, n_cores=30){
    require(DESeq2)
    require(BiocParallel)
    register(MulticoreParam(n_cores))
    dds <- DESeqDataSetFromMatrix(countsTable, colData, design=as.formula(formula)) 
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds, parallel = TRUE)
    assays(dds)$normcounts <- counts(dds, normalized = TRUE)
    return(dds)
}

#make split violin plot using ggplot2
require(ggplot2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    # Original function by Jan Gleixner (@jan-glx)
    # Adjustments by Wouter van der Bijl (@Axeman)
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
                             
#join two SummarizedExperiment objects by cbind, using a rowData column provided
joinSEsByRowData <- function(se1, se2, join_by_rowData_col){
    overlap_rows <- rowData(se1)[,join_by_rowData_col] %>% 
        subset(. %in% rowData(se2)[,join_by_rowData_col]) %>% 
        as.character()
    se1_overlap <- se1[match(overlap_rows, rowData(se1)[,join_by_rowData_col]),]
    se2_overlap <- se2[match(overlap_rows, rowData(se2)[,join_by_rowData_col]),]
    rowData <- rowData(se1_overlap)
    rowData(se1_overlap) <- rowData(se2_overlap) <- NULL
    joined_se <- cbind(se1_overlap, se2_overlap)
    rowData(joined_se) <- rowData
    return(joined_se)
}
                             
#data summary function
data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE),
         sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}
                             
#function for picking element(s) from a string separated by a symbol
pickElement <- function(string, select, splitBy = '_'){
    string %>%
    strsplit(split = splitBy) %>%
    lapply(function(i){
        new_string <- i[select] %>% paste(collapse = splitBy)
        return(new_string)
    }) %>%
    unlist()
}
                             
#get sequencing coverage in specified window in bedgraph
#window in GRanges format
getCoverageInWindow <- function(bedgraph_gr, window){
    bedgraph_gr$coverage <- as.numeric(bedgraph_gr$name)
    bedgraph_gr$weight <- bedgraph_gr$coverage * width(bedgraph_gr)
    bedgraph_in_window <- subsetByOverlaps(bedgraph_gr, window, ignore.strand = TRUE)
    coverage_in_window <- sum(bedgraph_in_window$weight) / sum(width(window))
    return(coverage_in_window)
}
                             
# RPKM versus TPM
# 
# RPKM and TPM are both normalized for library size and gene length.
#
# RPKM is not comparable across different samples.
#
# For more details, see: http://blog.nextgenetics.net/?e=51

rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}