#calculate scaled, mean-centered values for a matrix by row
scale_by_row <- function(matrix){
    matrix_scaled <- t(matrix) %>% apply(2, function(column){ scale(column) }) 
    matrix_scaled <- t(matrix_scaled)
    colnames(matrix_scaled) <- colnames(matrix)
    rownames(matrix_scaled) <- rownames(matrix)
    return(matrix_scaled)
}

#return mean of top samples per row. Number of top samples is user-defined
#rows: elements; columns: samples
returnMeanTopPerRow <- function(matrix, n_top_samples=3){
    apply(matrix, 1, function(row){
        row_sorted <- row[order(row, decreasing = TRUE)]
        mean_top_samples <- mean( head(row_sorted, n_top_samples) )
        return(mean_top_samples)
    })
}

#correlate two matrices by row. The two matrices must have identical dimensions with identical row orders and column orders
correlateMatricesByRow <- function(matrix1, matrix2){
    stopifnot( dim(matrix1)==dim(matrix2) )

    cor_test <- lapply(seq_len(nrow(matrix1)), function(row){
        skip_boolean <- FALSE
        tryCatch({

        cor_test <- cor.test(matrix1[row,],
                             matrix2[row,],
                             use = "complete.obs")

        out <- c('correlation'=cor_test$estimate,
                 'p_value'=cor_test$p.value)

        }, error = function(cond) {skip_boolean <<- TRUE} )
        if (skip_boolean) { 
            out <- c('correlation'=NA, 'p_value'=NA)
        }
        return(out)        
    })
    
    if ( is.null(rownames(matrix1)) ) { 
        matrix1_rownames <- rep(NA, nrow(matrix1)) 
    } else { matrix1_rownames <- rownames(matrix1) }
    
    if ( is.null(rownames(matrix2)) ) { 
        matrix2_rownames <- rep(NA, nrow(matrix2)) 
    } else { matrix2_rownames <- rownames(matrix2) }

    df <- data.frame(matrix1=matrix1_rownames,
                     matrix2=matrix2_rownames,
                     correlation=sapply(cor_test, function(i){i['correlation.cor']}),
                     p_value=sapply(cor_test, function(i){i['p_value']}))
    return(df)
}

#get motif matches from a specified column of a motif matrix. 
#Each column is a TF and each row is a peak/element, TRUE=motif match, FALSE=no match
#takes matrix from output of matchMotifs() from motifmatchr
#motif_subset_matrix is a matrix with subsetted rows representing the sample pool of interest for enrichment analysis
#hypergeometric test of enrichment using phyper
calculateMotifEnrichment <- function(motif_matrix, motif_subset_matrix, motif_column){
    motifmatches_in_subset <- sum( motif_subset_matrix[,motif_column] )
    n_peaks_in_subset <- nrow(motif_subset_matrix)
    motifmatches_in_peakset <- sum( motif_matrix[,motif_column] )
    n_peaks_in_peakset <- nrow(motif_matrix)

    motif_enrichment_p_value <- phyper(motifmatches_in_subset-1, 
                                       motifmatches_in_peakset, 
                                       n_peaks_in_peakset-motifmatches_in_peakset, 
                                       n_peaks_in_subset,
                                       lower.tail = FALSE)
    
    fraction_matched_subset <- motifmatches_in_subset / n_peaks_in_subset
    fraction_matched_peakset <- motifmatches_in_peakset / n_peaks_in_peakset
   
    out <- c(motifmatches_in_subset=motifmatches_in_subset,
             motifmatches_in_peakset=motifmatches_in_peakset,
             n_peaks_in_peakset=n_peaks_in_peakset,
             n_peaks_in_subset=n_peaks_in_subset,
             fraction_matched_subset=fraction_matched_subset,
             fraction_matched_peakset=fraction_matched_peakset,
             motif_enrichment_p_value=motif_enrichment_p_value)
    return(out)
}

returnMotifEnrichmentDF <- function(motif_enrichment_list){
    df <- data.frame(
        TF=names(motif_enrichment_list),
        motifmatches_in_subset=sapply(motif_enrichment_list, function(i){i['motifmatches_in_subset']}),
        motifmatches_in_peakset=sapply(motif_enrichment_list, function(i){i['motifmatches_in_peakset']}),
        n_peaks_in_peakset=sapply(motif_enrichment_list, function(i){i['n_peaks_in_peakset']}),
        n_peaks_in_subset=sapply(motif_enrichment_list, function(i){i['n_peaks_in_subset']}),
        fraction_matched_subset=sapply(motif_enrichment_list, function(i){i['fraction_matched_subset']}),
        fraction_matched_peakset=sapply(motif_enrichment_list, function(i){i['fraction_matched_peakset']}),
        motif_enrichment_p_value=sapply(motif_enrichment_list, function(i){i['motif_enrichment_p_value']}))
    rownames(df) <- NULL
    df <- df %>% dplyr::arrange(motif_enrichment_p_value)  
    return(df)
}