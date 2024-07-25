getPromotersOfGene <- function(gene, upstream = 100, downstream = 100, genome = 'hg38'){
    require(org.Hs.eg.db)
    
    if (genome == 'hg38'){
        require(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    } else if (genome == 'hg19'){
        require(TxDb.Hsapiens.UCSC.hg19.knownGene)
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    } else { stop('genome does not exist') }
        
    promoters <- suppressWarnings( promoters(txdb, upstream = upstream, downstream = downstream) )
    promoters <- GenomeInfoDb::keepStandardChromosomes(promoters, pruning.mode = "coarse")
    eid <- suppressMessages( AnnotationDbi::select(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL') )
    transcript_id <- suppressMessages( AnnotationDbi::select(txdb, eid$ENTREZID, "TXNAME", "GENEID") )
    gene_key <- merge(eid, transcript_id, by.x = 'ENTREZID', by.y = 'GENEID')
    promoters_of_gene <- promoters %>% subset(tx_name %in% gene_key$TXNAME)
    promoters_of_gene$symbol <- gene_key[match(promoters_of_gene$tx_name, gene_key$TXNAME),'SYMBOL']
    return(promoters_of_gene)   
}

dumpJuicer <- function(hic_files, chr, resolution, juicer_outdir, force_overwrite=FALSE){
    stopifnot( dir.exists(juicer_outdir) )

    resolution_char <- format(resolution, scientific=FALSE)
    chr_ncbi <- chr %>% gsub('chr','',.)
    dump_files <- c()
    chr_resolution_run <- paste(chr, chr, 'res', resolution_char, sep = '_')
    logfile <- paste0(juicer_outdir, chr_resolution_run, '.log')
    cat('', file = logfile, sep = '\n')
    
    timecheck_1 <- Sys.time()
    for (filename in names(hic_files)){
        dump_file <- paste(filename, chr_resolution_run, 'H3K27ac.juicer.dump', sep='_')
        cat(paste0('Dumping juice into: ', dump_file), file = logfile, sep = '\n', append = TRUE)

        if ( file.exists( paste0(juicer_outdir, dump_file) ) & !force_overwrite ) {
            cat('file already exists; set force_overwrite = TRUE to overwrite file', file = logfile, sep = '\n', append = TRUE)
        } else{
            juicer_dump_command <- paste('java -jar /storage/kdriest/tools/juicer/scripts/common/juicer_tools.jar',
                                          'dump observed NONE',
                                          hic_files[filename],
                                          chr_ncbi,
                                          chr_ncbi,
                                          'BP',
                                          resolution_char,
                                          paste0(juicer_outdir, dump_file),
                                          sep = ' ')
            system(juicer_dump_command)
        }
        dump_files[filename] <- list.files(juicer_outdir, pattern = dump_file, full.names = TRUE)
    }    
    timecheck_2 <- Sys.time()
    mins_spent <- difftime(timecheck_2, timecheck_1, units = 'min')
    cat(paste0('Time spent juicing: ', mins_spent, ' mins'), file = logfile, sep = '\n', append = TRUE)    

    return(dump_files)
}

readDumpFiles <- function(dump_files, sample_dictionary){
    sparse_matrix <- sapply(names(dump_files), function(filename){
        sample_name <- as.character( sample_dictionary[filename] )
        df <- suppressMessages(readr::read_delim(dump_files[filename],delim = "\t", col_names = FALSE))
        return(df)
    }, simplify=FALSE, USE.NAMES=TRUE)
    names(sparse_matrix) <- sample_dictionary[ names(sparse_matrix) ]
    return(sparse_matrix)
}

getGIFromSparseMatrix <- function(sparse_matrix, chr, resolution){
    left_bin <- GRanges(seqnames = chr, IRanges(start = sparse_matrix[,1,drop=TRUE] + 1, 
                                                end = (sparse_matrix[,1,drop=TRUE] + resolution)))
    right_bin <- GRanges(seqnames = chr, IRanges(start = sparse_matrix[,2,drop=TRUE] + 1, 
                                                 end = (sparse_matrix[,2,drop=TRUE] + resolution)))
    gi <- GInteractions(left_bin, right_bin)
    gi$counts <- sparse_matrix[,3,drop=TRUE]
    return(gi)    
}

normalizeGICounts <- function(gi, scale_to = 1e+6){
    gi$counts <- gi$counts / sum( gi$counts )
    gi$counts <- gi$counts * scale_to
    return(gi)
}      

getGIsFromHiCFiles <- function(hic_files, 
                               juicer_outdir, 
                               sample_dictionary, 
                               chr, 
                               resolution=10000, 
                               normalize_counts=TRUE,
                               force_overwrite=FALSE){
    dump_files <- dumpJuicer(hic_files = hic_files, 
                             chr = chr, 
                             resolution = resolution, 
                             juicer_outdir = juicer_outdir, 
                             force_overwrite = force_overwrite)   
    sparse_matrix_list <- readDumpFiles(dump_files = dump_files, sample_dictionary = sample_dictionary)
    gi_list <- lapply(sparse_matrix_list, function(sparse_matrix){
        gi <- getGIFromSparseMatrix(sparse_matrix = sparse_matrix, chr = chr, resolution = resolution)
        if (normalize_counts){ gi <- normalizeGICounts(gi = gi) }
        return(gi)
    })
    return(gi_list)
}

calculateOverlapWithBins <- function(bins, target_regions){
    hits <- findOverlaps(bins, target_regions)
    overlaps <- pintersect(bins[queryHits(hits)], target_regions[subjectHits(hits)])
    overlapping_bins <- bins[queryHits(hits)]
    overlapping_bins$percent_overlap <- width(overlaps) / width(target_regions[subjectHits(hits)])
    overlapping_bins <- unique(overlapping_bins)
    return(overlapping_bins)
}

getGIWithAnchor <- function(gi, 
                            anchor_gene = NULL, 
                            chr = NULL, 
                            anchor = NULL, 
                            bp_around_anchor = 1e+6, 
                            window_start = NULL, 
                            window_end = NULL){
    stopifnot( length(anchor_gene)<=1 )
    require(InteractionSet)
    require(GenomicInteractions)
    
    if ( is.null(chr) | is.null(anchor) ){ 
        if ( is.null(anchor_gene) ){ 
            stop('Error: either anchor_gene or chr and anchor coordinate must be provided') 
        } 
        target_regions <- getPromotersOfGene(gene=anchor_gene) %>% reduce()
        chr <- seqnames(target_regions) %>% unique() %>% as.character()
    } else { 
        target_regions <- GRanges(chr,anchor) 
    }

    overlapping_bins <- calculateOverlapWithBins(bins = regions(gi), target_regions = target_regions)
    anchor_gr <- overlapping_bins[ which.max(overlapping_bins$percent_overlap) ]    #which.max?
    anchor_center <- round( (start(anchor_gr) + end(anchor_gr)) / 2 )
    if ( is.null(window_start) | is.null(window_end) ){
        window_start <- anchor_center - bp_around_anchor
        window_end <- anchor_center + bp_around_anchor
        window <- GRanges(chr, IRanges(window_start, window_end))
    }
    
    gi_with_anchor <- viewPoint(gi, anchor_gr, window) %>% unique()
    return(gi_with_anchor)
}

plotVirtual4C <- function(gi_with_anchor, bp_around_anchor = 1e+6){
    
    to_plot <- data.frame( gi_with_anchor )
    to_plot$midbin2 <- (to_plot$start2 + to_plot$end2) / 2
    vantage_point <- unique( (to_plot$start1 + to_plot$end1) / 2 )

    p <- ggplot(to_plot, aes(x=midbin2, y=counts)) + 
    geom_line() + 
    geom_vline(xintercept = vantage_point, color = '#999999', linetype = 'longdash') +
    labs(x=element_blank(), y='normalized EIS') +
    scale_x_continuous(expand = c(0,0), limits = c(vantage_point - bp_around_anchor, vantage_point + bp_around_anchor)) +
    theme_classic() 
    return(p)
}

df2Iset <- function(df, leftbin_cols = 1:5, rightbin_cols = 6:10, counts_cols = 11:ncol(df)){
    leftbin <- df[,leftbin_cols] %>% makeGRangesFromDataFrameInColumnOrder()
    rightbin <- df[,rightbin_cols] %>% makeGRangesFromDataFrameInColumnOrder()
    gi <- GInteractions(leftbin, rightbin)
    counts <- df[,counts_cols] %>% as.matrix()
    colData <- DataFrame(sample_name=colnames(counts))
    iset <- InteractionSet(assays = list(counts=counts), gi, colData = colData)
    return(iset)    
}

##########
#unused
subsetIset <- function(iset, chr, window_start, window_end, anchor_gr){
    window_gr <- GRanges(chr, IRanges(start = window_start, end = window_end))
    window_gi <- GInteractions(window_gr, window_gr)
    iset_subset <- subsetByOverlaps(iset, window_gi, ignore.strand = TRUE, minoverlap = 2)
    iset_with_anchor <- subsetByOverlaps(iset_subset, anchor_gr, ignore.strand=TRUE)
    return(iset_with_anchor)
}
###########