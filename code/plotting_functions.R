peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

getArchDF <- function(lp, r = 100){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)
    return(df)
}

H3K27Ac_1D_track <- function(fragment, cellFragsRegion=NULL,region=NULL,tileSize=500,pal=NULL,uniqueGroups = unique(fragment$RG), tss=NULL)
{
	options(scipen = 999)
	#calculate for .regionSumArrows manually
	cellNames <- uniqueGroups
	cellGroups <- uniqueGroups
	cellFragsRegion =  subsetByOverlaps(fragment, region, ignore.strand = TRUE)
	regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize

	#Starts
	ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
	ids <- which(ts > 0)
  
	#Ends
	te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
	ide <- which(te > 0)

	#Match
	matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
	#Sparse Matrix
	mat <- Matrix::sparseMatrix(
	i = c(ts[ids], te[ide]),
	j = c(matchID[ids], matchID[ide]),
	x = rep(1,  length(ids) + length(ide)),
	dims = c(length(regionTiles), length(cellNames))
	)
	colnames(mat) <- cellNames
  
	#mat@x[mat@x > 1] <- 1

	#Create Group Matrix
	groupMat <- matrix(0, nrow = length(regionTiles), ncol = length(uniqueGroups))
	colnames(groupMat) <- uniqueGroups
	uniqueGroups <- uniqueGroups[uniqueGroups %in% unique(cellGroups)]
	for(i in seq_along(uniqueGroups)){
    groupMat[,uniqueGroups[i]] <- Matrix::rowSums(mat[,which(cellGroups == uniqueGroups[i]),drop=FALSE])
	}

	#Plot DF
	df <- data.frame(which(groupMat > 0, arr.ind=TRUE))
	df$y <- groupMat[cbind(df[,1], df[,2])]

	#Minus 1 Tile Size
	dfm1 <- df
	dfm1$row <- dfm1$row - 1
	dfm1$y <- 0

	#Plus 1 Size
	dfp1 <- df
	dfp1$row <- dfp1$row + 1
	dfp1$y <- 0

	#Create plot DF
	df <- rbind(df, dfm1, dfp1)
	df <- df[!duplicated(df[,1:2]),]
	df <- df[df$row > 0,]
	df$x <- regionTiles[df$row]
	df$group <- uniqueGroups[df$col]

	#Add In Ends
	dfs <- data.frame(
    col = seq_along(uniqueGroups), 
    row = 1, 
    y = 0,
    x = start(region),
    group = uniqueGroups
	)

	dfe <- data.frame(
    col = seq_along(uniqueGroups),
    row = length(regionTiles),
    y = 0,
    x = end(region),
    group = uniqueGroups
	)
  
	#Final output
	plotDF <- rbind(df,dfs,dfe)
	plotDF <- df[order(df$group,df$x),]
	plotDF <- df[,c("x", "y", "group")]
	
  
	#Quantify the tss reads
	tmpinf = tss
	info = read.table(tmpinf,sep="\t",quote=NULL)
	info$Sta = info$V2 - 1500
	info$End = info$V3 + 1500
	info$Sta = as.integer(info$Sta)
	info$End = as.integer(info$End)
	TSSregion = peakDF2GRanges(info)
	TSSfragment =  subsetByOverlaps(fragment, TSSregion, ignore.strand = TRUE)
	groupNormFactors <- table(TSSfragment$RG)

	#Scale with Norm Factors
	scaleFactors <- 10^4 / groupNormFactors
	matchGroup <- match(paste0(plotDF$group), names(scaleFactors))
	plotDF$y <- plotDF$y * as.vector(scaleFactors[matchGroup])
	df <- plotDF

	######################################################
	# Plot Track
	######################################################
	ylim = NULL
	useGroups = NULL
	title = ""
	pal = pal
	baseSize = 7
	scTileSize = 0.5
	scCellsMax = 100
	tickWidth = 0.4
	facetbaseSize = 7
	
	df$y <- zoo::rollmean(df$y, k = 3, fill = 0)
	
  
	if(!is.null(ylim)){
    	ylim <- quantile(df$y, ylim)
    	df$y[df$y < ylim[1]] <- ylim[1]
    	df$y[df$y > ylim[2]] <- ylim[2]
	}else{
    	ylim <- c(0, max(df$y)*1.05)
    	df$y[df$y < ylim[1]] <- ylim[1]
    	df$y[df$y > ylim[2]] <- ylim[2]
	}
  
	uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  
	if(!is.null(useGroups)){
    uniqueGroups <- unique(useGroups)
	}

	df$group <- factor(df$group, levels = uniqueGroups)
	title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)

	normMethod <- "ReadInTSS"
	allGroups <- uniqueGroups

	if(is.null(pal)){
    	pal <- suppressWarnings(paletteDiscrete(values = allGroups))
	}
  
	#Plot Track
	p <- ggplot(df, aes_string("x","y", color = "group", fill = "group")) + 
    geom_area(stat = "identity") + 
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    ylab(sprintf("Coverage\n(Norm. 1D H3K27ac Signal Range (%s-%s) by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    theme_ArchR(baseSize = baseSize,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_blank()) +
    guides(fill = "none", colour = "none") + ggtitle(title)

	p
}

bigWig_track <- function(bigwig_gr, region=NULL,tileSize=500,pal=NULL,uniqueGroups = unique(bigwig_gr$RG), tss = NULL)
{
	options(scipen = 999)
	regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
	regionTiles_gr <- GRanges(seqnames = as.character(unique(seqnames(region))), ranges = IRanges(regionTiles, regionTiles + tileSize-1))
	bigwig_region <- subsetByOverlaps(bigwig_gr, region, ignore.strand = TRUE)
	
	df <- data.frame(matrix(nrow = 0, ncol = 0))
	for (i in uniqueGroups) {
		bigwig_region_i <- bigwig_region[bigwig_region$RG == i]
		df_i <- data.frame(bigwig_region_i %>% group_by_overlaps(x = regionTiles_gr) %>% summarise(y = sum(score)))
		df_i$x <- regionTiles[df_i$query]
		df_i$group <- i
		df_i <- df_i[,c("x", "y", "group")]
		df <- rbind(df, df_i)
	}
	
	df <- df[order(df$group,df$x),]
	df <- df[,c("x", "y", "group")]
	
	
	#Quantify the tss reads
	if (!("groupNormFactors" %in% names(bigwig_gr@metadata))) {
		tmpinf = tss
		info = read.table(tmpinf,sep="\t",quote=NULL)
		info$Sta = info$V2 - 1500
		info$End = info$V3 + 1500
		info$Sta = as.integer(info$Sta)
		info$End = as.integer(info$End)
		TSSregion = peakDF2GRanges(info)
		TSSfragment =  suppressWarnings(subsetByOverlaps(bigwig_gr, TSSregion, ignore.strand = TRUE))
		groupNormFactors <- TSSfragment %>% group_by(RG) %>% summarise(ReadsInTSS = sum(score))
		groupNormFactors <- setNames(groupNormFactors$ReadsInTSS, groupNormFactors$RG)
	} else {
		groupNormFactors <- bigwig_gr@metadata$groupNormFactors
	}


	#Scale with Norm Factors
	scaleFactors <- 10^4 / groupNormFactors
	matchGroup <- match(paste0(df$group), names(scaleFactors))
	df$y <- df$y * as.vector(scaleFactors[matchGroup])
	
	######################################################
	# Plot Track
	######################################################
	ylim = NULL
	useGroups = NULL
	title = ""
	pal = pal
	baseSize = 7
	tickWidth = 0.4
	facetbaseSize = 7
	
	df$y <- zoo::rollmean(df$y, k = 3, fill = 0)
	
	if(!is.null(ylim)){
    	ylim <- quantile(df$y, ylim)
    	df$y[df$y < ylim[1]] <- ylim[1]
    	df$y[df$y > ylim[2]] <- ylim[2]
	}else{
    	ylim <- c(0, max(df$y)*1.05)
    	df$y[df$y < ylim[1]] <- ylim[1]
    	df$y[df$y > ylim[2]] <- ylim[2]
	}
  
	uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  
	if(!is.null(useGroups)){
    uniqueGroups <- unique(useGroups)
	}

	df$group <- factor(df$group, levels = uniqueGroups)
	title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)

	allGroups <- uniqueGroups

	if(is.null(pal)){
    	pal <- suppressWarnings(paletteDiscrete(values = allGroups))
	}
	
	#Plot Track
	p <- ggplot(df, aes_string("x","y", color = "group", fill = "group")) + 
    geom_area(stat = "identity") + 
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    ylab(sprintf("Coverage\n(Norm. Signal Range (%s-%s))", round(min(ylim),2), round(max(ylim),2))) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    theme_ArchR(baseSize = baseSize,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_blank()) +
    guides(fill = "none", colour = "none") + ggtitle(title)

	suppressWarnings(p)
}

virtual4C_plot <- function(hic_files, names, norm, bedpe, region = NULL, window_size = NULL, anchor = NULL
                     , anchor_gene = NULL, pal = NULL, color = NULL, prmtrs = NULL, ymax = NULL,res = 100000, points = 50, bed = NULL, java_path = '~/anaconda3/bin/java', juicer_path = '~/tools/juicer/scripts/common/juicer_tools.jar') {
    options(scipen = 999)
    loops <- read.delim(bedpe)
    loops$dist <- abs(rowMeans(loops[,c("x1", "x2")]) - rowMeans(loops[,c("y1", "y2")]))
    if (!grepl("chr", loops$chr1[1])) {
        loops$chr1 <- paste0("chr", loops$chr1)
        loops$chr2 <- paste0("chr", loops$chr2)
    }
	if (!is.null(region)) {
		anchor_chr <- as.character(unique(seqnames(region)))
	}
    if (!is.null(anchor_gene) & is.null(prmtrs)) {
        prmtrs <- getPromotersOfGene(anchor_gene, genome = "hg38")
        anchor_chr <- unique(as.character(seqnames(prmtrs)))
        if (as.character(strand(prmtrs))[1] == "-") {
            anchor <- max(round((end(prmtrs) + start(prmtrs))/2))
        } else {
            anchor <- min(round((end(prmtrs) + start(prmtrs))/2))
        }
    }

    if (!is.null(window_size)) {
        chr <- as.character(seqnames(prmtrs))
        chr <- names(sort(table(chr), decreasing = TRUE))[1]
		region <- GRanges(seqnames = chr, ranges = IRanges(start = anchor - window_size - res*10, end = anchor + window_size))
    }
    
    vp_plot <- data.frame(matrix(nrow = 0, ncol = 4))
	chr <- as.character(unique(seqnames(region)))

    for (i in 1:length(hic_files)) {
		hic_file <- hic_files[i]
        dump_file <- paste0("~/tmp/", basename(hic_file), "_", chr,"_",res, ".dump")
        system(command = paste0(java_path, " -jar ", juicer_path,
                                " dump observed NONE ", hic_file, " ", gsub("chr", "", anchor_chr), ":", 
                                plyr::round_any(anchor, res, f = floor),":",
                                plyr::round_any(anchor, res, f = ceiling),
                                " ", gsub("chr", "", chr), " BP ", res, 
                                " ", dump_file))
        
        chr_norm <- suppressMessages(readr::read_delim(dump_file,delim = "\t", col_names = FALSE))
        anchor_bin <- unique(chr_norm$X1[anchor >= chr_norm$X1 & anchor < chr_norm$X1 + res])
        if( length(anchor_bin) == 0 ) {
            anchor_bin <- unique(chr_norm$X2[anchor >= chr_norm$X2 & anchor < chr_norm$X2 + res])
        }
        loops_plot <- loops[loops$chr1 == anchor_chr & ((loops$x2 - anchor) < res*2 & (anchor - loops$x1) < res*2 |
			(loops$y2 - anchor) < res*2 & (anchor - loops$y1) < res*2),]
		if (min(loops_plot$Q.Value_Bias) ==0 ) {
			loops_plot$Q.Value_Bias[loops_plot$Q.Value_Bias == 0] <- min(loops_plot$Q.Value_Bias[
				loops_plot$Q.Value_Bias != 0], na.rm = TRUE)
		}
        loops_plot <- loops_plot[order(loops_plot$Q.Value_Bias, decreasing = TRUE),]
        
        loops_gr <- GRanges(seqnames = loops_plot$chr1,
                            ranges = IRanges(start = rowMeans(loops_plot[,c("x1", 'x2')])
                                             , end = rowMeans(loops_plot[,c("y1", 'y2')])))
        loops_gr$value <- -log10(loops_plot$Q.Value_Bias)
        loops_list <- SimpleList(Loops = loops_gr)
        
        loopO <- lapply(seq_along(loops_list), function(x){
            subLoops <- subsetByOverlaps(loops_list[[x]]
                                         , region, ignore.strand = TRUE, type = "within") 
            if(length(subLoops)>0){
                dfx <- getArchDF(subLoops)
                dfx$name <- Rle(paste0(names(loops_list)[x]))
                dfx
            }else{
                NULL
            }})
        chr_norm_vp <- chr_norm[(chr_norm$X1 == anchor_bin & 
                                 chr_norm$X2 > start(region) & chr_norm$X2 < end(region)) ,]
        tmp <- chr_norm[(chr_norm$X2 == anchor_bin & 
                         chr_norm$X1 > start(region) & chr_norm$X1 < end(region)) ,]
        colnames(tmp) <- c("X2", "X1", "X3")
        chr_norm_vp <- rbind(chr_norm_vp, tmp)
        chr_norm_vp <- chr_norm_vp[,c("X1", "X2","X3")]
        if (!(start(region) %in% chr_norm_vp$X2)) {
            chr_norm_vp[nrow(chr_norm_vp) + 1,] = list(anchor_bin,start(region), 0)
        }
        
        if (!(end(region) %in% chr_norm_vp$X2)) {
            chr_norm_vp[nrow(chr_norm_vp) + 1,] = list(anchor_bin,end(region), 0)
        }
        chr_norm_vp <- chr_norm_vp[order(chr_norm_vp$X2),]
        chr_norm_vp_plot <- chr_norm_vp
        chr_norm_vp_plot$sample <- names[i]
        anchor_reads <- sum(chr_norm_vp_plot$X3[chr_norm_vp_plot$X2 == anchor_bin])
        chr_norm_vp_plot$X3 <- chr_norm_vp_plot$X3/ (norm[i]#*anchor_reads
                                                    ) *100000
		chr_norm_vp_plot$X3 <- zoo::rollmean(chr_norm_vp_plot$X3, k = 3, fill = 0)
        vp_plot <- rbind(vp_plot, chr_norm_vp_plot)
    }
    vp_plot$sample <- factor(vp_plot$sample, levels = names)
    
    if (is.null(ymax)) {
        ymax <- max(vp_plot$X3[abs(vp_plot$X2 - vp_plot$X1) > res*5])
    }
	
	baseSize = 7
	tickWidth = 0.4
	facetbaseSize = 7
	
    p1 <- ggplot() + 
        geom_line(data = vp_plot, aes(x = X2, y = X3, color = sample)) + theme_classic() + 
        coord_cartesian(ylim=c(0, ymax)
                        , xlim = c(start(region), end(region))) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        ggtitle(label = paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region),"; ", res /1000, " kb resolution")) + 
		theme_ArchR(baseSize = baseSize,
		                baseLineSize = tickWidth,
		                legendPosition = "right",
		                axisTickCm = 0.1) +
		    theme(panel.spacing= unit(0, "lines"),
		          axis.title.x=element_blank(),
		          axis.text.y=element_blank(),
		          axis.ticks.y=element_blank(),
		          axis.text.x=element_blank(),
		          axis.ticks.x=element_blank(),
		          strip.text = element_text(
		            size = facetbaseSize, 
		            color = "black", 
		            margin = margin(0,0.35,0,0.35, "cm")),
		            strip.text.y = element_text(angle = 0),
		          strip.background = element_blank()) +
		    guides(fill = "none", colour = "none") + 
			  ylab("Norm EIS") + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 2))
    if (chr == anchor_chr){
        p1 <- p1 + geom_vline(xintercept = anchor, color = "grey", linetype = "longdash", size = 1)
    }
    if (is.null(pal)) {
        p1 <- p1 + scale_color_manual(values = paletteDiscrete(values = levels(vp_plot$sample)))
    } else {
        p1 <- p1 + scale_color_manual(values = pal)
    }
    p1 <- p1 + facet_grid(sample ~ .) + theme(legend.position = "none",  
                                                     strip.background = element_blank())
    if (is.null(color)) {
        color <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
    } else {
        color <- colorRampPalette(c("white",color))(100)
    }
        if (!is.null(loopO[[1]])) {
            valueMin <- min(loopO[[1]]$value)
            valueMax <- max(loopO[[1]]$value)
    p2 <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line() +
        facet_grid(name ~ .) + theme_classic() +
        ylab("") + xlab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = color, limits = c(valueMin, valueMax)
                              , name = "-log10(Adjusted P-value)") +
							   theme_ArchR(baseSize = baseSize,
							                  baseLineSize = tickWidth,
							                  legendPosition = "right",
							                  axisTickCm = 0.1) +
							      theme(panel.spacing= unit(0, "lines"),
							            axis.title.x=element_blank(),
							            axis.text.y=element_blank(),
							            axis.ticks.y=element_blank(),
							            strip.text = element_text(
							              size = facetbaseSize, 
							              color = "black", 
							              margin = margin(0,0.35,0,0.35, "cm")),
							              strip.text.y = element_text(angle = 0),
							            strip.background = element_blank()) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3)) 
    }
    else {
        df <- data.frame(facet = "Loops", start = 0, end = 0, strand = "*", symbol = "none")
        p2 <- ggplot(data = df, aes(start, end)) + 
            geom_point() +
            facet_grid(facet~.) + theme_classic() +
            scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0))+
    theme_ArchR(baseSize = baseSize,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_blank()) +
    guides(fill = "none", colour = "none")
        }
    
    height <- 1
    if (length(levels(vp_plot$sample)) > 3) {
        height <- min(1*length(levels(vp_plot$sample)), 6) 
    }
    plot_return <- patchwork::wrap_plots(p1,p2, ncol =1, heights = c(height,0.5))
    plot_return
}

plot_genes <- function(region = NULL, gene_list = NULL, window_size = NULL, anchor_gene = NULL, res = 100000, all_tx = FALSE, coding = TRUE){
	options(dplyr.summarise.inform = FALSE)
	if (!is.null(anchor_gene)) {
        prmtrs <- getPromotersOfGene(anchor_gene, genome = "hg38")
        anchor_chr <- unique(as.character(seqnames(prmtrs)))
        if (as.character(strand(prmtrs))[1] == "-") {
            anchor <- max(round((end(prmtrs) + start(prmtrs))/2))
        } else {
            anchor <- min(round((end(prmtrs) + start(prmtrs))/2))
        }
    }

    if (!is.null(window_size)) {
        chr <- as.character(seqnames(prmtrs))
        chr <- names(sort(table(chr), decreasing = TRUE))[1]
		region <- GRanges(seqnames = chr, ranges = IRanges(start = anchor - window_size - res*10, end = anchor + window_size))
    }
	
	chr <- as.character(unique(seqnames(region)))
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	species <- Homo.sapiens
	genes <- suppressMessages(data.frame(transcripts(txdb, columns = c('gene_id', "tx_id", "tx_name", "cds_id"))))
	genes$gene_id <- as.character(genes$gene_id)
	genes <- genes[genes$seqnames == chr,]
	if (coding) {
		genes <- genes[!unlist(lapply(genes$cds_id, FUN = function(x){all(is.na(x))})),]
	}
	genenames <- suppressMessages(AnnotationDbi::select(species, keys = unique(genes$gene_id)
		, columns = 'SYMBOL', keytype = 'GENEID'))
	genes <- merge(genes, genenames, by.x = "gene_id", by.y = "GENEID")
	genes <- genes[!duplicated(genes[,c('start','end', 'SYMBOL')]),]
	if (!all_tx) {
		message("collapsing")
		genes <- genes %>% group_by(gene_id, SYMBOL) %>% summarise(seqnames = seqnames[1], start = min(start), end = max(end)
			, strand = unique(strand))
	} 
	genes <- genes[(genes$start > start(region) & genes$start < end(region)) | 
		(genes$end > start(region) & genes$end < end(region)),]
	genes$start <- ifelse(genes$start < start(region), start(region), genes$start)
	genes$end <- ifelse(genes$end > end(region), end(region), genes$end)
	genes <- genes[!duplicated(genes[,c('start','end', 'SYMBOL')]),]
	genes <- genes[order(genes$gene_id, genes$start, decreasing = TRUE),]
	if (!is.null(gene_list)) {
	        genes <- genes[genes$SYMBOL %in% gene_list & genes$seqnames == chr,]
	}
	if (all_tx) {
		genes$SYMBOL <- paste0(genes$SYMBOL, '\n', genes$tx_name)
	}
	geneinfobed <- genes[,c("seqnames", "start", "end", "SYMBOL", "strand")]
	colnames(geneinfobed) <- c("chrom", "start", "stop", "gene", "strand")
	geneinfobed$score <- "."
	geneinfobed$pos <- seq(1:nrow(geneinfobed))
	geneinfobed$labels <- "Genes"
	
	baseSize = 7
	tickWidth = 0.4
	facetbaseSize = 7
	
	p <- ggplot(data = geneinfobed) + 
		geom_segment(data = geneinfobed, aes(x = start, xend = stop, y = pos, yend = pos
				, color = strand),size = 2) + 
		geom_label(data = geneinfobed, aes(x = start, y = pos, label = gene, color = strand),nudge_y = -0.5) +
		facet_grid(labels ~ .) + theme_classic() +
		scale_color_manual(values = paletteDiscrete(values = c("+","-")), guide = FALSE) + 
		ylab("") + xlab("") + 
		ylim(min(geneinfobed$pos) -1, max(geneinfobed$pos) +1) + 
		scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
	    theme_ArchR(baseSize = baseSize,
	                baseLineSize = tickWidth,
	                legendPosition = "right",
	                axisTickCm = 0.1) +
	    theme(panel.spacing= unit(0, "lines"),
	          axis.title.x=element_blank(),
	          axis.text.y=element_blank(),
	          axis.ticks.y=element_blank(),
	          strip.text = element_text(
	            size = facetbaseSize, 
	            color = "black", 
	            margin = margin(0,0.35,0,0.35, "cm")),
	            strip.text.y = element_text(angle = 0),
	          strip.background = element_blank()) +
	    guides(fill = "none", colour = "none")
	p
}

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