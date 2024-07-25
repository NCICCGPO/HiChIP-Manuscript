oncogenes_subset <- read.delim('intermediates/oncogenes_subset.txt')
oncogene_names <- oncogenes_subset$Gene %>% unique() %>% as.character()
sampleid <- read.delim('intermediates/sampleid.txt')

formatATACColnames <- function(column_names, sampleid){
    sapply(column_names, function(column_name){
        if(grepl('_',column_name)){
            sample_name <- column_name %>% gsub('_T1.*','',.) %>% gsub('_T2.*','',.)
            submitter_ID <- sampleid %>% 
                subset(grepl(sample_name,Library_Name)) %>% 
                .$Tissue_Barcode %>% 
                unique() %>% 
                as.character() %>%
                strsplit(split = '-') %>% 
                unlist() %>% 
                .[1:4] %>% 
                paste(collapse='.')
            return(submitter_ID)
        } else { return(column_name) }
    })
}

LibraryNames2Submitters <- function(library_names, sampleid){
    sapply(library_names, function(library_name){
        if(grepl('-',library_name)){
            library_name <- library_name %>% 
                strsplit(split = '-') %>%
                unlist() %>%
                .[1:min(6, length(.))] %>%
                paste(collapse='_')
            submitter_ID <- sampleid %>% 
                subset(grepl(library_name,Library_Name)) %>% 
                .$Tissue_Barcode %>% 
                unique() %>% 
                as.character() %>%
                strsplit(split = '-') %>% 
                unlist() %>% 
                .[1:4] %>% 
                paste(collapse='.')
            return(submitter_ID)
        } else { return(library_name) }
    })
}

Submitters2ProjectIDs <- function(submitter_IDs, sampleid){
    sapply(as.character(submitter_IDs), function(submitter){
        submitter <- submitter %>% strsplit(split='\\.') %>% unlist() %>% .[1:3] %>% paste(collapse='-')
        project_ID <- sampleid %>% 
            subset(grepl(submitter,submitter_id)) %>% 
            .$cohort %>% 
            as.character() %>% 
            unique() %>% 
            paste0('TCGA-',.)
        return(project_ID)
    })    
}

makeGRangesFromDataFrameInColumnOrder <- function(df){
    makeGRangesFromDataFrame(df,
                             seqnames.field = colnames(df)[1],
                             start.field = colnames(df)[2],
                             end.field = colnames(df)[3])
}