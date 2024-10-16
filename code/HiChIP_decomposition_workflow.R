#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to identify the sample specific decomposition analysis by integrating HiChIP and scATAC data
#Input : (1) sample specific loop annotation (2) sample specific h3k27ac peaks (3) matched scATAC data with cluster annotation
#Output : sample specific loop decomposition
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Sample specific decomposition
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)
library(ChIPpeakAnno)
library(ArchR)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$
#GR convert function
#$$$$$$$$$$$$$$$$$$$$$$$$$$$
#define the function
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

#load patient specific loop annotation
tmpinf1 =  "Patient_loops_gene_annotated_knownGene.Rda"
load(tmpinf1)
sam = unique(final_final_tmp$ID)
raw_final_final_tmp = final_final_tmp

#load matched scATAC data
mydir = "/scATAC/Data/"
files = list.files(mydir)
tag = grep("immune_working",files)
files = files[tag]
label = sapply(files,function(x) paste0(strsplit(x,"_")[[1]][2],"-",strsplit(x,"_")[[1]][3],"-",strsplit(x,"_")[[1]][4],"-",strsplit(x,"_")[[1]][5],"-",strsplit(x,"_")[[1]][6]))
label = as.vector(label)

#identify matched samples
sam = intersect(sam, label)

#load patient specific 1d signal
peak_folder = "/TCGA/H3K27Ac_peaks/"
peak_files = list.files(peak_folder)
peak_nam = sapply(peak_files,function(x) paste0(strsplit(x,"-")[[1]][1],"-",strsplit(x,"-")[[1]][2],"-",strsplit(x,"-")[[1]][3],"-",strsplit(x,"-")[[1]][4],"-",strsplit(x,"-")[[1]][5]))
peak_nam = as.vector(peak_nam)

for(i in 1 : length(sam))
{

	cat("\r",i)
	
	######################################
	#identify the loops with h3k27ac peaks
	######################################
	final_final_tmp = raw_final_final_tmp
	tag = grep(sam[i],final_final_tmp$ID)
	final_final_tmp = final_final_tmp[tag,]
	res1 = final_final_tmp
	row.names(res1) = paste0(res1$anchor_1_seqnames,"_",res1$anchor_1_start,"_",res1$anchor_1_end,"_",res1$anchor_2_seqnames,"_",res1$anchor_2_start,"_",res1$anchor_2_end)

	anchor_1 = res1[,c("anchor_1_seqnames","anchor_1_start","anchor_1_end")]
	anchor_2 = res1[,c("anchor_2_seqnames","anchor_2_start","anchor_2_end")]

	anchor_1 = peakDF2GRanges(anchor_1)
	anchor_2 = peakDF2GRanges(anchor_2)
	
	tag = grep(sam[i],peak_nam)
	tmpinf = paste0(peak_folder,peak_files[tag])
	
	peaks = read.table(tmpinf,sep="\t",quote=NULL)
	peaks = peakDF2GRanges(peaks)

	op_1 = findOverlaps(anchor_1, peaks)
	op_2 = findOverlaps(anchor_2, peaks)

	op_1 = as.data.frame(op_1)
	op_2 = as.data.frame(op_2)

	xx1 = unique(op_1$queryHits)
	xx2 = unique(op_2$queryHits)

	xx = intersect(xx1,xx2)
	res1 = res1[xx,]

	anchor_1 = res1[,c("anchor_1_seqnames","anchor_1_start","anchor_1_end")]
	anchor_2 = res1[,c("anchor_2_seqnames","anchor_2_start","anchor_2_end")]

	anchor_1 = peakDF2GRanges(anchor_1)
	anchor_2 = peakDF2GRanges(anchor_2)

	data = res1[,c("anchor_1_seqnames","anchor_1_start","anchor_1_end","anchor_2_seqnames","anchor_2_start","anchor_2_end")]
	
	######################################
	#load immune and tumor peaks
	######################################
	tag = grep(sam[i],label)
	tmpinf = paste0(mydir,files[tag],"/")
	sc_immune_folder = paste0(tmpinf,"PeakCalls/")
	sc_immune_folder = unique(sc_immune_folder)
	
	
	tumor_files = gsub("immune","tumor",files[tag])
	tmpinf = paste0(mydir,tumor_files,"/")
	sc_tumor_folder = paste0(tmpinf,"PeakCalls/")
	sc_tumor_folder = unique(sc_tumor_folder)
	

	immune_peak_files = list.files(sc_immune_folder)
	tumor_peak_files = list.files(sc_tumor_folder)
	
	tag = grep("rds",immune_peak_files)
	immune_peak_files = immune_peak_files[tag]
	Cluster = gsub("-reproduciblePeaks.gr.rds","",immune_peak_files)
	Cluster = unique(Cluster)
	
	tumor_peaks = gsub("_tumor_working","",tumor_files)
	tumor_peak_files = paste0(tumor_peaks,"-reproduciblePeaks.gr.rds")
	
	######################################
	#identify loops with immune peaks
	######################################
	for(k in 1 : length(Cluster))
	{
		tmp_immune_peaks = paste0(sc_immune_folder,immune_peak_files[k])
		C5_peaks = readRDS(tmp_immune_peaks)

		op_1 = findOverlaps(anchor_1, C5_peaks)
		op_2 = findOverlaps(anchor_2, C5_peaks)

		op_1 = as.data.frame(op_1)
		op_2 = as.data.frame(op_2)

		xx1 = unique(op_1$queryHits)
		xx2 = unique(op_2$queryHits)
	
		xx = intersect(xx1,xx2)
		data[,Cluster[k]] = rep(0,nrow(data))
		data[xx,Cluster[k]] = 1
		
		Union = c(xx1,xx2)
		Union = unique(Union)
		tag_union = which(Union %in% xx)
		Union = Union[-tag_union]
		label_union = paste0(Cluster[k],"_Union")
		data[,label_union] = rep(0, nrow(data))
		data[Union,label_union] = 1
	}
	
	######################################
	#identify loops with tumor peaks
	######################################
	tmp_tumor_peaks = paste0(sc_tumor_folder,tumor_peak_files)
	tumor_peaks = readRDS(tmp_tumor_peaks)

	op_1 = findOverlaps(anchor_1, tumor_peaks)
	op_2 = findOverlaps(anchor_2, tumor_peaks)

	op_1 = as.data.frame(op_1)
	op_2 = as.data.frame(op_2)

	xx1 = unique(op_1$queryHits)
	xx2 = unique(op_2$queryHits)

	xx = intersect(xx1,xx2)
	data$Tumor = rep(0,nrow(data))
	data[xx,"Tumor"] = 1
	
	Union = c(xx1,xx2)
	Union = unique(Union)
	tag_union = which(Union %in% xx)
	Union = Union[-tag_union]
	label_union = "Tumor_Union"
	data[,label_union] = rep(0, nrow(data))
	data[Union,label_union] = 1
	
	chose_label = c("anchor_1_seqnames","anchor_1_start","anchor_1_end","anchor_2_seqnames","anchor_2_start","anchor_2_end")
	tag_chose = which(colnames(data) %in% chose_label)
	data = data[,-tag_chose]
	
	#output result
	myoutf = paste0("/Sample_loop_annotation_detail/",sam[i],"_loops_annotation.txt")
	write.table(data,myoutf,sep="\t",quote=F)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to summarize the patient specific decomposition and discover a pre-set of cell-type specific interactions
#Input : sample specific loop annotation 
#Output : pre-set of cell-type specific interactions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Pre-set of cell-type specific interactions
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

tmpinf1 = "/Sample_loop_annotation_detail/BRCA-14AD76EE-12F9-40B3-8DCD_loops_annotation.txt"
tmpinf2 = "/Sample_loop_annotation_detail/BRCA-7C6A3AE4-E2EA-42B3-B3F1_loops_annotation.txt"
tmpinf3 = "/Sample_loop_annotation_detail/BRCA-8D1E6006-85CB-484A-8B5C_loops_annotation.txt"
tmpinf4 = "/Sample_loop_annotation_detail/BRCA-94AF19F0-1F2A-41EC-8CB6_loops_annotation.txt"
tmpinf5 = "/Sample_loop_annotation_detail/BRCA-C147AAD5-A8F1-41D5-8709_loops_annotation.txt"
tmpinf6 = "/Sample_loop_annotation_detail/BRCA-C9C8D426-A3FD-4455-89A9_loops_annotation.txt"
tmpinf7 = "/Sample_loop_annotation_detail/COAD-0914606C-2CA1-4287-B530_loops_annotation.txt"
tmpinf8 = "/Sample_loop_annotation_detail/GBMx-09C0DCE7-D669-4D28-980D_loops_annotation.txt"
tmpinf9 = "/Sample_loop_annotation_detail/GBMx-6BEE2CB6-9AFD-42A6-9C26_loops_annotation.txt"
tmpinf10 = "/Sample_loop_annotation_detail/GBMx-ED12A6C9-D96E-49C4-B882_loops_annotation.txt"
tmpinf11 = "/Sample_loop_annotation_detail/KIRC-7D6A394E-01EC-4C58-A010_loops_annotation.txt"
tmpinf12 = "/Sample_loop_annotation_detail/LUAD-9B9C5C6D-1755-41CD-9BC1_loops_annotation.txt"
tmpinf13 = "/Sample_loop_annotation_detail/LUAD-A1ABC0B2-D7F8-45B5-B431_loops_annotation.txt"
tmpinf14 = "/Sample_loop_annotation_detail/LUAD-CA6F245B-30E0-48DE-AB15_loops_annotation.txt"
tmpinf15 = "/Sample_loop_annotation_detail/SKCM-211D9CF4-3348-4DCD-8A01_loops_annotation.txt"
tmpinf16 = "/Sample_loop_annotation_detail/SKCM-FDA487D2-5293-4315-9212_loops_annotation.txt"

res1 = read.table(tmpinf1,sep="\t",quote=NULL)
res2 = read.table(tmpinf2,sep="\t",quote=NULL)
res3 = read.table(tmpinf3,sep="\t",quote=NULL)
res4 = read.table(tmpinf4,sep="\t",quote=NULL)
res5 = read.table(tmpinf5,sep="\t",quote=NULL)
res6 = read.table(tmpinf6,sep="\t",quote=NULL)
res7 = read.table(tmpinf7,sep="\t",quote=NULL)
res8 = read.table(tmpinf8,sep="\t",quote=NULL)
res9 = read.table(tmpinf9,sep="\t",quote=NULL)
res10 = read.table(tmpinf10,sep="\t",quote=NULL)
res11 = read.table(tmpinf11,sep="\t",quote=NULL)
res12 = read.table(tmpinf12,sep="\t",quote=NULL)
res13 = read.table(tmpinf13,sep="\t",quote=NULL)
res14 = read.table(tmpinf14,sep="\t",quote=NULL)
res15 = read.table(tmpinf15,sep="\t",quote=NULL)
res16 = read.table(tmpinf16,sep="\t",quote=NULL)

######################################
#BRCA_14AD76EE
######################################
tag1 = res1$C3 == 1 & res1$C5 == 0 & res1$Tumor == 0
tag2 = res1$C3 == 0 & res1$C5 == 1 & res1$Tumor == 0
tag3 = res1$C3 == 0 & res1$C5 == 0 & res1$Tumor == 1
label = apply(res1,1,sum)
tag4 = label >= 2
tag5 = label == 0

C3_loop = row.names(res1)[tag1]
C5_loop = row.names(res1)[tag2]
Tumor_loop = row.names(res1)[tag3]

C3_loop_res1 = C3_loop
C5_loop_res1 = C5_loop
Tumor_loop_res1 = Tumor_loop

#######################################
#BRCA_7C6A3AE4
#######################################
tag1 = res2$C3 == 1 & res2$C5 == 0 & res2$Tumor == 0
tag2 = res2$C3 == 0 & res2$C5 == 1 & res2$Tumor == 0
tag3 = res2$C3 == 0 & res2$C5 == 0 & res2$Tumor == 1
label = apply(res2,1,sum)
tag4 = label >= 2
tag5 = label == 0

C3_loop = row.names(res2)[tag1]
C5_loop = row.names(res2)[tag2]
Tumor_loop = row.names(res2)[tag3]

C3_loop_res2 = C3_loop
C5_loop_res2 = C5_loop
Tumor_loop_res2 = Tumor_loop

######################################
#BRCA_8D1E6006
#######################################
tag1 = res3$C5 == 1 & res3$Tumor == 0
tag2 = res3$C5 == 0 & res3$Tumor == 1

label = apply(res3,1,sum)
tag4 = label >= 2
tag5 = label == 0

C5_loop = row.names(res3)[tag1]
Tumor_loop = row.names(res3)[tag2]

C5_loop_res3 = C5_loop
Tumor_loop_res3 = Tumor_loop

######################################
#BRCA_94AF19F0
######################################
tag1 = res4$C3 == 1 & res4$C5 == 0 & res4$Tumor == 0
tag2 = res4$C3 == 0 & res4$C5 == 1 & res4$Tumor == 0
tag3 = res4$C3 == 0 & res4$C5 == 0 & res4$Tumor == 1
label = apply(res4,1,sum)
tag4 = label >= 2
tag5 = label == 0

C3_loop = row.names(res4)[tag1]
C5_loop = row.names(res4)[tag2]
Tumor_loop = row.names(res4)[tag3]

C3_loop_res4 = C3_loop
C5_loop_res4 = C5_loop
Tumor_loop_res4 = Tumor_loop

######################################
#BRCA_C147AAD5
######################################
tag1 = res5$C4 == 1 & res5$C5 == 0 & res5$Tumor == 0
tag2 = res5$C4 == 0 & res5$C5 == 1 & res5$Tumor == 0
tag3 = res5$C4 == 0 & res5$C5 == 0 & res5$Tumor == 1
label = apply(res5,1,sum)
tag4 = label >= 2
tag5 = label == 0

C4_loop = sum(tag1)
C5_loop = sum(tag2)
Tumor_loop = sum(tag3)

C4_loop = row.names(res5)[tag1]
C5_loop = row.names(res5)[tag2]
Tumor_loop = row.names(res5)[tag3]

C4_loop_res5 = C4_loop
C5_loop_res5 = C5_loop
Tumor_loop_res5 = Tumor_loop

######################################
#BRCA_C9C8D426
######################################
tag1 = res6$C5 == 1 & res6$Tumor == 0
tag2 = res6$C5 == 0 & res6$Tumor == 1

label = apply(res6,1,sum)
tag4 = label >= 2
tag5 = label == 0

C5_loop = row.names(res6)[tag1]
Tumor_loop = row.names(res6)[tag2]

C5_loop_res6 = C5_loop
Tumor_loop_res6 = Tumor_loop

#######################################
#COAD_0914606C
#######################################
tag1 = res7$C1 == 1 & res7$C5 == 0 & res7$Tumor == 0
tag2 = res7$C1 == 0 & res7$C5 == 1 & res7$Tumor == 0
tag3 = res7$C1 == 0 & res7$C5 == 0 & res7$Tumor == 1
label = apply(res7,1,sum)
tag4 = label >= 2
tag5 = label == 0

C1_loop = sum(tag1)
C5_loop = sum(tag2)
Tumor_loop = sum(tag3)

C1_loop = row.names(res7)[tag1]
C5_loop = row.names(res7)[tag2]
Tumor_loop = row.names(res7)[tag3]

C1_loop_res7 = C1_loop
C5_loop_res7 = C5_loop
Tumor_loop_res7 = Tumor_loop


#######################################
#GBMx_09C0DCE7
#######################################
tag1 = res8$C6 == 1 & res8$Tumor == 0
tag2 = res8$C6 == 0 & res8$Tumor == 1

label = apply(res8,1,sum)
tag4 = label >= 2
tag5 = label == 0

C6_loop = row.names(res8)[tag1]
Tumor_loop = row.names(res8)[tag2]

C6_loop_res8 = C6_loop
Tumor_loop_res8 = Tumor_loop


#######################################
#GBMx_6BEE2CB6
#######################################
tag1 = res9$C6 == 1 & res9$Tumor == 0
tag2 = res9$C6 == 0 & res9$Tumor == 1

label = apply(res9,1,sum)
tag4 = label >= 2
tag5 = label == 0

C6_loop = row.names(res9)[tag1]
Tumor_loop = row.names(res9)[tag2]

C6_loop_res9 = C6_loop
Tumor_loop_res9 = Tumor_loop

#######################################
#GBMx_ED12A6C9
#######################################
tag1 = res10$C6 == 1 & res10$C5 == 0 & res10$Tumor == 0
tag2 = res10$C6 == 0 & res10$C5 == 1 & res10$Tumor == 0
tag3 = res10$C6 == 0 & res10$C5 == 0 & res10$Tumor == 1
label = apply(res10,1,sum)
tag4 = label >= 2
tag5 = label == 0

C6_loop = row.names(res10)[tag1]
C5_loop = row.names(res10)[tag2]
Tumor_loop = row.names(res10)[tag3]

C6_loop_res10 = C6_loop
C5_loop_res10 = C5_loop
Tumor_loop_res10 = Tumor_loop

#######################################
#KIRC_7D6A394E
#######################################
tag1 = res11$C5 == 1 & res11$Tumor == 0
tag2 = res11$C5 == 0 & res11$Tumor == 1

label = apply(res11,1,sum)
tag4 = label >= 2
tag5 = label == 0

C5_loop = row.names(res11)[tag1]
Tumor_loop = row.names(res11)[tag2]

C5_loop_res11 = C5_loop
Tumor_loop_res11 = Tumor_loop

#######################################
#LUAD_9B9C5C6D
#######################################
tag1 = res12$C1 == 1 & res12$C5 == 0 & res12$Tumor == 0
tag2 = res12$C1 == 0 & res12$C5 == 1 & res12$Tumor == 0
tag3 = res12$C1 == 0 & res12$C5 == 0 & res12$Tumor == 1
label = apply(res12,1,sum)
tag4 = label >= 2
tag5 = label == 0

C1_loop = row.names(res12)[tag1]
C5_loop = row.names(res12)[tag2]
Tumor_loop = row.names(res12)[tag3]

C1_loop_res12 = C1_loop
C5_loop_res12 = C5_loop
Tumor_loop_res12 = Tumor_loop

#######################################
#LUAD_A1ABC0B2
#######################################
tag1 = res13$C1 == 1 & res13$Tumor == 0
tag2 = res13$C1 == 0 & res13$Tumor == 1

label = apply(res13,1,sum)
tag4 = label >= 2
tag5 = label == 0

C1_loop = row.names(res13)[tag1]
Tumor_loop = row.names(res13)[tag2]

C1_loop_res13 = C1_loop
Tumor_loop_res13 = Tumor_loop

#######################################
#LUAD_CA6F245B
#######################################
tag1 = res14$C3 == 1 & res14$C5 == 0 & res14$Tumor == 0
tag2 = res14$C3 == 0 & res14$C5 == 1 & res14$Tumor == 0
tag3 = res14$C3 == 0 & res14$C5 == 0 & res14$Tumor == 1
label = apply(res14,1,sum)
tag4 = label >= 2
tag5 = label == 0

C3_loop = row.names(res14)[tag1]
C5_loop = row.names(res14)[tag2]
Tumor_loop = row.names(res14)[tag3]

C3_loop_res14 = C3_loop
C5_loop_res14 = C5_loop
Tumor_loop_res14 = Tumor_loop

#######################################
#SKCM_211D9CF4
#######################################
tag1 = res15$C3 == 1 & res15$C2 == 0 & res15$Tumor == 0
tag2 = res15$C3 == 0 & res15$C2 == 1 & res15$Tumor == 0
tag3 = res15$C3 == 0 & res15$C2 == 0 & res15$Tumor == 1
label = apply(res15,1,sum)
tag4 = label >= 2
tag5 = label == 0

C3_loop = row.names(res15)[tag1]
C2_loop = row.names(res15)[tag2]
Tumor_loop = row.names(res15)[tag3]

C3_loop_res15 = C3_loop
C2_loop_res15 = C2_loop
Tumor_loop_res15 = Tumor_loop

#######################################
#SKCM_FDA487D2
#######################################
tag1 = res16$C5 == 1 & res16$Tumor == 0
tag2 = res16$C5 == 0 & res16$Tumor == 1

label = apply(res16,1,sum)
tag4 = label >= 2
tag5 = label == 0

C5_loop = row.names(res16)[tag1]
Tumor_loop = row.names(res16)[tag2]

C5_loop_res16 = C5_loop
Tumor_loop_res16 = Tumor_loop 

#Group the pre-set cell-type specific loop together
C1_loop_final = c(C1_loop_res12,C1_loop_res13,C1_loop_res7)
C2_loop_final = c(C2_loop_res15)
C3_loop_final = c(C3_loop_res1,C3_loop_res14,C3_loop_res15,C3_loop_res2,C3_loop_res4)
C4_loop_final = c(C4_loop_res5)
C5_loop_final = c(C5_loop_res1,C5_loop_res10,C5_loop_res11,C5_loop_res12,C5_loop_res14,C5_loop_res16,C5_loop_res2,C5_loop_res3,C5_loop_res4,C5_loop_res5,C5_loop_res6,C5_loop_res7)
C6_loop_final = c(C6_loop_res10,C6_loop_res8,C6_loop_res9)
Tumor_loop_final = c(Tumor_loop_res1,Tumor_loop_res2,Tumor_loop_res3,Tumor_loop_res4,Tumor_loop_res5,Tumor_loop_res6,Tumor_loop_res7,Tumor_loop_res8,Tumor_loop_res9,Tumor_loop_res10,Tumor_loop_res11,Tumor_loop_res12,Tumor_loop_res13,Tumor_loop_res14,Tumor_loop_res15,Tumor_loop_res16)

C1_loop_final = unique(C1_loop_final)
C2_loop_final = unique(C2_loop_final)
C3_loop_final = unique(C3_loop_final)
C4_loop_final = unique(C4_loop_final)
C5_loop_final = unique(C5_loop_final)
C6_loop_final = unique(C6_loop_final)
Tumor_loop_final = unique(Tumor_loop_final)

myoutf1 = "/Immune_loop_full/C1_loop_final.Rda"
myoutf2 = "/Immune_loop_full/C2_loop_final.Rda"
myoutf3 = "/Immune_loop_full/C3_loop_final.Rda"
myoutf4 = "/Immune_loop_full/C4_loop_final.Rda"
myoutf5 = "/Immune_loop_full/C5_loop_final.Rda"
myoutf6 = "/Immune_loop_full/C6_loop_final.Rda"
myoutf7 = "/Immune_loop_full/C7_loop_final.Rda"

save(C1_loop_final,file=myoutf1)
save(C2_loop_final,file=myoutf2)
save(C3_loop_final,file=myoutf3)
save(C4_loop_final,file=myoutf4)
save(C5_loop_final,file=myoutf5)
save(C6_loop_final,file=myoutf6)
save(Tumor_loop_final,file=myoutf7)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to map the sample specific pre-set of cell-type specific loops to the union loop set
#Input : pre-set sample specific pre-set of cell-type specific loops
#Output : pre-set cell-type specific union loops
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Map to union loopset
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)
library(ChIPpeakAnno)
library(ArchR)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$
#GR convert function
#$$$$$$$$$$$$$$$$$$$$$$$$$$$
#define the function
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

#load pre-set loops
myoutf1 = "/Immune_loop_full/C1_loop_final.Rda"
myoutf2 = "/Immune_loop_full/C2_loop_final.Rda"
myoutf3 = "/Immune_loop_full/C3_loop_final.Rda"
myoutf4 = "/Immune_loop_full/C4_loop_final.Rda"
myoutf5 = "/Immune_loop_full/C5_loop_final.Rda"
myoutf6 = "/Immune_loop_full/C6_loop_final.Rda"
myoutf7 = "/Immune_loop_full/C7_loop_final.Rda"

load(myoutf1)
load(myoutf2)
load(myoutf3)
load(myoutf4)
load(myoutf5)
load(myoutf6)
load(myoutf7)

########################################
#C1 interaction create
########################################
chr1 = sapply(C1_loop_final,function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(C1_loop_final,function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(C1_loop_final,function(x) strsplit(x,"_")[[1]][3])

chr2 = sapply(C1_loop_final,function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(C1_loop_final,function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(C1_loop_final,function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))
sta1 = sta1 - 20000
end1 = end1 + 20000

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))
sta2 = sta2 - 20000
end2 = end2 + 20000

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
C1 = GenomicInteractions(anchor_1, anchor_2)

########################################
#C2 interaction create
########################################
chr1 = sapply(C2_loop_final,function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(C2_loop_final,function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(C2_loop_final,function(x) strsplit(x,"_")[[1]][3])

chr2 = sapply(C2_loop_final,function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(C2_loop_final,function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(C2_loop_final,function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))
sta1 = sta1 - 20000
end1 = end1 + 20000

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))
sta2 = sta2 - 20000
end2 = end2 + 20000

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
C2 = GenomicInteractions(anchor_1, anchor_2)

########################################
#C3 interaction create
########################################
chr1 = sapply(C3_loop_final,function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(C3_loop_final,function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(C3_loop_final,function(x) strsplit(x,"_")[[1]][3])
chr2 = sapply(C3_loop_final,function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(C3_loop_final,function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(C3_loop_final,function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))
sta1 = sta1 - 20000
end1 = end1 + 20000

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))
sta2 = sta2 - 20000
end2 = end2 + 20000

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
C3 = GenomicInteractions(anchor_1, anchor_2)

########################################
#C4 interaction create
########################################
chr1 = sapply(C4_loop_final,function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(C4_loop_final,function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(C4_loop_final,function(x) strsplit(x,"_")[[1]][3])
chr2 = sapply(C4_loop_final,function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(C4_loop_final,function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(C4_loop_final,function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))
sta1 = sta1 - 20000
end1 = end1 + 20000

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))
sta2 = sta2 - 20000
end2 = end2 + 20000

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
C4 = GenomicInteractions(anchor_1, anchor_2)

########################################
#C5 interaction create
########################################
chr1 = sapply(C5_loop_final,function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(C5_loop_final,function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(C5_loop_final,function(x) strsplit(x,"_")[[1]][3])
chr2 = sapply(C5_loop_final,function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(C5_loop_final,function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(C5_loop_final,function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))
sta1 = sta1 - 20000
end1 = end1 + 20000

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))
sta2 = sta2 - 20000
end2 = end2 + 20000

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
C5 = GenomicInteractions(anchor_1, anchor_2)

########################################
#C7 interaction create
########################################
chr1 = sapply(Tumor_loop_final,function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(Tumor_loop_final,function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(Tumor_loop_final,function(x) strsplit(x,"_")[[1]][3])
chr2 = sapply(Tumor_loop_final,function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(Tumor_loop_final,function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(Tumor_loop_final,function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))
sta1 = sta1 - 20000
end1 = end1 + 20000

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))
sta2 = sta2 - 20000
end2 = end2 + 20000

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
C7 = GenomicInteractions(anchor_1, anchor_2)

########################################
#union interaction create
########################################
tmpinf = "/Loop_immune_correlation_Spearman.txt"
data = read.table(tmpinf,sep="\t",quote=NULL)
tag = is.na(data)
data[tag]=0

chr1 = sapply(row.names(data),function(x) strsplit(x,"_")[[1]][1])
sta1 = sapply(row.names(data),function(x) strsplit(x,"_")[[1]][2])
end1 = sapply(row.names(data),function(x) strsplit(x,"_")[[1]][3])
chr2 = sapply(row.names(data),function(x) strsplit(x,"_")[[1]][4])
sta2 = sapply(row.names(data),function(x) strsplit(x,"_")[[1]][5])
end2 = sapply(row.names(data),function(x) strsplit(x,"_")[[1]][6])

chr1 = as.vector(chr1)
sta1 = as.numeric(as.vector(sta1))
end1 = as.numeric(as.vector(end1))

chr2 = as.vector(chr2)
sta2 = as.numeric(as.vector(sta2))
end2 = as.integer(as.vector(end2))

anchor_1 = data.frame(V1=chr1,V2=sta1,V3=end1)
anchor_2 = data.frame(V1=chr2,V2=sta2,V3=end2)
anchor_1 = peakDF2GRanges(anchor_1)
anchor_2 = peakDF2GRanges(anchor_2)
All = GenomicInteractions(anchor_1, anchor_2)

#map to the union loop
op_1 = findOverlaps(C1,All)
op_2 = findOverlaps(C2,All)
op_3 = findOverlaps(C3,All)
op_4 = findOverlaps(C4,All)
op_5 = findOverlaps(C5,All)
op_6 = findOverlaps(C7,All)

op_1 = as.data.frame(op_1)
tag1 = unique(op_1$subjectHits)
C1_op = All[tag1]
C1_op = as.data.frame(C1_op)
C1_loop = paste0(C1_op$seqnames1,"_",C1_op$start1,"_",C1_op$end1,"_",C1_op$seqnames2,"_",C1_op$start2,"_",C1_op$end2)
myoutf = "/Specific_loop_General_Full/C1_Full.Rda"
save(C1_loop,file=myoutf)

op_2 = as.data.frame(op_2)
tag2 = unique(op_2$subjectHits)
C2_op = All[tag2]
C2_op = as.data.frame(C2_op)
C2_loop = paste0(C2_op$seqnames1,"_",C2_op$start1,"_",C2_op$end1,"_",C2_op$seqnames2,"_",C2_op$start2,"_",C2_op$end2)
myoutf = "/Specific_loop_General_Full/C2_Full.Rda"
save(C2_loop,file=myoutf)

op_3 = as.data.frame(op_3)
tag3 = unique(op_3$subjectHits)
C3_op = All[tag3]
C3_op = as.data.frame(C3_op)
C3_loop = paste0(C3_op$seqnames1,"_",C3_op$start1,"_",C3_op$end1,"_",C3_op$seqnames2,"_",C3_op$start2,"_",C3_op$end2)
myoutf = "/Specific_loop_General_Full/C3_Full.Rda"
save(C3_loop,file=myoutf)

op_4 = as.data.frame(op_4)
tag4 = unique(op_4$subjectHits)
C4_op = All[tag4]
C4_op = as.data.frame(C4_op)
C4_loop = paste0(C4_op$seqnames1,"_",C4_op$start1,"_",C4_op$end1,"_",C4_op$seqnames2,"_",C4_op$start2,"_",C4_op$end2)
myoutf = "/Specific_loop_General_Full/C4_Full.Rda"
save(C4_loop,file=myoutf)

op_5 = as.data.frame(op_5)
tag5 = unique(op_5$subjectHits)
C5_op = All[tag5]
C5_op = as.data.frame(C5_op)
C5_loop = paste0(C5_op$seqnames1,"_",C5_op$start1,"_",C5_op$end1,"_",C5_op$seqnames2,"_",C5_op$start2,"_",C5_op$end2)
myoutf = "/Specific_loop_General_Full/C5_Full.Rda"
save(C5_loop,file=myoutf)

op_6 = as.data.frame(op_6)
tag6 = unique(op_6$subjectHits)
C7_op = All[tag6]
C7_op = as.data.frame(C7_op)
C7_loop = paste0(C7_op$seqnames1,"_",C7_op$start1,"_",C7_op$end1,"_",C7_op$seqnames2,"_",C7_op$start2,"_",C7_op$end2)
myoutf = "/Specific_loop_General_Full/C7_Full.Rda"
save(C7_loop,file=myoutf)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to filtered the union cell-type specific loop set based on the correlation between interaction signal and actual cell fraction in tumor microenvironment
#Input : pre-set union cell-type specific loops
#Output : filtered union cell-type specific loops
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step1 filter of cell-type specific union loop
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

################################################
#load the pre-set union cell-type specific loop
################################################
load("/Specific_loop_General_Full/C1_Full.Rda")
load("/Specific_loop_General_Full/C2_Full.Rda")
load("/Specific_loop_General_Full/C3_Full.Rda")
load("/Specific_loop_General_Full/C4_Full.Rda")
load("/Specific_loop_General_Full/C5_Full.Rda")
load("/Specific_loop_General_Full/C7_Full.Rda")


C1_loop = paste0(C1_op$seqnames1,"_",C1_op$start1,"_",C1_op$end1,"_",C1_op$seqnames2,"_",C1_op$start2,"_",C1_op$end2)
C2_loop = paste0(C2_op$seqnames1,"_",C2_op$start1,"_",C2_op$end1,"_",C2_op$seqnames2,"_",C2_op$start2,"_",C2_op$end2)
C3_loop = paste0(C3_op$seqnames1,"_",C3_op$start1,"_",C3_op$end1,"_",C3_op$seqnames2,"_",C3_op$start2,"_",C3_op$end2)
C4_loop = paste0(C4_op$seqnames1,"_",C4_op$start1,"_",C4_op$end1,"_",C4_op$seqnames2,"_",C4_op$start2,"_",C4_op$end2)
C5_loop = paste0(C5_op$seqnames1,"_",C5_op$start1,"_",C5_op$end1,"_",C5_op$seqnames2,"_",C5_op$start2,"_",C5_op$end2)
C7_loop = paste0(C7_op$seqnames1,"_",C7_op$start1,"_",C7_op$end1,"_",C7_op$seqnames2,"_",C7_op$start2,"_",C7_op$end2)

################################################
#load the correlation between union loop and actual cell fraction
################################################
tmpinf = "/Loop_immune_correlation_Spearman.txt"
data = read.table(tmpinf,sep="\t",quote=NULL)
tag = is.na(data)
data[tag]=0

tag1 = data$C1 >= 0.3
tag2 = data$C2 >= 0.3
tag3 = data$C3 >= 0.3
tag4 = data$C4 >= 0.3
tag5 = data$C5 >= 0.3
tag7 = data$C7 >= 0.3

################################################
#C1 filter
################################################
C1_All = All[tag1]
op = findOverlaps(C1,C1_All)
op = as.data.frame(op)

C1_op = C1_All[unique(op$subjectHits)]
C1_op = as.data.frame(C1_op)

myoutf1 = "/Specific_loop_General/C1_General_loop.Rda"
save(C1_op,file=myoutf1)

################################################
#C2 filter
################################################
C2_All = All[tag2]
op = findOverlaps(C2,C2_All)
op = as.data.frame(op)

C2_op = C2_All[unique(op$subjectHits)]
C2_op = as.data.frame(C2_op)

myoutf2 = "/Specific_loop_General/C2_General_loop.Rda"
save(C2_op,file=myoutf2)

################################################
#C3 filter
################################################
C3_All = All[tag3]
op = findOverlaps(C3,C3_All)
op = as.data.frame(op)

C3_op = C3_All[unique(op$subjectHits)]
C3_op = as.data.frame(C3_op)

myoutf2 = "/Specific_loop_General/C3_General_loop.Rda"
save(C3_op,file=myoutf2)

################################################
#C4 filter
################################################
C4_All = All[tag4]
op = findOverlaps(C4,C4_All)
op = as.data.frame(op)

C4_op = C4_All[unique(op$subjectHits)]
C4_op = as.data.frame(C4_op)

myoutf2 = "/Specific_loop_General/C4_General_loop.Rda"
save(C4_op,file=myoutf2)

################################################
#C5 filter
################################################
C5_All = All[tag5]
op = findOverlaps(C5,C5_All)
op = as.data.frame(op)

op_C5 = findOverlaps(C5,All)
op_C5 = as.data.frame(op_C5)

C5_op = C5_All[unique(op$subjectHits)]
C5_op = as.data.frame(C5_op)

myoutf2 = "/Specific_loop_General/C5_General_loop.Rda"
save(C5_op,file=myoutf2)

################################################
#C7 filter
################################################

C7_All = All[tag7]
op = findOverlaps(C7,C7_All)
op = as.data.frame(op)

C7_op = C7_All[unique(op$subjectHits)]
C7_op = as.data.frame(C7_op)

myoutf2 = "/Specific_loop_General/C7_General_loop.Rda"
save(C7_op,file=myoutf2)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to further filtered the union cell-type specific loop set based on the correlation between interaction signal and RNA-seq estimated immune cell fraction in tumor microenvironment
#Input : filtered pre-set union cell-type specific loops
#Output : final union cell-type specific loops
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step2 filter of cell-type specific union loop
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

library(GenomicInteractions)
library(GenomicRanges)
library(InteractionSet)
library(ChIPpeakAnno)
library(ArchR)

################################################
#load the step1 filtered pre-set union cell-type specific loop
################################################

myinf1 = "/Specific_loop_General/C1_General_loop.Rda"
myinf2 = "/Specific_loop_General/C2_General_loop.Rda"
myinf3 = "/Specific_loop_General/C3_General_loop.Rda"
myinf4 = "/Specific_loop_General/C4_General_loop.Rda"
myinf5 = "/Specific_loop_General/C5_General_loop.Rda"
myinf7 = "/Specific_loop_General/C7_General_loop.Rda"

load(myinf1)
load(myinf2)
load(myinf3)
load(myinf4)
load(myinf5)
load(myinf7)

C1_loop = paste0(C1_op$seqnames1,"_",C1_op$start1,"_",C1_op$end1,"_",C1_op$seqnames2,"_",C1_op$start2,"_",C1_op$end2)
C2_loop = paste0(C2_op$seqnames1,"_",C2_op$start1,"_",C2_op$end1,"_",C2_op$seqnames2,"_",C2_op$start2,"_",C2_op$end2)
C3_loop = paste0(C3_op$seqnames1,"_",C3_op$start1,"_",C3_op$end1,"_",C3_op$seqnames2,"_",C3_op$start2,"_",C3_op$end2)
C4_loop = paste0(C4_op$seqnames1,"_",C4_op$start1,"_",C4_op$end1,"_",C4_op$seqnames2,"_",C4_op$start2,"_",C4_op$end2)
C5_loop = paste0(C5_op$seqnames1,"_",C5_op$start1,"_",C5_op$end1,"_",C5_op$seqnames2,"_",C5_op$start2,"_",C5_op$end2)
C7_loop = paste0(C7_op$seqnames1,"_",C7_op$start1,"_",C7_op$end1,"_",C7_op$seqnames2,"_",C7_op$start2,"_",C7_op$end2)

################################################
#load the correlation between union loop and RNA-seq estimated immune cell fraction
################################################

tmpinf = "/Loop_infiltration_correlation.txt"
info = read.table(tmpinf,sep="\t",quote=NULL)
info = round(info,2)

################################################
#Step2 filter on correlation
################################################

tag1 = info[C1_loop,"Leukocyte.Fraction"] >= 0.25 
tag2 = info[C2_loop,"Leukocyte.Fraction"] >= 0.25
tag3 = info[C3_loop,"Leukocyte.Fraction"] >= 0.25
tag4 = info[C4_loop,"Leukocyte.Fraction"] >= 0.25
tag5 = info[C5_loop,"Leukocyte.Fraction"] >= 0.25
tag7 = info[C7_loop,"Leukocyte.Fraction"] <= (-1)*0.25

C1_loop_final = C1_loop[tag1]
C2_loop_final = C2_loop[tag2]
C3_loop_final = C3_loop[tag3]
C4_loop_final = C4_loop[tag4]
C5_loop_final = C5_loop[tag5]
C7_loop_final = C7_loop[tag7]

#Final loopset output
myoutf1 = "/Specific_loop_General_Final/C1_General_loop.Rda"
myoutf2 = "/Specific_loop_General_Final/C2_General_loop.Rda"
myoutf3 = "/Specific_loop_General_Final/C3_General_loop.Rda"
myoutf4 = "/Specific_loop_General_Final/C4_General_loop.Rda"
myoutf5 = "/Specific_loop_General_Final/C5_General_loop.Rda"
myoutf7 = "/Specific_loop_General_Final/C7_General_loop.Rda"

save(C1_loop_final,file = myoutf1)
save(C2_loop_final,file = myoutf2)
save(C3_loop_final,file = myoutf3)
save(C4_loop_final,file = myoutf4)
save(C5_loop_final,file = myoutf5)
save(C7_loop_final,file = myoutf7)

