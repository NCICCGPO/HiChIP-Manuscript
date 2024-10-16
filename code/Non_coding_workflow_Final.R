#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to mileup the sequencing information from the HiChIP bam file
#Input : (1) sample specific HiChIP bam (2) matched somatic mutation profile from WGS (3) fasta file from the reference genome
#Output : sample specific bcf files
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Sample specific HiChIP pile up
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

#create the directory for bash job submission
mydir = "/Global_VAF_calling_HiChIP/"
dir.create(mydir,recursive = TRUE)
setwd(mydir)

#load the hichip bam files
mydir = "/merge_sortp_global_bam/"
folders = list.files(mydir)
tag = grep(".bai",folders)
folders = folders[-tag]
file_name = gsub("HiChIP-","",folders)
file_name = gsub(".bam","",file_name)
file_pt_name = sapply(file_name,function(x) paste0(strsplit(x,"-")[[1]][1],"-",strsplit(x,"-")[[1]][2]))
file_pt_name = as.vector(file_pt_name)

#load sample information
tmpinf = "/Merged.FitHiChIP.interactions_Q0.1_MergeNearContacts_normcounts.txt"
data = read.table(tmpinf,sep="\t",quote=NULL,header=T)
sam = colnames(data)
sam = sam[8:length(sam)]

#load sample annotation information
tmpinf = "/TCGA_sample_annotation.csv"
pt = read.csv(tmpinf,header=T)
pt_id = pt$Library_Name
pt_id = sapply(pt_id,function(x) paste0(strsplit(x,"_")[[1]][1],"_",strsplit(x,"_")[[1]][2]))
pt$pt_id = pt_id
tag = which(pt_id %in% sam)
pt = pt[tag,]
hichip_sample = unique(pt$submitter_id)

#load somatic mutation information
tmpinf = "/gdc_tcga_atac_jamboree_AWG_metadata.tsv"
map = read.table(tmpinf,sep="\t",quote=NULL,header=T)
tag = map$data_type == "raw simple somatic mutation"
map = map[tag,]
tag = which(map$case_submitter_id %in% hichip_sample)
map = map[tag,]
tag = which(map$workflow_type %in% "CaVEMan")
map = map[tag,]
xx = map$file_name
sam = unique(map$case_submitter_id)

#bash job create
for(i in 1 : length(sam))
{
	cat("\r",i)
	
	myoutf1 = paste("job_", sam[i], ".sh", sep="")
	conOut = file(myoutf1, "w")
	
	curLine = c("#!/bin/bash")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = paste0("#SBATCH --job-name=",myoutf1)
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --nodes=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --mem-per-cpu=32GB")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --ntasks-per-node=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --time=20:00:00")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --output=%x.%j.out")
	writeLines(curLine, conOut)
	curLine = c("#SBATCH --error=%x.%j.err")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --partition=howchang,owners,normal,sfgf")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("module load biology")
	writeLines(curLine, conOut)
	
	curLine = c("module load bcftools")
	writeLines(curLine, conOut)
	
	tag = pt$submitter_id == sam[i]
	hichip_id = unique(pt$pt_id[tag])
	hichip_id = gsub("_","-",hichip_id)
	tag = which(file_pt_name %in% hichip_id)
	input_file = folders[tag]
	input_file = paste0(mydir,input_file)
	
	TCGA_ID = sam[i]
	mutation_input = "/WGS_SNV_POS.txt"
	
	output_file = paste0("/Global_HiChIP_pileup_BCF/",TCGA_ID,".bcf")
	
	#use bcftools for pileup
	curLine = paste0("bcftools mpileup -f /GRCh38.primary_assembly.genome.fa ",input_file," -I -R",mutation_input," -Ob -o ",output_file)
	
	writeLines(curLine, conOut)	
	close(conOut)
	
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to call the somatic mutation the sequencing information from the HiChIP pile up files
#Input : sample specific bcf files
#Output : sample specific vcf files
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Sample specific HiChIP vcf calling
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())

#create the directory for bash job submission
mydir = "/Global_VCF_calling_HiChIP/"
dir.create(mydir)
setwd(mydir)

#load the bcf files
mydir = "/Global_HiChIP_pileup_BCF/"
files = list.files(mydir)
sam = gsub(".bcf","",files)

for(i in 1 : length(sam))
{
	cat("\r",i)
	
	myoutf1 = paste("job_", sam[i], ".sh", sep="")
	conOut = file(myoutf1, "w")
	
	curLine = c("#!/bin/bash")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = paste0("#SBATCH --job-name=",myoutf1)
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --nodes=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --mem-per-cpu=32GB")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --ntasks-per-node=1")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --time=20:00:00")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --output=%x.%j.out")
	writeLines(curLine, conOut)
	curLine = c("#SBATCH --error=%x.%j.err")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	curLine = c("#SBATCH --partition=howchang,owners,normal,sfgf")
	writeLines(curLine, conOut)
	writeLines("", conOut)
	
	input_file = paste0(mydir,sam[i],".bcf")
	
	curLine = c("module load biology")
	writeLines(curLine, conOut)
	
	curLine = c("module load bcftools")
	writeLines(curLine, conOut)
	
	
	output_file = paste0("/Global_HiChIP_pileup_VCF/",sam[i],".vcf.gz")
	curLine = paste0("bcftools call -mO z -o ",output_file, " ", input_file)
	
	writeLines(curLine, conOut)	
	close(conOut)
	
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This section is used to examine the sequencing depth and the allele frequency from HiChIP vcf files
#Input : sample specific vcf files
#Output : sample specific mutation summary
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Sample specific mutation summary
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
library(vcfR)

#load sample id
tmpinf = "/Merged.FitHiChIP.interactions_Q0.1_MergeNearContacts_normcounts.txt"
data = read.table(tmpinf,sep="\t",quote=NULL,header=T)
sam = colnames(data)
sam = sam[8:length(sam)]

#load annotation file for hichip
tmpinf = "/TCGA_sample_annotation.csv"
pt = read.csv(tmpinf,header=T)
pt_id = pt$Library_Name
pt_id = sapply(pt_id,function(x) paste0(strsplit(x,"_")[[1]][1],"_",strsplit(x,"_")[[1]][2]))
tag = which(pt_id %in% sam)
pt = pt[tag,]
pt = unique(pt$submitter_id)

#load annotation file for wgs
tmpinf = "/gdc_tcga_atac_jamboree_AWG_metadata.tsv"
map = read.table(tmpinf,sep="\t",quote=NULL,header=T)
tag = map$data_type == "raw simple somatic mutation"
map = map[tag,]
tag = which(map$case_submitter_id %in% pt)
map = map[tag,]
tag = which(map$workflow_type %in% "CaVEMan")
map = map[tag,]
xx = map$file_name
sam = unique(map$case_submitter_id)

tmpinf = "/gdc_sample_sheet.2022-07-19.tsv"
info = read.table(tmpinf,sep="\t",quote=NULL,header=T)
tag = which(info$Case.ID %in% sam)
info = info[tag,]

#load hichip vcf files
mydir1 = "/Global_HiChIP_pileup_VCF/"
#load snv fles
mydir2 = "/all_snv/"

files1 = list.files(mydir1)
files2 = list.files(mydir2)

for(i in 1 : length(sam))
{
	cat("\r",i)
	tmpinf1 = paste0(mydir1,sam[i],".vcf.gz")
	
	tag = map$case_submitter_id == sam[i]
	tmp = map$file_name[tag]
	tmpinf2 = paste0(mydir2,tmp)
	
	res1 = read.vcfR(tmpinf1)
	
	for(k in 1 : length(tmpinf2))
	{
	
		res2 = read.vcfR(tmpinf2[k])
	
		HiChIP_SNP = res1@fix
		DP4 = sapply(HiChIP_SNP[,"INFO"],function(x) strsplit(x,";")[[1]][grep("DP4=",strsplit(x,";")[[1]],fix=T)])
		
		#HiChIP depth information
		DP = HiChIP_SNP[,"INFO"]
		DP = sapply(DP,function(x) strsplit(x,";")[[1]][grep("DP=",strsplit(x,";")[[1]],fix=T)])
		DP = as.vector(DP)
		DP = gsub("DP=","",DP)
		DP = as.numeric(DP)
		
		#HiChIP mutation read and reference read information
		DP = gsub("DP=","",DP)
		DP4 = gsub("DP4=","",DP4)
		
		Ref_Forward = sapply(DP4,function(x) strsplit(x,",")[[1]][1])
		Ref_Reverse = sapply(DP4,function(x) strsplit(x,",")[[1]][2])
		Alt_Forward = sapply(DP4,function(x) strsplit(x,",")[[1]][3])
		Alt_Reverse = sapply(DP4,function(x) strsplit(x,",")[[1]][4])
		
		Ref_Forward = as.numeric(Ref_Forward)
		Ref_Reverse = as.numeric(Ref_Reverse)
		Alt_Forward = as.numeric(Alt_Forward)
		Alt_Reverse = as.numeric(Alt_Reverse)
		
		Ref = Ref_Forward + Ref_Reverse
		Alt = Alt_Forward + Alt_Reverse
		
		DP = as.numeric(DP)
		Qual_DP = Ref + Alt
		AF = Alt/Qual_DP
		
		HiChIP_SNP = HiChIP_SNP[,c("CHROM","POS","REF","ALT")]
		HiChIP_SNP = as.data.frame(HiChIP_SNP)
		HiChIP_SNP[,"AF"] = AF
		HiChIP_SNP[,"DP"] = DP
		HiChIP_SNP[,"Qual_DP"] = Qual_DP
		HiChIP_SNP[,"Ref_Count"] = Ref
		HiChIP_SNP[,"Alt_Count"] = Alt
		HiChIP_SNP = unique(HiChIP_SNP)
		
		WGS_AF = res2@gt
		WGS_AF = as.data.frame(WGS_AF)
		WGS_AF = sapply(WGS_AF$TUMOR,function(x) strsplit(x,":")[[1]][10])
		WGS_AF = as.numeric(WGS_AF)
		WGS_SNP = res2@fix
		
		#wgs depth
		DP = WGS_SNP[,"INFO"]
		DP_tag = sapply(DP,function(x) strsplit(x,";")[[1]][grep("DP",strsplit(x,";")[[1]])])
		DP_tag = gsub("DP=","",DP_tag)
		DP_tag = as.numeric(DP_tag)
		
		WGS_SNP = WGS_SNP[,c("CHROM","POS","REF","ALT")]
		WGS_SNP = as.data.frame(WGS_SNP)
		WGS_SNP[,"AF"] = WGS_AF
		WGS_SNP[,"DP"] = DP_tag
		
		row.names(HiChIP_SNP) = paste0(HiChIP_SNP$CHROM,"-",HiChIP_SNP$POS,"-",HiChIP_SNP$REF,"-",HiChIP_SNP$ALT)
		row.names(WGS_SNP) = paste0(WGS_SNP$CHROM,"-",WGS_SNP$POS,"-",WGS_SNP$REF,"-",WGS_SNP$ALT)
		com = intersect(row.names(HiChIP_SNP),row.names(WGS_SNP))
	
		HiChIP_SNP = HiChIP_SNP[com,]
		WGS_SNP = WGS_SNP[com,]
	
		tmp = HiChIP_SNP[,c("CHROM","POS","REF","ALT")]
		tmp[,"HiChIP_AF"] = HiChIP_SNP[,"AF"]
		tmp[,"WGS_AF"] = WGS_SNP[,"AF"]
		tmp[,"HiChIP_DP"] = HiChIP_SNP[,"DP"]
		tmp[,"WGS_DP"] = WGS_SNP[,"DP"]
		tmp[,"Qual_DP"] = HiChIP_SNP[,"Qual_DP"]
		tmp[,"Ref_Count"] = HiChIP_SNP[,"Ref_Count"]
		tmp[,"Alt_Count"] = HiChIP_SNP[,"Alt_Count"]
		
		tmp[,"ID"] = rep(sam[i],nrow(tmp))
	
		if(i==1)
		{
			final = tmp
		}else{
	
			final = rbind(final,tmp)
		}
	}
}

final = unique(final)
#output
myoutf = "/Global_HiChIP_WGS_VAF_DP.Rda"
save(final,file=myoutf)

