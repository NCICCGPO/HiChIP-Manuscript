#custom palettes
cancer_types <- c('#FF7D7D','#FF0000','#A30D0D','#602E01','#A37857','#FF6700'
	,'#FFB55C','#FCE017','#C1B62B','#9AC48A','#82FC56','#036D03','#1FBF82'
	,'#83FFFF','#43B9F9','#002BFF','#262C6B','#A983F2','#7E2AD8','#AA059F'
	,'#E26DDF','#FAC0FF','#757575')
names(cancer_types) <- c('ACC','BRCA','SKCM','CESC','GBM','COAD','PCPG','BLCA'
	,'HNSC','MESO','KIRP','ESCA','LUSC','TGCT','LIHC','LUAD','PRAD','LGG'
	,'KIRC','THCA','STAD','UCEC','CHOL')
	
ic10_ic11 <- c('#FF7D7D','#D51F26','#F97E2B','#FFE600','#FFE600','#7BE561','#89288F'
	,'#75F6FC','#83A4FF','#208A42','#DB65D2','#272E6A','#AAAAAA')
names(ic10_ic11) <- c('ic1','ic2','ic3','ic4','ic4ER-','ic4ER+','ic5','ic6','ic7'
	,'ic8','ic9','ic10','Unknown')
	
schmod2 <- c('#208A42','#D51F26','#272E6A','#89288F','#AAAAAA')
names(schmod2) <- c('ERpos_HER2neg_Low_Prolif','ERpos_HER2neg_High_Prolif'
	,'ERneg_HER2neg','HER2pos','Unknown')
	
pam50 <- c('#208A42','#D51F26','#272E6A','#89288F','#F47D2B','#AAAAAA')
names(pam50) <- c('LumA','LumB','Basal','HER2pos','Normal_like','Unknown')

tcga_discrete_palettes <- list(cancer_types = cancer_types, 
                           ic10_ic11 = ic10_ic11, 
                           schmod2 = schmod2,
						   pam50 = pam50)
						   
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2","#2DB7A3"
	,"#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_methylation <- colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC"
	,"#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC'
	, '#EAD397', '#FDB31A','#E42A2A', '#A31D1D'))(100)
pal_atac_bw_heatmap <- colorRampPalette(c("white","#2488F0","#7F3F98","#E22929"
	,"#FCB31A"))(100)
	
tcga_continuous_palettes <- list(pal_rna = pal_rna
	, pal_methylation = pal_methylation
	, pal_atac = pal_atac
	, pal_atac_bw_heatmap = pal_atac_bw_heatmap)