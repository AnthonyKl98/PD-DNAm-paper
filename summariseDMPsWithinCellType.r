##---------------------------------------------------------------------#
##
## Title: Classify & Compare EWAS Results
##
## Purpose of script: Characterize EWAS results and produce summary plots
##
## Author: Anthony Klokkaris (adapted from Eilis Hannon's script) 
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE CHARACTERISATION FUNCTION
#----------------------------------------------------------------------#

#pThres either 9e-8 or 1e-6

classifyDMPs<-function(res, pThres = 9e-8){
	# identify DMPs for given cell type
	dmp <- res[,"PD_P"] < pThres
	res<-cbind(res, dmp)
	colnames(res)[ncol(res)]<-c("DMP")
	return(res)
}

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)
library(dplyr)
library(QCEWAS)
library(Haplin)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- #[path to data directory]
normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS")
pos <- position_dodge(0.9)

#----------------------------------------------------------------------#
# LOAD AND PREPARE RESULTS
#----------------------------------------------------------------------#

setwd(dataDir)

load(file.path(resPath, "NeuN+_lm_results.rdata"))
res.lm_neun<-outtab
nameLM_neun <- "NeuN+_lm_results"
remove(outtab)

load(file.path(resPath, "SOX10+_lm_results.rdata"))
res.lm_sox10<-outtab
nameLM_sox10 <- "SOX10+_lm_results"
remove(outtab)

load(file.path(resPath, "Double neg_lm_results.rdata"))
res.lm_dn<-outtab
nameLM_dn <- "Double neg_lm_results"
remove(outtab)

#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#

## filter out cross hybridizing & SNP probes using lists e.g. in reference directory

crosshyb<-read.table(file.path("[REFDIR]/CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter<-read.csv(file.path("[REFDIR]/EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes<-read.table(file.path("[REFDIR]/SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
crosshyb2<-read.csv(file.path("[REFDIR]/Pidsley_SM1.csv"), stringsAsFactors = FALSE)
snpProbes2<-read.csv(file.path("[REFDIR]/Pidsley_SM4.csv"), stringsAsFactors = FALSE)
snpProbes3<-read.csv(file.path("[REFDIR]/Pidsley_SM5.csv"), stringsAsFactors = FALSE)
snpProbes4<-read.csv(file.path("[REFDIR]/Pidsley_SM6.csv"), stringsAsFactors = FALSE)

snpProbes<-snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01),]
snpProbes2<-snpProbes2[which(snpProbes2$AF > 0.01),]
snpProbes3<-snpProbes3[which(snpProbes3$AF > 0.01),]
snpProbes4<-snpProbes4[which(snpProbes4$AF > 0.01),]

dist<-cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist<-apply(dist, 1, min)
snpProbes4<-snpProbes4[which(dist <=10),]

remove<-unique(c(tofilter$IlmnID, crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))
remove<-intersect(remove, rownames(res.lm_neun))
remove<-intersect(remove, rownames(res.lm_sox10))
remove<-intersect(remove, rownames(res.lm_dn))

if(length(remove) > 0){
	res.lm_neun<-res.lm_neun[-match(remove, rownames(res.lm_neun)),]
  res.lm_sox10<-res.lm_sox10[-match(remove, rownames(res.lm_sox10)),]
  res.lm_dn<-res.lm_dn[-match(remove, rownames(res.lm_dn)),]
}

#----------------------------------------------------------------------#
# ANNOTATE RESULTS
#----------------------------------------------------------------------#

# add columns to indicate cell type dmps
res.lm_neun<-classifyDMPs(res.lm_neun)
res.lm_sox10<-classifyDMPs(res.lm_sox10)
res.lm_dn<-classifyDMPs(res.lm_dn)

# add gene annotation
annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))

newfData <- newfData[rownames(res.lm_neun), ]
newfData <- newfData[rownames(res.lm_sox10), ]
newfData <- newfData[rownames(res.lm_dn), ]

res.lm_neun<-cbind(res.lm_neun, newfData[,c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Regulatory_Feature_Group", "Regulatory_Feature_Name", "Relation_to_Island")])
res.lm_sox10<-cbind(res.lm_sox10, newfData[,c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Regulatory_Feature_Group", "Regulatory_Feature_Name", "Relation_to_Island")])
res.lm_dn<-cbind(res.lm_dn, newfData[,c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Regulatory_Feature_Group", "Regulatory_Feature_Name", "Relation_to_Island")])

#----------------------------------------------------------------------#
# WRITE TABLES
#----------------------------------------------------------------------#

dmps.lm_neun<-res.lm_neun[rowSums(as.matrix(res.lm_neun[,"DMP"])) > 0,]
dmps.lm_sox10<-res.lm_sox10[rowSums(as.matrix(res.lm_sox10[,"DMP"])) > 0,]
dmps.lm_dn<-res.lm_dn[rowSums(as.matrix(res.lm_dn[,"DMP"])) > 0,]

save(dmps.lm_neun, file = file.path(resPath, paste0(nameLM_neun, "DMPs_9e-8.rdata")))
save(dmps.lm_sox10, file = file.path(resPath, paste0(nameLM_sox10, "DMPs_9e-8.rdata")))
save(dmps.lm_dn, file = file.path(resPath, paste0(nameLM_dn, "DMPs_9e-8.rdata")))

write.csv(dmps.lm_neun, file.path(resPath, paste0(nameLM_neun, "DMPs_9e-8.csv")))
write.csv(dmps.lm_sox10, file.path(resPath, paste0(nameLM_sox10, "DMPs_9e-8.csv")))
write.csv(dmps.lm_dn, file.path(resPath, paste0(nameLM_dn, "DMPs_9e-8.csv")))

#----------------------------------------------------------------------#
# QQ PLOTS
#----------------------------------------------------------------------#

pdf(file.path(resPath, paste0(nameLM, "QQplots_withinNeun_SOX10_DN.pdf")), width = 12, height = 4)
par(mfrow = c(1,3))
#options(bitmapType="cairo")
#par(family = "serif")
qq(res.lm_neun[,"PD_P"], main = "NeuN+")
qq(res.lm_sox10[,"PD_P"], main = "SOX10+")
qq(res.lm_dn[,"PD_P"], main = "Double neg")
dev.off()

P_lambda(res.lm_neun[,"PD_P"])
P_lambda(res.lm_sox10[,"PD_P"])
P_lambda(res.lm_dn[,"PD_P"])

#----------------------------------------------------------------------#
# MANHATTAN PLOTS
#----------------------------------------------------------------------#

res.lm_neun[,"chr"][which(res.lm_neun[,"chr"] == "chrX")]<-"chr23"
res.lm_neun[,"chr"][which(res.lm_neun[,"chr"] == "chrY")]<-"chr24"
res.lm_neun[,"chr"]<-gsub("chr", "", res.lm_neun[,"chr"])
res.lm_neun[,"chr"]<-as.numeric(as.character(res.lm_neun[,"chr"]))

res<-data.frame("SNP" = rownames(res.lm_neun), "P" = res.lm_neun[,"PD_P"], "CHR" = res.lm_neun[,"chr"],"BP" = res.lm_neun[,"pos"])
res<-na.omit(res)
options(bitmapType="cairo")
jpeg(file.path(resPath, paste0(nameLM_neun, "ManhattanPlot.jpeg")), width = 1100, height = 600, quality = 100)
par(mar = c(4,4,0.5,0.5))
par(family = "serif")
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(1e-6), ylim = c(0, 10), cex.lab = 1.5, cex.axis = 1.3)
dev.off()

res.lm_sox10[,"chr"][which(res.lm_sox10[,"chr"] == "chrX")]<-"chr23"
res.lm_sox10[,"chr"][which(res.lm_sox10[,"chr"] == "chrY")]<-"chr24"
res.lm_sox10[,"chr"]<-gsub("chr", "", res.lm_sox10[,"chr"])
res.lm_sox10[,"chr"]<-as.numeric(as.character(res.lm_sox10[,"chr"]))

res<-data.frame("SNP" = rownames(res.lm_sox10), "P" = res.lm_sox10[,"PD_P"], "CHR" = res.lm_sox10[,"chr"],"BP" = res.lm_sox10[,"pos"])
res<-na.omit(res)
options(bitmapType="cairo")
jpeg(file.path(resPath, paste0(nameLM_sox10, "ManhattanPlot.jpeg")), width = 1100, height = 600, quality = 100)
par(mar = c(4,4,0.5,0.5))
par(family = "serif")
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(1e-6), ylim = c(0, 10), cex.lab = 1.5, cex.axis = 1.3)
dev.off()

res.lm_dn[,"chr"][which(res.lm_dn[,"chr"] == "chrX")]<-"chr23"
res.lm_dn[,"chr"][which(res.lm_dn[,"chr"] == "chrY")]<-"chr24"
res.lm_dn[,"chr"]<-gsub("chr", "", res.lm_dn[,"chr"])
res.lm_dn[,"chr"]<-as.numeric(as.character(res.lm_dn[,"chr"]))

res<-data.frame("SNP" = rownames(res.lm_dn), "P" = res.lm_dn[,"PD_P"], "CHR" = res.lm_dn[,"chr"],"BP" = res.lm_dn[,"pos"])
res<-na.omit(res)
options(bitmapType="cairo")
jpeg(file.path(resPath, paste0(nameLM_dn, "ManhattanPlot.jpeg")), width = 1100, height = 600, quality = 100)
par(mar = c(4,4,0.5,0.5))
par(family = "serif")
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(1e-6), ylim = c(0, 10), cex.lab = 1.5, cex.axis = 1.3)
dev.off()

#----------------------------------------------------------------------#
# PLOT DMPS
#----------------------------------------------------------------------#

load(normData)
QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "PD_GBA" | Phenotype == "PD_nonGBA", "PD", Phenotype))
QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "Control_GBA" | Phenotype == "Control_nonGBA", "Control", Phenotype))

#Exclude samples which failed CETYGO check
QCmetrics <- QCmetrics[!(QCmetrics$Sample_ID %in% c("P94/05_DOUBLE NEG", "P11/17_SOX10 +", "P47/11_SOX10 +", "P24/17_SOX10 +", "P40/17_SOX10 +", "P2/11_SOX10 +", "P79/10_SOX10 +", "P72/12_DOUBLE NEG", "P16/12_SOX10 +", "P73/15_SOX10 +", "19870835_SOX10 +", "20000117_SOX10 +", "20050096_SOX10 +", "20174934_SOX10 +", "20174929_SOX10 +", "20040076_SOX10 +", "19990275_DOUBLE NEG")),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]

# Plot NeuN+ DMPs 

QCmetrics_neun<-QCmetrics[which(QCmetrics$Cell_Type == "NeuN+"),]
celltypeNormbeta_neun<-celltypeNormbeta[,QCmetrics_neun$Basename]

pdf(file.path(resPath, paste0(nameLM_neun, "Final DMPs_NeuNonly_boxplot.pdf")), width = 6, height = 5)
for(i in order(dmps.lm_neun[,"PD_P"], decreasing = FALSE)){
  tmpDat<-data.frame("Phenotype" = QCmetrics_neun$Phenotype, "DNAm" = celltypeNormbeta_neun[rownames(dmps.lm_neun)[i],])
  p<-ggplot(tmpDat, aes(x=Phenotype, y=DNAm, fill = Phenotype)) + 
  geom_boxplot(position = pos, width = 0.6, outlier.shape = NA) +
#     stat_summary(fun = "mean", 
#                geom = "point", 
#                position = pos,
# 			   shape = 18, size = 3, colour = "black") + 
	ggtitle(paste(rownames(dmps.lm_neun)[i], dmps.lm_neun$UCSC_RefGene_Name[i])) +
    geom_point(position = pos, scale = 'width') +
    scale_fill_discrete(name = "PD Status") + 
	  theme(
      axis.title = element_text(size = 18),  # Axis labels size
      axis.text = element_text(size = 14)    # Axis numbers size
    )  	
  print(p)
}
dev.off()

#----------------------------------------------------------------------#
# VOLCANO PLOT
#----------------------------------------------------------------------#

#for NeuN+ EWAS

res.lm_neun$colour <- "black"  # Default color
res.lm_neun$colour[res.lm_neun[,"PD_P"] < 9e-8] <- "red"    # Points above both lines
res.lm_neun$colour[res.lm_neun[,"PD_P"] <= 1e-6 & res.lm_neun[,"PD_P"] > 9e-8] <- "blue"  # Points between the lines

pdf(file.path(resPath, paste0(date, "_Volcano_NeuN_1e-6.pdf")), width = 6, height = 5)
ggplot(res.lm_neun, aes(x = PD_coeff, y = -log10(PD_P), color = colour)) +  # Change to color
  geom_point(shape = 19, size = 1) +  # Use color aesthetic
  geom_hline(yintercept = -log10(1e-6), colour = "blue", linetype = "solid") +  # Blue line at PD_P = 1e-6
  geom_hline(yintercept = -log10(9e-8), colour = "red", linetype = "solid") +
  xlab("Estimate") + ylab("-log10(P)") +
  theme_minimal() +  # Optional: Use a minimal theme
  scale_color_identity() + 
	  theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12))
dev.off()  