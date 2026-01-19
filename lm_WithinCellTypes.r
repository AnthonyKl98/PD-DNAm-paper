##---------------------------------------------------------------------#
##
## Title: EWAS within each cell type
##
## Purpose of script: perform DNA methylation association analysis of Parkinson's cases vs controls in cell types separately
##
## Author: Anthony Klokkaris (adapted from Eilis Hannon's script) 
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS<-function(row,QCmetrics){
  
  nullCT<-lm(row ~ QCmetrics$Phenotype + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Institute)
  
  # extract case control main effect
  return(c(summary(nullCT)$coefficients["QCmetrics$PhenotypePD", c(1,2,4)]))  
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)
library(dplyr)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- [path to data directory]

#Cell type: NeuN+, SOX10+ or double negative
cellType <- "NeuN+"    

normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS") 

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData) 

QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "PD_GBA" | Phenotype == "PD_nonGBA", "PD", Phenotype))
QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "Control_GBA" | Phenotype == "Control_nonGBA", "Control", Phenotype))

#Exclude samples with failed CETYGO check
QCmetrics <- QCmetrics[!(QCmetrics$Sample_ID %in% c("P94/05_DOUBLE NEG", "P11/17_SOX10 +", "P47/11_SOX10 +", "P24/17_SOX10 +", "P40/17_SOX10 +", "P2/11_SOX10 +", "P79/10_SOX10 +", "P72/12_DOUBLE NEG", "P16/12_SOX10 +", "P73/15_SOX10 +", "19870835_SOX10 +", "20000117_SOX10 +", "20050096_SOX10 +", "20174934_SOX10 +", "20174929_SOX10 +", "20040076_SOX10 +", "19990275_DOUBLE NEG")),]


print(paste0("running EWAS on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]

# subset beta matrix to analysis samples

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]
 
QCmetrics$Phenotype <- factor(QCmetrics$Phenotype)
QCmetrics$Cell_Type <- factor(QCmetrics$Cell_Type)
QCmetrics$Age <- as.numeric(QCmetrics$Age)
QCmetrics$Sex <- factor(QCmetrics$Sex)
QCmetrics$Institute <- factor(QCmetrics$Institute)
QCmetrics$Array_Plate <- factor(QCmetrics$Array_Plate)


#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 3, byrow = TRUE)

rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("PD_coeff", "PD_SE", "PD_P")  

save(outtab, file = file.path(resPath, paste0(cellType,"lm_results.rdata")))

write.csv(outtab, file = file.path(resPath, paste0(cellType,"lm_results.csv")))
