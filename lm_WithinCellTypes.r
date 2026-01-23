##---------------------------------------------------------------------#
##
## Title: EWAS within each cell type
##
## Purpose of script: perform DNA methylation association analysis of Parkinson's cases vs controls in cell types separately
##
## Author: Anthony Klokkaris (adapted from Eilis Hannon) 
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
normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS") 

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData) 

QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "PD_GBA" | Phenotype == "PD_nonGBA", "PD", Phenotype))
QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "Control_GBA" | Phenotype == "Control_nonGBA", "Control", Phenotype))

#Exclude any samples which failed CETYGO check - add samples to exclude here
QCmetrics <- QCmetrics[!(QCmetrics$Sample_ID %in% c("...")),]

QCmetrics$Phenotype <- factor(QCmetrics$Phenotype)
QCmetrics$Age <- as.numeric(QCmetrics$Age)
QCmetrics$Sex <- factor(QCmetrics$Sex)
QCmetrics$Institute <- factor(QCmetrics$Institute)

#Cell type: NeuN+, SOX10+ or Double negative
cellType <- "NeuN+"    
print(paste0("running EWAS on ", cellType, " cell type..."))

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]
 
#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 3, byrow = TRUE)

#----------------------------------------------------------------------#
# LABEL AND SAVE RESULTS
#----------------------------------------------------------------------#

rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("PD_coeff", "PD_SE", "PD_P")  

save(outtab, file = file.path(resPath, paste0(cellType,"_lm_results.rdata")))
write.csv(outtab, file = file.path(resPath, paste0(cellType,"_lm_results.csv")))


