##---------------------------------------------------------------------#
##
## Title: EWAS within each cell type using interaction model to assess effect of GBA1 status
##
## Purpose of script: perform DNA methylation association analysis to assess PD and GBA1 main and interaction effects (in cell types separately)
##
## Author: Anthony Klokkaris (adapted from Eilis Hannon's script) 
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#


diagGBA <- function(row, QCmetrics){
  
  
  full.model <- lm(row ~ Phenotype*GBA_Status + Phenotype + GBA_Status  + Age + Sex + Institute, data=QCmetrics)
  null.model <- lm(row ~ Phenotype + GBA_Status  + Age + Sex + Institute, data=QCmetrics) 
  
  stats <- c(
    coef(summary(full.model))['PhenotypePD:GBA_StatusGBA1',c('Estimate','Std. Error','Pr(>|t|)')], 
    coef(summary(full.model))['PhenotypePD',c('Estimate','Std. Error','Pr(>|t|)')], 
    coef(summary(full.model))['GBA_StatusGBA1',c('Estimate','Std. Error','Pr(>|t|)')]
  )
  stats.combo <- stats 
   
  return(stats.combo)
}

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(stats)
library(doParallel)
library(dplyr)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- #[path to data directory] 
normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS/") 

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

#Exclude samples with failed CETYGO check
QCmetrics <- QCmetrics[!(QCmetrics$Sample_ID %in% c("P94/05_DOUBLE NEG", "P11/17_SOX10 +", "P47/11_SOX10 +", "P24/17_SOX10 +", "P40/17_SOX10 +", "P2/11_SOX10 +", "P79/10_SOX10 +", "P72/12_DOUBLE NEG", "P16/12_SOX10 +", "P73/15_SOX10 +", "19870835_SOX10 +", "20000117_SOX10 +", "20050096_SOX10 +", "20174934_SOX10 +", "20174929_SOX10 +", "20040076_SOX10 +", "19990275_DOUBLE NEG")),]

colnames(QCmetrics)[colnames(QCmetrics) == "Phenotype"] <- "GBA_Status"
colnames(QCmetrics)[colnames(QCmetrics) == "Disease_Status"] <- "Phenotype"

QCmetrics <- QCmetrics %>% mutate(GBA_Status = ifelse(GBA_Status == "PD_GBA" | GBA_Status == "Control_GBA", "GBA1", GBA_Status))
QCmetrics <- QCmetrics %>% mutate(GBA_Status = ifelse(GBA_Status == "PD_nonGBA" | GBA_Status == "Control_nonGBA", "non-GBA1", GBA_Status))
QCmetrics$GBA_Status <- factor(QCmetrics$GBA_Status, levels=c('non-GBA1', 'GBA1'))

QCmetrics$Phenotype <- factor(QCmetrics$Phenotype)
QCmetrics$Age <- as.numeric(QCmetrics$Age)
QCmetrics$Sex <- factor(QCmetrics$Sex)
QCmetrics$Institute <- factor(QCmetrics$Institute)

#Cell type: NeuN+, SOX10+ or double negative
cellType <- "NeuN+"
print(paste0("running EWAS on ", cellType, " cell type..."))

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]
dim(celltypeNormbeta) 

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

mod <- 'diagGBA' 
ANOVA.model <- get(mod)

print("Making clusters")
cl<-makeCluster(10)
registerDoParallel(cl)
clusterExport(cl, list(mod))

## Run the function on all probes 
res <- foreach(i=1:nrow(celltypeNormbeta), .combine=rbind, .packages=c('stats'), .verbose=F) %dopar%{
  ANOVA.model(row=as.numeric(celltypeNormbeta[i,]), QCmetrics)}

#----------------------------------------------------------------------#
# LABEL AND SAVE RESULTS
#----------------------------------------------------------------------#

rownames(res) <- rownames(celltypeNormbeta)
columns <- c('_Estimate','_SE','_P')
output.colnames <- c(paste0('PD:GBA', columns), paste0('PD', columns), paste0('GBA', columns)) #keep in same order as extracted from stats
colnames(res) <- output.colnames

save(res, file = file.path(resPath, paste0(cellType,"lm_GBAModel_results.rdata")))
write.csv(res, file = file.path(resPath, paste0(cellType,"lm_GBAModel_results.csv")))