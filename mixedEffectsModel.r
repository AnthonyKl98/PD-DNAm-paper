##---------------------------------------------------------------------#
##
## Title: EWAS with mixed effects linear regression model
##
## Purpose of script: perform DNA methylation association analysis of Parkinson's vs controls testing for main and cell-specific effects, 
## simulataneously testing for cell type differences
##
## Author: Anthony Klokkaris (adapted from Eilis Hannon's script) 
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS<-function(row,QCmetrics){

	pheno<-cbind(row,QCmetrics)
	modelMLM<-lmer(row ~ Phenotype * Cell_Type + Age + Sex + (1 | Institute) + (1 | Individual_ID), REML = FALSE, data = pheno)
	nullMLM<-lmer(row ~ Phenotype + Cell_Type + Age + Sex + (1 | Institute) +  (1 | Individual_ID), REML = FALSE, data = pheno)
	nullCT<-lmer(row ~ Phenotype + Age + Sex + (1 | Institute) + (1 | Individual_ID), REML = FALSE, data = pheno)
	
	# extract case control main effect
	return(c(summary(modelMLM)$coefficients["PhenotypePD",c(1,2,5)],
	
	# extract cell specific case control effect
	anova(modelMLM,nullMLM)[2,8], 
	summary(modelMLM)$coefficients["PhenotypePD:Cell_TypeNeuN+",c(1,2,5)],
	summary(modelMLM)$coefficients["PhenotypePD:Cell_TypeSOX10+",c(1,2,5)],
	
	# extract cell type effect
	anova(nullMLM, nullCT)[2,8],
	summary(nullMLM)$coefficients["Cell_TypeNeuN+",c(1,2,5)],
	summary(nullMLM)$coefficients["Cell_TypeSOX10+",c(1,2,5)]))
}

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(dplyr)
library(lme4)
library(lmerTest)
library(doParallel)
library(devtools)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- #[path to data directory]
normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "PD_GBA" | Phenotype == "PD_nonGBA", "PD", Phenotype))
QCmetrics <- QCmetrics %>% mutate(Phenotype = ifelse(Phenotype == "Control_GBA" | Phenotype == "Control_nonGBA", "Control", Phenotype))

#Exclude samples which failed CETYGO check
QCmetrics <- QCmetrics[!(QCmetrics$Sample_ID %in% c("P94/05_DOUBLE NEG", "P11/17_SOX10 +", "P47/11_SOX10 +", "P24/17_SOX10 +", "P40/17_SOX10 +", "P2/11_SOX10 +", "P79/10_SOX10 +", "P72/12_DOUBLE NEG", "P16/12_SOX10 +", "P73/15_SOX10 +", "19870835_SOX10 +", "20000117_SOX10 +", "20050096_SOX10 +", "20174934_SOX10 +", "20174929_SOX10 +", "20040076_SOX10 +", "19990275_DOUBLE NEG")),]

QCmetrics$Phenotype <- factor(QCmetrics$Phenotype)
QCmetrics$Cell_Type <- factor(QCmetrics$Cell_Type)
QCmetrics$Age <- as.numeric(QCmetrics$Age)
QCmetrics$Sex <- factor(QCmetrics$Sex)
QCmetrics$Institute <- factor(QCmetrics$Institute)

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
cellTypes<-unique(QCmetrics$Cell_Type)

celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]
dim(celltypeNormbeta) 

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "lmer"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 17, byrow = TRUE)

#----------------------------------------------------------------------#
# LABEL AND SAVE RESULTS
#----------------------------------------------------------------------#

rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("PD_coeff", "PD_SE", "PD_P", "CellType_PD_P", "NeuN_PD_coeff", "NeuN_PD_SE","NeuN_PD_P", "SOX10_PD_coeff", "SOX10_PD_SE","SOX10_PD_P", "CellType_P","NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P")  

save(outtab, file = file.path(resPath, "mlm_results.rdata")) 
write.csv(outtab, file = file.path(resPath, "mlm_results.csv"))