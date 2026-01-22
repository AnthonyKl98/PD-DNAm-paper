## Study parameters
projectTitle<-"FANS PD AK"
processedBy<-"Complex Disease Epigenetic Group, University of Exeter Medical School"
arrayVersion<-"Illumina EPIC V1 microarray"
tissueType<-"Brain"
cellSorted<-"TRUE"

arrayType<-"V1" #This must be either V1 or V2

## technical variables
techVar<-c("Sentrix_ID","Sentrix_Position")

## biological variables
bioVar<-c("Individual_ID", "Cell_Type", "Sex", "Disease_Status", "Phenotype", "Age", "Institute")

predDistinctCT<-c("NeuN+", "SOX10+", "Double neg")
neunCT<-c("NeuN+")

## QC thresholds
thresBS<-80
intenThres<-500
nvThres<-0.1
sexCheck<-TRUE
snpCheck<-TRUE
ctCheck<-TRUE
perMiss<-2
studentThres<-1.5
nSDThres<-3
pnthres<-0.05
perc<-1
