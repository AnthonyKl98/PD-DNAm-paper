# PD-DNAm-paper

This folder contains the scripts used for performing epigenome-wide association study analysis. These should be used following the data quality control pipeline on https://github.com/ejh243/BrainFANS/tree/master/array/DNAm/preprocessing. A config.r file for the quality control is provided here to show the parameters used. The quality control pipeline
generates the normalised beta file (normalised.rdata), which is the input for the EWAS analysis.

The EWAS analysis scripts are:

lm_WithinCellTypes.r - DNA methylation association analysis using linear (lm) model in each cell type separately, comparing Parkinson's cases vs controls 
lm_GBAinteractionModel.r -  DNA methylation association analysis using linear (lm) model in each cell type separately, testing for PD and GBA1 main and interaction effects 
summariseDMPsWithinCellType.r - characterising and plotting EWAS results, including DMP identification
mixedEffectsModel.r - DNA methylation association analysis using mixed-effects linear model, testing for main and cell type-specific effects


