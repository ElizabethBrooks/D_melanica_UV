#!/usr/bin/env Rscript

# R script to perform variance partitioning

# turn off scientific notation
options(scipen = 999)

# import libraries
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(tibble)

# set working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_14Oct2026/variance partitioning"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv", row.names="gene")[ ,1:24]

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/Repos/D_melanica_UV/inputs/expDesign_OlympicsGenotypesGroup.csv", row.names="sample")

# trim the data table to remove summary statistics
countsTable <- head(inputTable, - 5)

# specify the reference level
targets$treatment <- relevel(factor(targets$treatment), ref = "VIS")
targets$group <- relevel(factor(targets$group), ref = "LT")
targets$genotype <- relevel(factor(targets$genotype), ref = "R2")

# retain genes only if it is expressed at a minimum level
#list <- DGEList(countsTable)
#keep <- filterByExpr(list)
#list <- list[keep, , keep.lib.sizes=FALSE]

# filter genes by number of counts
isexpr <- rowSums(cpm(countsTable) > 0.1) >= 5
list <- DGEList(countsTable[isexpr, ])

# how many genes are lost from filtering
nrow(countsTable) - nrow(list)

# use TMM normalization to eliminate composition biases between libraries
#dge <- calcNormFactors(list)
dge <- normLibSizes(list)
# write normalized counts to file
normList <- cpm(dge, normalized.lib.sizes=TRUE)
# add gene row name tag
normList <- as_tibble(normList, rownames = "gene")

# write log transformed normalized counts to file
normListLog <- cpm(dge, normalized.lib.sizes=TRUE, log=TRUE)
normListLog <- as_tibble(normListLog, rownames = "gene")

# setup a design matrix
sample_group_pca <- factor(paste(targets$treatment,targets$genotype,sep="."))

## variance partitioning

# specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ treatment + group + group + treatment:group

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ treatment + group + treatment:group
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart_fixed.jpg")
plotVarPart(sortCols(vp))
dev.off()

# The variable to be tested must be a fixed effect
form <- ~ treatment + group + group:(1|genotype) + treatment:group

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ (1|treatment) + (1|group) + (1|group:genotype) + (1|treatment:group)
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart_random.jpg")
plotVarPart(sortCols(vp))
dev.off()

