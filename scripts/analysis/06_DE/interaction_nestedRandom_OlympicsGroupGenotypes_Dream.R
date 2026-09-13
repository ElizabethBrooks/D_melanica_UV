#!/usr/bin/env Rscript

# R script to perform DE analysis with DREAM

# turn off scientific notation
options(scipen = 999)

# import libraries
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(tibble)
library(rcartocolor)
library(ggVennDiagram)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_9Sep2026/interaction_nested_random_group_genotype"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv", row.names="gene")[ ,1:24]

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/Repos/D_melanica_UV/inputs/expDesign_OlympicsGenotypesGroup.csv", row.names="sample")

# set input FDR cutoff
cutFDR <- 0.05

# set LFC cut
#cutFC <- 1.2
cutLFC <- log2(1.2)

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
write.table(normList, file="normalizedCounts.csv", sep=",", row.names=FALSE, quote=FALSE)

# write log transformed normalized counts to file
normListLog <- cpm(dge, normalized.lib.sizes=TRUE, log=TRUE)
normListLog <- as_tibble(normListLog, rownames = "gene")
write.table(normListLog, file="normalizedCounts_logTransformed.csv", sep=",", row.names=FALSE, quote=FALSE)

# setup a design matrix
sample_group <- factor(paste(targets$treatment,targets$genotype,sep="."))

# setup points and colors for PCA
points <- c(0,1,2,3,15,16,17,18)
#colors <- rep(c(plotColors[4], plotColors[5], plotColors[6], plotColors[11]), 2)
colors <-  rep(c("#D55E00", plotColors[6], plotColors[4], plotColors[5]), 2)

# create a PCA plot with a legend
png("plotPCA.png", units="in", width=5, height=4, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group], pch=points[sample_group], gene.selection="common")
legend("top", inset=c(0,-0.47), legend=levels(sample_group), pch=points, col=colors, ncol=2)
dev.off()

## Dream analysis

# specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
# note: adjust reference levels
form <- ~ treatment + group + group:(1|genotype) + treatment:group

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ (1|treatment) + (1|group) + (1|group:genotype) + (1|treatment:group)
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart.jpg")
plotVarPart(sortCols(vp))
dev.off()

# estimate regression coefficients
fitmm <- dream(vobjDream, form, targets)

# apply empirical Bayes shrinkage on linear mixed models
fitmm <- eBayes(fitmm)

# examine design matrix
colnames(fitmm$design)

# get all results of hypothesis test on treatment
dge_results <- topTable(fitmm, coef = "treatmentUV", number = nrow(dge))

# export table of DE genes
dge_results_tbl <- as_tibble(dge_results, rownames = "gene")
write.table(dge_results_tbl, file="treatmentUV.csv", sep=",", row.names=FALSE, quote=FALSE)

# get genes < FDR and LFC cutoffs
dge_sig <- topTable(fitmm, coef = "treatmentUV", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

# export table of sig DE genes
dge_sig_tbl <- as_tibble(dge_sig, rownames = "gene")
out_file <- paste("treatmentUV", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
out_file <- paste(out_file, "csv", sep = ".")
write.table(dge_sig_tbl, file=out_file, sep=",", row.names=FALSE, quote=FALSE)

# check the number of sig DE genes
nrow(dge_sig)

# add column for identifying direction of DE gene expression
dge_sig_tbl$colorDE <- plotColors[5]
dge_sig_tbl$alphaDE <- 0.75
# identify significantly up DE genes
dge_sig_tbl$colorDE[dge_sig_tbl$logFC > cutLFC & dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
dge_sig_tbl$alphaDE[dge_sig_tbl$logFC > cutLFC & dge_sig_tbl$adj.P.Val < cutFDR] <- 1
# identify significantly down DE genes
dge_sig_tbl$colorDE[dge_sig_tbl$logFC < (-1*cutLFC) & dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
dge_sig_tbl$alphaDE[dge_sig_tbl$logFC < (-1*cutLFC) & dge_sig_tbl$adj.P.Val < cutFDR] <- 1
# add column with -log10(adj.P.Val) values
dge_sig_tbl$negLog10FDR <- -log10(dge_sig_tbl$adj.P.Val)
# create volcano plot
jpeg("treatmentUV_volcano.jpg")
ggplot(data=dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()

# identify significantly DE genes by FDR & LFC
#dge_sig_tbl_keep <- dge_sig_tbl$adj.P.Val < cutFDR & (dge_sig_tbl$logFC > cutLFC | dge_sig_tbl$logFC < (-1*cutLFC))
# create filtered results table of DE genes
#dge_sig_tbl.filtered <- dge_sig_tbl[dge_sig_tbl_keep,]
# subset counts table by DE gene set
DGESubset_treatment <- dge_sig_tbl[!grepl(plotColors[5], dge_sig_tbl$colorDE),]
DGESubset_treatment.keep <- normListLog$gene %in% DGESubset_treatment$gene
logcountsSubset_treatment <- normListLog[DGESubset_treatment.keep, ]
# format for plotting
gene_names <- logcountsSubset_treatment$gene
logcountsSubset_treatment$gene <- NULL
rownames(logcountsSubset_treatment) <- gene_names
# heatmap of results
jpeg("VISvsUV_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_treatment), margins = c(8, 1), labRow = FALSE)
dev.off()

# retrieve previous results
treatment_edgeR <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes/glmQLF_2WayANOVA_treatment_topTags.csv")
# identify significantly DE genes by FDR and LFC
treatment_edgeR_keep <- treatment_edgeR$FDR < cutFDR & (treatment_edgeR$logFC > cutLFC | treatment_edgeR$logFC < (-1*cutLFC))
# create filtered results table of DE genes
treatment_edgeR.filtered <- treatment_edgeR[treatment_edgeR_keep,]
# create combined list of DE gene names
list_venn <- list(Dream = dge_sig_tbl$gene, 
                  edgeR = treatment_edgeR.filtered$gene)
# create venn diagram
jpeg("treatment_results_venn.jpg")
ggVennDiagram(list_venn, label_alpha=0.25, category.names = c("Dream","edgeR")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()
