#!/usr/bin/env Rscript

# R script to perform DE analysis with limma
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7873980/
# https://diseaseneurogenomics.github.io/variancePartition/articles/dream.html

# turn off scientific notation
options(scipen = 999)

# import libraries
library(variancePartition)
library(edgeR)
#library(BiocParallel)
library(tibble)
library(rcartocolor)
library(ggVennDiagram)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_14Sep2026/mixed_interaction_group"
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

# setup a design matrix
#sample_group <- factor(paste(targets$treatment,targets$group,sep="."))

# specify the reference level
#sample_group <- relevel(sample_group, ref = "VIS.LT")

# specify the reference level
targets$treatment <- relevel(factor(targets$treatment), ref = "VIS")
targets$group <- relevel(factor(targets$group), ref = "LT")
targets$genotype <- relevel(factor(targets$genotype), ref = "R2")

# experimental design as a one-way layout, where one coefficient is assigned to each group
#design <- model.matrix(~ 0 + sample_group)
#colnames(design) <- levels(sample_group)

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
normLog <- cpm(dge, normalized.lib.sizes=TRUE, log=TRUE)
normListLog <- as_tibble(normLog, rownames = "gene")
write.table(normListLog, file="normalizedCounts_logTransformed.csv", sep=",", row.names=FALSE, quote=FALSE)

# setup a design matrix
sample_group_pca <- factor(paste(targets$treatment,targets$genotype,sep="."))

# setup points and colors for PCA
points <- c(0,1,2,3,15,16,17,18)
#colors <- rep(c(plotColors[4], plotColors[5], plotColors[6], plotColors[11]), 2)
colors <-  rep(c("#D55E00", plotColors[6], plotColors[4], plotColors[5]), 2)

# create a PCA plot with a legend
png("plotPCA.png", units="in", width=5, height=4, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.47), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2)
dev.off()


## limma analysis

# estimate the correlation between repeated measurements
#cor <- duplicateCorrelation(dge, design, block=targets$genotype)
#cor$consensus.correlation

# fit the random effects
#fit <- lmFit(object=dge, design=design, 
#             block=targets$genotype, correlation=cor$consensus.correlation)

# make contrasts
#contrasts <- makeContrasts(
#  LT_VISvsUV=VIS.LT-UV.LT, 
#  HT_VISvsUV=VIS.HT-UV.HT, 
#  VISvsUV=(VIS.LT-VIS.HT)-(UV.LT-UV.HT), 
#  levels=colnames(design)) 
#contrasts

# incorporate contrast matrix into the limma pipeline
#fit <- contrasts.fit(fit, contrasts)

# setup design matrix
design <- model.matrix(~ treatment + group + treatment:group, targets)

# apply duplicateCorrelation is two rounds
vobj_tmp <- voom(dge, design, plot = FALSE)
dupcor <- duplicateCorrelation(vobj_tmp, design, block=targets$genotype)

# view correlation between repeated measurements
dupcor$consensus.correlation

# run voom again to compute more accurate precision weights
vobj <- voom(dge, design, plot = FALSE, block = targets$genotype, correlation = dupcor$consensus)

# estimate linear mixed model with a single variance component and fit the model for each gene
dupcor <- duplicateCorrelation(vobj, design, block = targets$genotype)

# this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj, design, block = targets$genotype, correlation = dupcor$consensus)

# fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes(fitDupCor)

# get names of available coefficients and contrasts for testing
colnames(fitDupCor)

# get all results of hypothesis test on treatment
treatment_dge_results <- topTable(fitDupCor, coef = "treatmentUV", number = nrow(dge))
tolerance_dge_results <- topTable(fitDupCor, coef = "groupHT", number = nrow(dge))
interaction_dge_results <- topTable(fitDupCor, coef = "treatmentUV:groupHT", number = nrow(dge))

# export table of DE genes
treatment_dge_results_tbl <- as_tibble(treatment_dge_results, rownames = "gene")
write.table(treatment_dge_results_tbl, file="treatmentUV.csv", sep=",", row.names=FALSE, quote=FALSE)
tolerance_dge_results_tbl <- as_tibble(tolerance_dge_results, rownames = "gene")
write.table(tolerance_dge_results_tbl, file="groupHT.csv", sep=",", row.names=FALSE, quote=FALSE)
interaction_dge_results_tbl <- as_tibble(interaction_dge_results, rownames = "gene")
write.table(interaction_dge_results_tbl, file="treatmentUV_groupHT.csv", sep=",", row.names=FALSE, quote=FALSE)

# get genes < FDR and LFC cutoffs
treatment_dge_sig <- topTable(fitDupCor, coef = "treatmentUV", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
tolerance_dge_sig <- topTable(fitDupCor, coef = "groupHT", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
interaction_dge_sig <- topTable(fitDupCor, coef = "treatmentUV:groupHT", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

# export table of sig DE genes
treatment_dge_sig_tbl <- as_tibble(treatment_dge_sig, rownames = "gene")
treatment_out_file <- paste("treatmentUV", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_out_file <- paste(treatment_out_file, "csv", sep = ".")
write.table(treatment_dge_sig_tbl, file=treatment_out_file, sep=",", row.names=FALSE, quote=FALSE)
tolerance_dge_sig_tbl <- as_tibble(tolerance_dge_sig, rownames = "gene")
tolerance_out_file <- paste("groupHT", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
tolerance_out_file <- paste(tolerance_out_file, "csv", sep = ".")
write.table(tolerance_dge_sig_tbl, file=tolerance_out_file, sep=",", row.names=FALSE, quote=FALSE)
interaction_dge_sig_tbl <- as_tibble(interaction_dge_sig, rownames = "gene")
interaction_out_file <- paste("treatmentUV_groupHT", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
interaction_out_file <- paste(interaction_out_file, "csv", sep = ".")
write.table(interaction_dge_sig_tbl, file=interaction_out_file, sep=",", row.names=FALSE, quote=FALSE)

# check the number of sig DE genes
nrow(treatment_dge_sig)
nrow(tolerance_dge_sig)
nrow(interaction_dge_sig)

# add column for identifying direction of DE gene expression
treatment_dge_sig_tbl$colorDE <- plotColors[5]
treatment_dge_sig_tbl$alphaDE <- 0.75
tolerance_dge_sig_tbl$colorDE <- plotColors[5]
tolerance_dge_sig_tbl$alphaDE <- 0.75
#interaction_dge_sig_tbl$colorDE <- plotColors[5]
#interaction_dge_sig_tbl$alphaDE <- 0.75

# identify significantly up DE genes
treatment_dge_sig_tbl$colorDE[treatment_dge_sig_tbl$logFC > cutLFC & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
treatment_dge_sig_tbl$alphaDE[treatment_dge_sig_tbl$logFC > cutLFC & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
tolerance_dge_sig_tbl$colorDE[tolerance_dge_sig_tbl$logFC > cutLFC & tolerance_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
tolerance_dge_sig_tbl$alphaDE[tolerance_dge_sig_tbl$logFC > cutLFC & tolerance_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
#interaction_dge_sig_tbl$colorDE[interaction_dge_sig_tbl$logFC > cutLFC & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
#interaction_dge_sig_tbl$alphaDE[interaction_dge_sig_tbl$logFC > cutLFC & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- 1

# identify significantly down DE genes
treatment_dge_sig_tbl$colorDE[treatment_dge_sig_tbl$logFC < (-1*cutLFC) & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
treatment_dge_sig_tbl$alphaDE[treatment_dge_sig_tbl$logFC < (-1*cutLFC) & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
tolerance_dge_sig_tbl$colorDE[tolerance_dge_sig_tbl$logFC < (-1*cutLFC) & tolerance_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
tolerance_dge_sig_tbl$alphaDE[tolerance_dge_sig_tbl$logFC < (-1*cutLFC) & tolerance_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
#interaction_dge_sig_tbl$colorDE[interaction_dge_sig_tbl$logFC < (-1*cutLFC) & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
#interaction_dge_sig_tbl$alphaDE[interaction_dge_sig_tbl$logFC < (-1*cutLFC) & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- 1

# add column with -log10(adj.P.Val) values
treatment_dge_sig_tbl$negLog10FDR <- -log10(treatment_dge_sig_tbl$adj.P.Val)
tolerance_dge_sig_tbl$negLog10FDR <- -log10(tolerance_dge_sig_tbl$adj.P.Val)
#interaction_dge_sig_tbl$negLog10FDR <- -log10(interaction_dge_sig_tbl$adj.P.Val)

# create volcano plot
jpeg("treatmentUV_volcano.jpg")
ggplot(data=treatment_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
jpeg("groupHT_volcano.jpg")
ggplot(data=tolerance_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
#jpeg("treatmentUV_groupHT_volcano.jpg")
#ggplot(data=interaction_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
#  geom_point() +
#  theme_minimal() +
#  scale_color_identity() +
#  scale_alpha(guide = 'none') +
#  xlab("LFC")
#dev.off()

# subset counts table by DE gene set
DGESubset_treatment <- treatment_dge_sig_tbl[!grepl(plotColors[5], treatment_dge_sig_tbl$colorDE),]
DGESubset_treatment.keep <- normListLog$gene %in% DGESubset_treatment$gene
logcountsSubset_treatment <- normListLog[DGESubset_treatment.keep, ]
DGESubset_tolerance <- tolerance_dge_sig_tbl[!grepl(plotColors[5], tolerance_dge_sig_tbl$colorDE),]
DGESubset_tolerance.keep <- normListLog$gene %in% DGESubset_tolerance$gene
logcountsSubset_tolerance <- normListLog[DGESubset_tolerance.keep, ]
#DGESubset_interaction <- interaction_dge_sig_tbl[!grepl(plotColors[5], interaction_dge_sig_tbl$colorDE),]
#DGESubset_interaction.keep <- normListLog$gene %in% DGESubset_interaction$gene
#logcountsSubset_interaction <- normListLog[DGESubset_interaction.keep, ]

# format for plotting
gene_names <- logcountsSubset_treatment$gene
logcountsSubset_treatment$gene <- NULL
rownames(logcountsSubset_treatment) <- gene_names
gene_names <- logcountsSubset_tolerance$gene
logcountsSubset_tolerance$gene <- NULL
rownames(logcountsSubset_tolerance) <- gene_names
#gene_names <- logcountsSubset_interaction$gene
#logcountsSubset_interaction$gene <- NULL
#rownames(logcountsSubset_interaction) <- gene_names

# heatmap of results
jpeg("treatmentUV_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_treatment), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("groupHT_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_tolerance), margins = c(8, 1), labRow = FALSE)
dev.off()
#jpeg("treatmentUV_groupHT_heatmap.jpg")
#heatmap(as.matrix(logcountsSubset_interaction), margins = c(8, 1), labRow = FALSE)
#dev.off()

# retrieve previous results
treatment_edgeR <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes/glmQLF_2WayANOVA_treatment_topTags.csv")
group_edgeR <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes/glmQLF_2WayANOVA_tolerance_topTags.csv")
#interaction_edgeR <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes/glmQLF_2WayANOVA_interaction_topTags.csv")

# identify significantly DE genes by FDR and LFC
treatment_edgeR_keep <- treatment_edgeR$FDR < cutFDR & (treatment_edgeR$logFC > cutLFC | treatment_edgeR$logFC < (-1*cutLFC))
group_edgeR_keep <- group_edgeR$FDR < cutFDR & (group_edgeR$logFC > cutLFC | group_edgeR$logFC < (-1*cutLFC))
#interaction_edgeR_keep <- interaction_edgeR$FDR < cutFDR & (interaction_edgeR$logFC > cutLFC | interaction_edgeR$logFC < (-1*cutLFC))

# create filtered results table of DE genes
treatment_edgeR.filtered <- treatment_edgeR[treatment_edgeR_keep,]
tolerance_edgeR.filtered <- group_edgeR[group_edgeR_keep,]
#interaction_edgeR.filtered <- interaction_edgeR[interaction_edgeR_keep,]

# create combined list of DE gene names
treatment_list_venn <- list(limma = treatment_dge_sig_tbl$gene, edgeR = treatment_edgeR.filtered$gene)
tolerance_list_venn <- list(limma = tolerance_dge_sig_tbl$gene, edgeR = tolerance_edgeR.filtered$gene)
#interaction_list_venn <- list(limma = interaction_dge_sig_tbl$gene, edgeR = interaction_edgeR.filtered$gene)

# create venn diagram
jpeg("treatment_results_venn.jpg")
ggVennDiagram(treatment_list_venn, label_alpha=0.25, category.names = c("limma","edgeR")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()
# create venn diagram
jpeg("tolerance_results_venn.jpg")
ggVennDiagram(tolerance_list_venn, label_alpha=0.25, category.names = c("limma","edgeR")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()# create venn diagram
#jpeg("interaction_results_venn.jpg")
#ggVennDiagram(interaction_list_venn, label_alpha=0.25, category.names = c("limma","edgeR")) +
#  scale_colour_discrete(type = plotColorSubset)
#dev.off()
