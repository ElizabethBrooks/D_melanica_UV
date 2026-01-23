#!/usr/bin/env Rscript
#Usage: Rscript glmQLF_OlympicsGenotypes_edgeR.r workingDir countsFile startColumn endColumn factorGroupingFile
# usage ex: Rscript glmQLF_OlympicsGenotypes_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/DEAnalysis/Genotypes /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/geneCounts_merged_genome_counted_htseq_run1.txt 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypes.csv
#Usage Ex: Rscript glmQLF_OlympicsGenotypes_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/GeneCountsAnalyzed_KAP4/Formatted/cleaned.csv 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypes.csv
#Usage Ex: Rscript glmQLF_OlympicsGenotypes_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypes.csv
#R script to perform statistical analysis of gene count tables using edgeR GLM

# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
# https://support.bioconductor.org/p/132926/
# https://support.bioconductor.org/p/106608/

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Turn off scientific notation
options(scipen = 999)

#Load the edgeR library
library(edgeR)
library(statmod)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(rcartocolor)
library(tibble)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1]
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes"
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_28May2025/Genotypes"
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes"
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes"
dir.create(workingDir)
setwd(workingDir)

#Import gene count data
#inputTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv", row.names="gene")[ ,1:24]

#Trim the data table
countsTable <- head(inputTable, - 5)

#Import grouping factor
#targets <- read.csv(file=args[5], row.names="sample")
#targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypes.csv", row.names="sample")
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/inputs/expDesign_OlympicsGenotypes_IDs.csv", row.names="sample")

# set input FDR cutoff
#cutFDR <- as.numeric(args[6])
cutFDR <- 0.05

# set LFC cut
#cutLFC <- as.numeric(args[7])
cutLFC <- log2(1.2)

#Setup a design matrix
sample_group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=sample_group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=sample_group)
colnames(list) <- rownames(targets)

#Plot the library sizes before normalization
jpeg("glmQLF_plotBarsBefore.jpg")
barplot(list$samples$lib.size*1e-6, names=1:24, ylab="Library size (millions)")
dev.off()

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
# add gene row name tag
normList <- as_tibble(normList, rownames = "gene")
write.table(normList, file="glmQLF_normalizedCounts.csv", sep=",", row.names=FALSE, quote=FALSE)

#Write log transformed normalized counts to file
normListLog <- cpm(list, normalized.lib.sizes=TRUE, log=TRUE)
normListLog <- as_tibble(normListLog, rownames = "gene")
write.table(normListLog, file="glmQLF_normalizedCounts_logTransformed.csv", sep=",", row.names=FALSE, quote=FALSE)

#Verify TMM normalization using a MD plot
#Write plot to file
jpeg("glmQLF_plotMDBefore.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col=plotColorSubset[3], lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
#colors <- rep(c(plotColors[4], plotColors[5], plotColors[6], plotColors[11]), 2)
colors <-  rep(c("#D55E00", plotColors[6], plotColors[4], plotColors[5]), 2)
#Write plot with legend to file
jpeg("glmQLF_plotMDS.jpg")
par(mar=c(5.1, 4.1, 4.1, 11.1), xpd=TRUE)
plotMDS(list, col=colors[sample_group], pch=points[sample_group])
legend("topright", inset=c(-0.8,0), legend=levels(sample_group), pch=points, col=colors, ncol=2)
#legend("topleft", legend=levels(sample_group), pch=points, col=colors, ncol=2)
dev.off()
#Write plot without legend to file
jpeg("glmQLF_plotMDS_noLegend.jpg")
plotMDS(list, col=colors[sample_group], pch=points[sample_group])
dev.off()

# Create a PCA plot with a legend
png("glmQLF_plotPCA.png", units="in", width=5, height=4, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(list, col=colors[sample_group], pch=points[sample_group], gene.selection="common")
legend("top", inset=c(0,-0.47), legend=levels(sample_group), pch=points, col=colors, ncol=2)
#legend("topleft", legend=levels(sample_group), pch=points, col=colors, ncol=2)
dev.off()

# Create a PCA plot without a legend
png("glmQLF_plotPCA_noLegend.png", units="in", width=5, height=4, res=300)
plotMDS(list, col=colors[sample_group], pch=points[sample_group], gene.selection="common")
dev.off()

##
#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + sample_group)
colnames(design) <- levels(sample_group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
#Write plot to file
jpeg("glmQLF_plotBCV.jpg")
plotBCV(list)
dev.off()

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)
#Write plot to file
jpeg("glmQLF_plotQLDisp.jpg")
plotQLDisp(fit)
dev.off()

# view column order
colnames(fit)


# testing explicit nested contrasts
con.all.nest <- makeContrasts(treatment = (UV.LT_C + UV.LT_A + UV.HT_I + UV.HT_K)/4
                              - (VIS.LT_C + VIS.LT_A + VIS.HT_I + VIS.HT_K)/4,
                              tolerance = (UV.HT_I + UV.HT_K + VIS.HT_I + VIS.HT_K)/4
                              - (UV.LT_C + UV.LT_A + VIS.LT_C + VIS.LT_A)/4,
                              levels=design)
# treatment
treat.anov.treatment <- glmTreat(fit, contrast=con.all.nest[,"treatment"], lfc=cutLFC)
summary(decideTests(treat.anov.treatment))
# tolerance
treat.anov.tolerance <- glmTreat(fit, contrast=con.all.nest[,"tolerance"], lfc=cutLFC)
summary(decideTests(treat.anov.tolerance))
# interaction
treat.anov.Inter <- glmTreat(fit, contrast=c(con.all.nest[,"treatment"]-con.all.nest[,"tolerance"]), lfc=cutLFC)
summary(decideTests(treat.anov.Inter))


# export tables of DE genes
#Write tags table of DE genes to file
tagsTblANOVATreatment <- topTags(treat.anov.treatment, n=nrow(treat.anov.treatment$table), adjust.method="fdr")$table
tagsTblANOVATreatment <- as_tibble(tagsTblANOVATreatment, rownames = "gene")
write.table(tagsTblANOVATreatment, file="glmQLF_2WayANOVA_treatment_topTags.csv", sep=",", row.names=FALSE, quote=FALSE)

#Write tags table of DE genes to file
tagsTblANOVATolerance <- topTags(treat.anov.tolerance, n=nrow(treat.anov.tolerance$table), adjust.method="fdr")$table
tagsTblANOVATolerance <- as_tibble(tagsTblANOVATolerance, rownames = "gene")
write.table(tagsTblANOVATolerance, file="glmQLF_2WayANOVA_tolerance_topTags.csv", sep=",", row.names=FALSE, quote=FALSE)

#Generate table of DE genes
tagsTblANOVAInter <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table), adjust.method="fdr")$table
tagsTblANOVAInter <- as_tibble(tagsTblANOVAInter, rownames = "gene")
write.table(tagsTblANOVAInter, file="glmQLF_2WayANOVA_interaction_topTags.csv", sep=",", row.names=FALSE, quote=FALSE)


# MD plots
#Write plot to file
jpeg("glmQLF_2WayANOVA_treatment_plotMD.jpg")
plotMD(treat.anov.treatment)
abline(h=c((-1*cutLFC), cutLFC), col="blue")
dev.off()

#Write plot to file
jpeg("glmQLF_2WayANOVA_tolerance_plotMD.jpg")
plotMD(treat.anov.tolerance)
abline(h=c((-1*cutLFC), cutLFC), col="blue")
dev.off()

#Write plot to file
jpeg("glmQLF_2WayANOVA_interaction_plotMD.jpg")
plotMD(treat.anov.Inter)
abline(h=c((-1*cutLFC), cutLFC), col="blue")
dev.off()


# Volcano plots
# add column for identifying direction of DE gene expression
tagsTblANOVATreatment$colorDE <- plotColors[5]
tagsTblANOVATreatment$alphaDE <- 0.75
# identify significantly up DE genes
tagsTblANOVATreatment$colorDE[tagsTblANOVATreatment$logFC > cutLFC & tagsTblANOVATreatment$FDR < cutFDR] <- plotColors[4]
tagsTblANOVATreatment$alphaDE[tagsTblANOVATreatment$logFC > cutLFC & tagsTblANOVATreatment$FDR < cutFDR] <- 1
#tagsTblANOVATreatment$colorDE[sign(tagsTblANOVATreatment$logFC) == 1 & tagsTblANOVATreatment$FDR < cutFDR] <- "Up"
# identify significantly down DE genes
tagsTblANOVATreatment$colorDE[tagsTblANOVATreatment$logFC < (-1*cutLFC) & tagsTblANOVATreatment$FDR < cutFDR] <- plotColors[6]
tagsTblANOVATreatment$alphaDE[tagsTblANOVATreatment$logFC < (-1*cutLFC) & tagsTblANOVATreatment$FDR < cutFDR] <- 1
#tagsTblANOVATreatment$colorDE[sign(tagsTblANOVATreatment$logFC) == -1 & tagsTblANOVATreatment$FDR < cutFDR] <- "Down"
# add column with -log10(FDR) values
tagsTblANOVATreatment$negLog10FDR <- -log10(tagsTblANOVATreatment$FDR)
# create volcano plot
jpeg("glmQLF_2WayANOVA_treatment_volcano.jpg")
ggplot(data=tagsTblANOVATreatment, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
# create volcano plot with labels
#labelSetTreatment <- tagsTblANOVATreatment[tagsTblANOVATreatment$colorDE == "UP" | tagsTblANOVATreatment$colorDE == "DOWN",]
#jpeg("glmQLF_2WayANOVA_treatment_volcanoLabeled.jpg")
#ggplot(data=tagsTblANOVATreatment, aes(x=logFC, y=-log10(FDR), color = colorDE)) + 
#  geom_point() +
#  ggrepel::geom_text_repel(data = labelSetTreatment, aes(label = row.names(labelSetTreatment))) +
#  theme_minimal() +
#  scale_colour_discrete(type = plotColorSubset, breaks = c("Up", "Down")) +
#  xlab("LFC")
#dev.off()

# add column for identifying direction of DE gene expression
tagsTblANOVATolerance$colorDE <- plotColors[5]
tagsTblANOVATolerance$alphaDE <- 0.75
# identify significantly up DE genes
tagsTblANOVATolerance$colorDE[tagsTblANOVATolerance$logFC > cutLFC & tagsTblANOVATolerance$FDR < cutFDR] <- plotColors[4]
tagsTblANOVATolerance$alphaDE[tagsTblANOVATolerance$logFC > cutLFC & tagsTblANOVATolerance$FDR < cutFDR] <- 1
#tagsTblANOVATolerance$colorDE[sign(tagsTblANOVATolerance$logFC) == 1 & tagsTblANOVATolerance$FDR < cutFDR] <- "Up"
# identify significantly down DE genes
tagsTblANOVATolerance$colorDE[tagsTblANOVATolerance$logFC < (-1*cutLFC) & tagsTblANOVATolerance$FDR < cutFDR] <- plotColors[6]
tagsTblANOVATolerance$alphaDE[tagsTblANOVATolerance$logFC < (-1*cutLFC) & tagsTblANOVATolerance$FDR < cutFDR] <- 1
#tagsTblANOVATolerance$colorDE[sign(tagsTblANOVATolerance$logFC) == -1 & tagsTblANOVATolerance$FDR < cutFDR] <- "Down"
# add column with -log10(FDR) values
tagsTblANOVATolerance$negLog10FDR <- -log10(tagsTblANOVATolerance$FDR)
# create volcano plot
jpeg("glmQLF_2WayANOVA_tolerance_volcano.jpg")
ggplot(data=tagsTblANOVATolerance, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
# create volcano plot with labels
#labelSetTolerance <- tagsTblANOVATolerance[tagsTblANOVATolerance$colorDE == "UP" | tagsTblANOVATolerance$colorDE == "DOWN",]
#jpeg("glmQLF_2WayANOVA_tolerance_volcanoLabeled.jpg")
#ggplot(data=tagsTblANOVATolerance, aes(x=logFC, y=-log10(FDR), color = colorDE)) + 
#  geom_point() +
#  ggrepel::geom_text_repel(data = labelSetTolerance, aes(label = row.names(labelSetTolerance)), max.overlaps=20) +
#  theme_minimal() +
#  scale_colour_discrete(type = plotColorSubset, breaks = c("Up", "Down")) +
#  xlab("LFC")
#dev.off()

# add column for identifying direction of DE gene expression
tagsTblANOVAInter$colorDE <- plotColors[5]
tagsTblANOVAInter$alphaDE <- 0.75
# identify significantly up DE genes
tagsTblANOVAInter$colorDE[tagsTblANOVAInter$logFC > cutLFC & tagsTblANOVAInter$FDR < cutFDR] <- plotColors[4]
tagsTblANOVAInter$alphaDE[tagsTblANOVAInter$logFC > cutLFC & tagsTblANOVAInter$FDR < cutFDR] <- 1
#tagsTblANOVAInter$colorDE[sign(tagsTblANOVAInter$logFC) == 1 & tagsTblANOVAInter$FDR < cutFDR] <- "Up"
# identify significantly down DE genes
tagsTblANOVAInter$colorDE[tagsTblANOVAInter$logFC < (-1*cutLFC) & tagsTblANOVAInter$FDR < cutFDR] <- plotColors[6]
tagsTblANOVAInter$alphaDE[tagsTblANOVAInter$logFC < (-1*cutLFC) & tagsTblANOVAInter$FDR < cutFDR] <- 1
#tagsTblANOVAInter$colorDE[sign(tagsTblANOVAInter$logFC) == -1 & tagsTblANOVAInter$FDR < cutFDR] <- "Down"
# add column with -log10(FDR) values
tagsTblANOVAInter$negLog10FDR <- -log10(tagsTblANOVAInter$FDR)
# create volcano plot
jpeg("glmQLF_2WayANOVA_interaction_volcano.jpg")
ggplot(data=tagsTblANOVAInter, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
# create volcano plot with labels
#labelSetInteraction <- tagsTblANOVAInter[tagsTblANOVAInter$colorDE == "UP" | tagsTblANOVAInter$colorDE == "DOWN",]
#jpeg("glmQLF_2WayANOVA_interaction_volcanoLabeled.jpg")
#ggplot(data=tagsTblANOVAInter, aes(x=logFC, y=-log10(FDR), color = colorDE)) + 
#  geom_point() +
#  ggrepel::geom_text_repel(data = labelSetInteraction, aes(label = row.names(labelSetInteraction)), max.overlaps=100) +
#  theme_minimal() +
#  scale_colour_discrete(type = plotColorSubset, breaks = c("Up", "Down")) +
#  xlab("LFC")
#dev.off()


# venn diagram
# identify significantly DE genes by FDR
tagsTblANOVATreatment.glm_keep <- tagsTblANOVATreatment$FDR < cutFDR
# create filtered results table of DE genes
tagsTblANOVATreatment.filtered <- tagsTblANOVATreatment[tagsTblANOVATreatment.glm_keep,]
# identify significantly DE genes by FDR
tagsTblANOVATolerance.glm_keep <- tagsTblANOVATolerance$FDR < cutFDR
# create filtered results table of DE genes
tagsTblANOVATolerance.filtered <- tagsTblANOVATolerance[tagsTblANOVATolerance.glm_keep,]
# identify significantly DE genes by FDR
tagsTblANOVAInter.glm_keep <- tagsTblANOVAInter$FDR < cutFDR
# create filtered results table of DE genes
tagsTblANOVAInter.filtered <- tagsTblANOVAInter[tagsTblANOVAInter.glm_keep,]
# retrieve set of DE gene names for hours contrast
geneSet_treatment <- rownames(tagsTblANOVATreatment.filtered)
# retrieve set of DE gene names for hours contrast
geneSet_tolerance <- rownames(tagsTblANOVATolerance.filtered)
# retrieve set of DE gene names for interaction contrast
geneSet_interaction <- rownames(tagsTblANOVAInter.filtered)
# create combined glm_list of DE gene names
glm_list_venn <- list(treatment = geneSet_treatment, 
                      tolerance = geneSet_tolerance,
                      interaction = geneSet_interaction)
# create venn diagram
jpeg("glmQLF_2WayANOVA_venn.jpg")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance","Interaction")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()
# create venn lists
#vennList <- venn(glm_list_venn, show.plot = FALSE)
# retrieve intersections
#listaAtt <- attributes(vennList)$intersections


# heatmaps
# heatmap data
logcounts = cpm(list, log=TRUE)
# view DGE genes
# subset counts table by DE gene set
DGESubset_treatment <- tagsTblANOVATreatment.filtered[!grepl(plotColors[5], tagsTblANOVATreatment.filtered$colorDE),]
# logcountsSubset_treatment <- subset(logcounts,
#                           grepl(
#                             paste0(rownames(DGESubset_treatment), collapse = "|"),
#                             rownames(logcounts),
#                             ignore.case = TRUE
#                           )
# )
DGESubset_treatment.keep <- rownames(logcounts) %in% DGESubset_treatment$gene
logcountsSubset_treatment <- logcounts[DGESubset_treatment.keep, ]
jpeg("glmQLF_treatment_heatmap.jpg")
heatmap(logcountsSubset_treatment, margins = c(8, 1), labRow = FALSE)
dev.off()


# view tolerance DGE genes
# subset counts table by DE gene set
DGESubset_tolerance <- tagsTblANOVATolerance.filtered[!grepl(plotColors[5], tagsTblANOVATolerance.filtered$colorDE),]
# logcountsSubset_tolerance <- subset(logcounts,
#                           grepl(
#                             paste0(rownames(DGESubset_tolerance), collapse = "|"),
#                             rownames(logcounts),
#                             ignore.case = TRUE
#                           )
# )
DGESubset_tolerance.keep <- rownames(logcounts) %in% DGESubset_tolerance$gene
logcountsSubset_tolerance <- logcounts[DGESubset_tolerance.keep, ]
jpeg("glmQLF_tolerance_heatmap.jpg")
heatmap(logcountsSubset_tolerance, margins = c(8, 1), labRow = FALSE)
dev.off()

# view interaction DGE genes
# subset counts table by DE gene set
DGESubset_interaction <- tagsTblANOVAInter.filtered[!grepl(plotColors[5], tagsTblANOVAInter.filtered$colorDE),]
# logcountsSubset_interaction <- subset(logcounts,
#                           grepl(
#                             paste0(rownames(DGESubset_interaction), collapse = "|"),
#                             rownames(logcounts),
#                             ignore.case = TRUE
#                           )
# )
DGESubset_interaction.keep <- rownames(logcounts) %in% DGESubset_interaction$gene
logcountsSubset_interaction <- logcounts[DGESubset_interaction.keep, ]
jpeg("glmQLF_interaction_heatmap.jpg")
heatmap(logcountsSubset_interaction, margins = c(8, 1), labRow = FALSE)
dev.off()
