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
library(pheatmap)
library(ggplotify)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
defaultColors <- palette.colors(palette = "Okabe-Ito")
blindColors <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safeColors <- carto_pal(12, "Safe")
plotColors <- c(safeColors, blindColors, defaultColors)
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6], plotColors[7])
colorList <- c(plotColors[1], plotColors[2])
# view colors
scales::show_col(plotColors)

# set working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_9Nov2026/interaction_nested_random_group_genotype_LFC_0_FDR_0.05"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv", row.names="gene")[ ,1:24]

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/Repos/D_melanica_UV/inputs/expDesign_OlympicsGenotypesGroup.csv", row.names="sample")

# update sample and gneotype IDs
colnames(inputTable) <- gsub("E05", "LT_C", colnames(inputTable))
colnames(inputTable) <- gsub("R2", "LT_A", colnames(inputTable))
colnames(inputTable) <- gsub("Y05", "HT_K", colnames(inputTable))
colnames(inputTable) <- gsub("Y023_5", "HT_I", colnames(inputTable))
rownames(targets) <- gsub("E05", "LT_C", rownames(targets))
rownames(targets) <- gsub("R2", "LT_A", rownames(targets))
rownames(targets) <- gsub("Y05", "HT_K", rownames(targets))
rownames(targets) <- gsub("Y023_5", "HT_I", rownames(targets))
targets$genotype <- gsub("E05", "LT_C", targets$genotype)
targets$genotype <- gsub("R2", "LT_A", targets$genotype)
targets$genotype <- gsub("Y05", "HT_K", targets$genotype)
targets$genotype <- gsub("Y023", "HT_I", targets$genotype)

# set input FDR cutoff
cutFDR <- 0.05
#cutFDR <- 0.1
#cutFDR <- 0.01
#cutFDR <- 0.025

# set LFC cut
#cutLFC <- 1.2
#cutLFC <- log2(1.2)
cutLFC <- 0

# trim the data table to remove summary statistics
countsTable <- head(inputTable, - 5)

# specify the reference level
targets$treatment <- relevel(factor(targets$treatment), ref = "VIS")
targets$group <- relevel(factor(targets$group), ref = "LT")
#targets$genotype <- relevel(factor(targets$genotype), ref = "R2")
targets$genotype <- relevel(factor(targets$genotype), ref = "LT_A")

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
treatment_dge_results <- topTable(fitmm, coef = "treatmentUV", number = nrow(dge))
tolerance_dge_results <- topTable(fitmm, coef = "groupHT", number = nrow(dge))
interaction_dge_results <- topTable(fitmm, coef = "treatmentUV:groupHT", number = nrow(dge))

# export table of DE genes
treatment_dge_results_tbl <- as_tibble(treatment_dge_results, rownames = "gene")
write.table(treatment_dge_results_tbl, file="treatmentUV.csv", sep=",", row.names=FALSE, quote=FALSE)
tolerance_dge_results_tbl <- as_tibble(tolerance_dge_results, rownames = "gene")
write.table(tolerance_dge_results_tbl, file="groupHT.csv", sep=",", row.names=FALSE, quote=FALSE)
interaction_dge_results_tbl <- as_tibble(interaction_dge_results, rownames = "gene")
write.table(interaction_dge_results_tbl, file="treatmentUV_groupHT.csv", sep=",", row.names=FALSE, quote=FALSE)

# get genes < FDR and LFC cutoffs
treatment_dge_sig <- topTable(fitmm, coef = "treatmentUV", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
tolerance_dge_sig <- topTable(fitmm, coef = "groupHT", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
interaction_dge_sig <- topTable(fitmm, coef = "treatmentUV:groupHT", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

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
interaction_edgeR <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes/glmQLF_2WayANOVA_interaction_topTags.csv")
positive <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/Pulex_Olympics_kaksResults_dNdS_positive.csv")

# identify significantly DE genes by FDR and LFC
treatment_edgeR_keep <- treatment_edgeR$FDR < cutFDR & (treatment_edgeR$logFC > cutLFC | treatment_edgeR$logFC < (-1*cutLFC))
group_edgeR_keep <- group_edgeR$FDR < cutFDR & (group_edgeR$logFC > cutLFC | group_edgeR$logFC < (-1*cutLFC))
interaction_edgeR_keep <- interaction_edgeR$FDR < cutFDR & (interaction_edgeR$logFC > cutLFC | interaction_edgeR$logFC < (-1*cutLFC))

# create filtered results table of DE genes
treatment_edgeR.filtered <- treatment_edgeR[treatment_edgeR_keep,]
tolerance_edgeR.filtered <- group_edgeR[group_edgeR_keep,]
interaction_edgeR.filtered <- interaction_edgeR[interaction_edgeR_keep,]

# output filtered results
treatment_out_file <- paste("treatment_edgeR", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_out_file <- paste(treatment_out_file, "csv", sep = ".")
write.table(treatment_edgeR.filtered, file=treatment_out_file, sep=",", row.names=FALSE, quote=FALSE)
tolerance_out_file <- paste("group_edgeR", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
tolerance_out_file <- paste(tolerance_out_file, "csv", sep = ".")
write.table(tolerance_edgeR.filtered, file=tolerance_out_file, sep=",", row.names=FALSE, quote=FALSE)
interaction_out_file <- paste("treatment_group_edgeR", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
interaction_out_file <- paste(interaction_out_file, "csv", sep = ".")
write.table(interaction_edgeR.filtered, file=interaction_out_file, sep=",", row.names=FALSE, quote=FALSE)

# list of DEGs
DEG_union <- union(treatment_dge_sig_tbl$gene, tolerance_dge_sig_tbl$gene)
dream_DEGs <- data.frame(geneID = DEG_union)
write.table(dream_DEGs, file="dream_results_DEGs.csv", sep=",", row.names=FALSE, quote=FALSE)

# import non-grey module gene IDs
module_genes <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/OLYM_30_modules/non_grey/geneIDs.csv", header = FALSE)

# create combined list of DE gene names
treatment_list_venn <- list(dream = treatment_dge_sig_tbl$gene, edgeR = treatment_edgeR.filtered$gene)
tolerance_list_venn <- list(dream = tolerance_dge_sig_tbl$gene, edgeR = tolerance_edgeR.filtered$gene)
interaction_list_venn <- list(dream = interaction_dge_sig_tbl$gene, edgeR = interaction_edgeR.filtered$gene)
dream_list_venn <- list(Treatment = treatment_dge_sig_tbl$gene, Tolerance = tolerance_dge_sig_tbl$gene)
edgeR_list_venn <- list(Treatment = treatment_edgeR.filtered$gene, Tolerance = tolerance_edgeR.filtered$gene, Interaction = interaction_edgeR.filtered$gene)
dream_positive_list_venn <- list(Treatment = treatment_dge_sig_tbl$gene, Tolerance = tolerance_dge_sig_tbl$gene, Positive = positive$geneID)
edgeR_positive_list_venn <- list(Treatment = treatment_edgeR.filtered$gene, Tolerance = tolerance_edgeR.filtered$gene, Interaction = interaction_edgeR.filtered$gene, Positive = positive$geneID)
dream_modules_positive_list_venn <- list(DEGs = dream_DEGs$geneID, Modules = module_genes$V1, Positive = positive$geneID)

# create venn diagrams

file <- "treatment_results_venn.png"
outVenn <- ggVennDiagram(treatment_list_venn, label_alpha=0.25, category.names = c("dream","edgeR")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "tolerance_results_venn.png"
outVenn <- ggVennDiagram(tolerance_list_venn, label_alpha=0.25, category.names = c("dream","edgeR")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "interaction_results_venn.png"
outVenn <- ggVennDiagram(interaction_list_venn, label_alpha=0.25, category.names = c("dream","edgeR")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "dream_results_venn.png"
outVenn <- ggVennDiagram(dream_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "edgeR_results_venn.png"
outVenn <- ggVennDiagram(edgeR_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance", "Interaction")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "dream_positive_venn.png"
outVenn <- ggVennDiagram(dream_positive_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance", "Positive")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "edgeR_positive_venn.png"
outVenn <- ggVennDiagram(edgeR_positive_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance", "Interaction", "Positive")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

file <- "dream_modules_positive_venn.png"
outVenn <- ggVennDiagram(dream_modules_positive_list_venn, label_alpha=0.25, category.names = c("DEGs", "Modules", "Positive")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradientn(colors = colorList)
ggsave(file, plot = outVenn, device = "png")

# get set intersections
edgeR_treatment_positive <- intersect(treatment_edgeR.filtered$gene, positive$geneID)
edgeR_tolerance_positive <- intersect(tolerance_edgeR.filtered$gene, positive$geneID)
edgeR_interaction_positive <- intersect(interaction_edgeR.filtered$gene, positive$geneID)
dream_treatment_positive <- intersect(treatment_dge_sig_tbl$gene, positive$geneID)

# output set intersections
write.table(edgeR_treatment_positive, file="edgeR_treatment_positive.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(edgeR_tolerance_positive, file="edgeR_tolerance_positive.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(edgeR_interaction_positive, file="edgeR_interaction_positive.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(dream_treatment_positive, file="dream_treatment_positive.csv", sep=",", row.names=FALSE, quote=FALSE)

# heatmap data
logcounts = cpm(list, normalized.lib.sizes=TRUE, log=TRUE)

# subset counts table by DE gene set
DGESubset_treatment <- treatment_dge_sig_tbl[!grepl(plotColors[5], treatment_dge_sig_tbl$colorDE),]
DGESubset_treatment.keep <- rownames(logcounts) %in% DGESubset_treatment$gene
logcountsSubset_treatment <- logcounts[DGESubset_treatment.keep, ]
DGESubset_tolerance <- tolerance_dge_sig_tbl[!grepl(plotColors[5], tolerance_dge_sig_tbl$colorDE),]
DGESubset_tolerance.keep <- rownames(logcounts) %in% DGESubset_tolerance$gene
logcountsSubset_tolerance <- logcounts[DGESubset_tolerance.keep, ]

# heatmap of treatment DEGs
jpeg("treatment_heatmap.jpg")
heatmap(logcountsSubset_treatment, margins = c(8, 1), labRow = FALSE)
dev.off()

# heatmap of tolerance DEGs
jpeg("tolerance_heatmap.jpg")
heatmap(logcountsSubset_tolerance, margins = c(8, 1), labRow = FALSE)
dev.off()

# update sample names for plotting


# sample plotting colors
exp_colour = list(
  Sample = c(HT_I_UV = "#CC79A7", HT_I_VIS = plotColors[6],
             HT_K_UV = "#E69F00", HT_K_VIS = "#D55E00",  
             LT_A_UV = "#56B4E9", LT_A_VIS = plotColors[5],
             LT_C_UV = plotColors[7], LT_C_VIS = plotColors[4]
  )
)

# create data frame with the experimental design layout
exp_factor <- data.frame(Sample = c(paste(targets$genotype, targets$treatment, sep="_")))
rownames(exp_factor) <- colnames(logcountsSubset_treatment)

# create a symmetric break vector centered at 0
n_colors <- 100
data_colors <- colorRampPalette(c(plotColors[16], "white", plotColors[2]))(n_colors)
break_length <- n_colors + 1

# Treatment Effect
# define the data range
data_min <- min(logcountsSubset_treatment)
data_max <- max(logcountsSubset_treatment)
# ensure that 0 is white
data_breaks <- c(
  seq(data_min, 0, length.out = ceiling(n_colors / 2) + 1),
  seq(data_max / n_colors, data_max, length.out = floor(n_colors / 2))
)
# create heatmap
as.ggplot(
  pheatmap(logcountsSubset_treatment, 
           annotation_col = exp_factor, 
           annotation_colors = exp_colour,
           #main="Heatmap of Treatment Effect DE Genes",
           #color = colorRampPalette(c(plotColors[16], "white", plotColors[2]))(100),
           color = data_colors,
           breaks = data_breaks,
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12
  )
)
ggsave("treatment_pheatmap.png", bg = "white", width = 30, height = 20, units = "cm")
# create scaled heatmap
as.ggplot(
  pheatmap(logcountsSubset_treatment, 
           scale="row", 
           annotation_col = exp_factor, 
           annotation_colors = exp_colour,
           #main="Heatmap of Treatment Effect DE Genes",
           #color = colorRampPalette(c(plotColors[16], "white", plotColors[2]))(100),
           color = data_colors,
           #breaks = data_breaks,
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12
  )
)
ggsave("treatment_scaled_pheatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# Tolerance Effect
# define the data range
data_min <- min(logcountsSubset_tolerance)
data_max <- max(logcountsSubset_tolerance)
# ensure that 0 is white
data_breaks <- c(
  seq(data_min, 0, length.out = ceiling(n_colors / 2) + 1),
  seq(data_max / n_colors, data_max, length.out = floor(n_colors / 2))
)
# create heatmap
as.ggplot(
  pheatmap(logcountsSubset_tolerance, 
           annotation_col = exp_factor,
           annotation_colors = exp_colour,
           #main="Heatmap of Tolerance Effect DE Genes", 
           #color = colorRampPalette(c(plotColors[16], "white", plotColors[2]))(100),
           color = data_colors,
           breaks = data_breaks,
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12
  )
)
ggsave("tolerance_pheatmap.png", bg = "white", width = 30, height = 20, units = "cm")
# create scaled heatmap
as.ggplot(
  pheatmap(logcountsSubset_tolerance, 
           scale="row", 
           annotation_col = exp_factor,
           annotation_colors = exp_colour,
           #main="Heatmap of Tolerance Effect DE Genes", 
           #color = colorRampPalette(c(plotColors[16], "white", plotColors[2]))(100),
           color = data_colors,
           #breaks = data_breaks,
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12
  )
)
ggsave("tolerance_scaled_pheatmap.png", bg = "white", width = 30, height = 20, units = "cm")

