#!/usr/bin/env Rscript

#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Load the libraries
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggplotify)
library(rcartocolor)

# color blind safe plotting palettes
plotColors <- carto_pal(12, "Safe")

#Set working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/OLYM_30_modules"
setwd(workingDir)

# set networks directory
netDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set counts directory
deDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes"

# set positively selected genes directory
posDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests"

# retrieve subsetTag tag
set <- "OLYM"

# set the minimum module size
minModSize <- 30

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
inFile <- paste(netDir, importFile, sep="/")
lnames1 = load(file = inFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
inFile <- paste(netDir, importFile, sep="/")
lnames2 = load(file = inFile)

# import normalized gene count data for the Olympics
inFile <- paste(deDir, "glmQLF_normalizedCounts_logTransformed.csv", sep="/")
#normList <- read.csv(file=inFile, row.names="gene")[ ,1:24]
normList <- read.csv(file=inFile, row.names="gene")
# names(normList) <- c("HT_K_VIS_Pool1", "HT_K_VIS_Pool2", "HT_K_VIS_Pool3",
#                      "HT_K_UV_Pool1", "HT_K_UV_Pool2", "HT_K_UV_Pool3",
#                      "HT_I_VIS_Pool1", "HT_I_VIS_Pool2", "HT_I_VIS_Pool3",
#                      "HT_I_UV_Pool1", "HT_I_UV_Pool2", "HT_I_UV_Pool3",
#                      "LT_C_VIS_Pool1", "LT_C_VIS_Pool2", "LT_C_VIS_Pool3",
#                      "LT_C_UV_Pool1", "LT_C_UV_Pool2", "LT_C_UV_Pool3",
#                      "LT_A_VIS_Pool1", "LT_A_VIS_Pool2", "LT_A_VIS_Pool3",
#                      "LT_A_UV_Pool1", "LT_A_UV_Pool2", "LT_A_UV_Pool3")
names(normList) <- c("HT_K_VIS_Pool1", "HT_K_VIS_Pool2", "HT_K_VIS_Pool3",
                     "HT_K_UV_Pool1", "HT_K_UV_Pool2", "HT_K_UV_Pool3",
                     "HT_I_VIS_Pool1", "HT_I_VIS_Pool2", "HT_I_VIS_Pool3",
                     "HT_I_UV_Pool1", "HT_I_UV_Pool2", "HT_I_UV_Pool3",
                     "LT_C_VIS_Pool1", "LT_C_VIS_Pool2", "LT_C_VIS_Pool3",
                     "LT_C_UV_Pool1", "LT_C_UV_Pool2", "LT_C_UV_Pool3",
                     "LT_A_VIS_Pool1", "LT_A_VIS_Pool2", "LT_A_VIS_Pool3",
                     "LT_A_UV_Pool1", "LT_A_UV_Pool2", "LT_A_UV_Pool3")

# import DE gene IDs
inFile <- paste(deDir, "glmQLF_2WayANOVA_DEGs_geneIDs_unique.csv", sep="/")
DEGs_list <- readLines(inFile)
DEGs_list <- DEGs_list[-1]
DEGs_normList <- normList[rownames(normList) %in% DEGs_list,]

# import treatment DE gene IDs
inFile <- paste(deDir, "glmQLF_2WayANOVA_treatment_DEGs_geneIDs.csv", sep="/")
DEGs_treatment_list <- readLines(inFile)
DEGs_treatment_list <- DEGs_treatment_list[-1]
DEGs_treatment_normList <- normList[rownames(normList) %in% DEGs_treatment_list,]

# import tolerance DE gene IDs
inFile <- paste(deDir, "glmQLF_2WayANOVA_tolerance_DEGs_geneIDs.csv", sep="/")
DEGs_tolerance_list <- readLines(inFile)
DEGs_tolerance_list <- DEGs_tolerance_list[-1]
DEGs_tolerance_normList <- normList[rownames(normList) %in% DEGs_tolerance_list,]

# import interaction DE gene IDs
inFile <- paste(deDir, "glmQLF_2WayANOVA_interaction_DEGs_geneIDs.csv", sep="/")
DEGs_interaction_list <- readLines(inFile)
DEGs_interaction_list <- DEGs_interaction_list[-1]
DEGs_interaction_normList <- normList[rownames(normList) %in% DEGs_interaction_list,]

# import positively selected gene IDs
inFile <- paste(posDir, "Pulex_Olympics_kaksResults_dNdS_positive_geneIDs.csv", sep="/")
pos_list <- readLines(inFile)
pos_list <- pos_list[-1]
pos_normList <- normList[rownames(normList) %in% pos_list,]

# find positively selected genes associated with each effect
pos_treatment_normList <- pos_normList[rownames(pos_normList) %in% DEGs_treatment_list,]
pos_tolerance_normList <- pos_normList[rownames(pos_normList) %in% DEGs_tolerance_list,]
pos_interaction_normList <- pos_normList[rownames(pos_normList) %in% DEGs_interaction_list,]
pos_tolerance_interaction_normList <- pos_tolerance_normList[rownames(pos_tolerance_normList) %in% DEGs_interaction_list,]
pos_DEGs_normList <- pos_normList[rownames(pos_normList) %in% DEGs_list,]
# output the lists of gene IDs
write(c("geneID", rownames(pos_treatment_normList)), "treatment_positive_geneIDs.csv")
write(c("geneID", rownames(pos_tolerance_normList)), "tolerance_positive_geneIDs.csv")
write(c("geneID", rownames(pos_interaction_normList)), "interaction_positive_geneIDs.csv")
write(c("geneID", rownames(pos_tolerance_interaction_normList)), "tolerance_interaction_positive_geneIDs.csv")
write(c("geneID", rownames(pos_DEGs_normList)), "DEGs_positive_geneIDs.csv")

# identify genes in both the interaction and treament that are positive
duprows <- rownames(pos_interaction_normList) %in% rownames(pos_treatment_normList)
pos_treatment_interaction_normList <- rbind(pos_treatment_normList, pos_interaction_normList[!duprows,])
write(c("geneID", rownames(pos_treatment_interaction_normList)), "treatment_interaction_positive_geneIDs.csv")

# loop over each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- c("geneID", names(datExpr)[moduleColors==unique(moduleColors)[modColor]])
  # output the list of gene IDs
  write(geneID_list, paste(unique(moduleColors)[modColor], "geneIDs.csv", sep="_"))
}

#Create data frame with the experimental design layout
exp_factor <- data.frame(sample = gsub('.{6}$', '', names(normList)))
rownames(exp_factor) <- names(normList)
exp_colour = list(
  sample = c(HT_K_VIS = "#E69F00", HT_K_UV = "#D55E00", HT_I_VIS = "#CC79A7", HT_I_UV = plotColors[6], 
             LT_C_VIS = plotColors[7], LT_C_UV = plotColors[4], LT_A_VIS = "#56B4E9", LT_A_UV = plotColors[5])
)

# create heatmaps for each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- names(datExpr)[moduleColors==unique(moduleColors)[modColor]]
  # retrieve expression data
  normData_module <- normList[rownames(normList) %in% geneID_list,]
  # create heatmap of normalized log expression
  as.ggplot(pheatmap(normData_module, annotation_col = exp_factor, annotation_colors = exp_colour,
                     #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
                     )
            )
  # save plot to a file
  ggsave(paste(unique(moduleColors)[modColor], "heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
}

# create heatmaps for DE genes in each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- names(datExpr)[moduleColors==unique(moduleColors)[modColor]]
  # retrieve expression data
  normData_module <- DEGs_normList[rownames(DEGs_normList) %in% geneID_list,]
  if (nrow(normData_module) >= 2) {
    # create heatmap of normalized log expression
    as.ggplot(pheatmap(normData_module, annotation_col = exp_factor, annotation_colors = exp_colour,
                       #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
    )
    )
    # save plot to a file
    ggsave(paste(unique(moduleColors)[modColor], "DEGs_heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  }
}

# create heatmaps for treatment DE genes in each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- names(datExpr)[moduleColors==unique(moduleColors)[modColor]]
  # retrieve expression data
  normData_module <- DEGs_treatment_normList[rownames(DEGs_treatment_normList) %in% geneID_list,]
  if (nrow(normData_module) >= 2) {
    # create heatmap of normalized log expression
    as.ggplot(pheatmap(normData_module, annotation_col = exp_factor, annotation_colors = exp_colour,
                       #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
    )
    )
    # save plot to a file
    ggsave(paste(unique(moduleColors)[modColor], "treatment_DEGs_heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  }
}

# create heatmaps for tolerance DE genes in each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- names(datExpr)[moduleColors==unique(moduleColors)[modColor]]
  # retrieve expression data
  normData_module <- DEGs_tolerance_normList[rownames(DEGs_tolerance_normList) %in% geneID_list,]
  if (nrow(normData_module) >= 2) {
    # create heatmap of normalized log expression
    as.ggplot(pheatmap(normData_module, annotation_col = exp_factor, annotation_colors = exp_colour,
                       #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
    )
    )
    # save plot to a file
    ggsave(paste(unique(moduleColors)[modColor], "tolerance_DEGs_heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  }
}

# create heatmaps for interaction DE genes in each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- names(datExpr)[moduleColors==unique(moduleColors)[modColor]]
  # retrieve expression data
  normData_module <- DEGs_interaction_normList[rownames(DEGs_interaction_normList) %in% geneID_list,]
  if (nrow(normData_module) >= 2) {
    # create heatmap of normalized log expression
    as.ggplot(pheatmap(normData_module, annotation_col = exp_factor, annotation_colors = exp_colour,
                       #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
    )
    )
    # save plot to a file
    ggsave(paste(unique(moduleColors)[modColor], "interaction_DEGs_heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  }
}

# create heatmap of treatment DE genes normalized log expression
as.ggplot(pheatmap(DEGs_treatment_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("treatment_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmap of tolerance DE genes normalized log expression
as.ggplot(pheatmap(DEGs_tolerance_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("tolerance_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmap of interaction DE genes normalized log expression
as.ggplot(pheatmap(DEGs_interaction_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("interaction_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmaps for interaction DE genes in each module
for (modColor in 1:length(unique(moduleColors))){
  # retrieve gene IDs
  geneID_list <- names(datExpr)[moduleColors==unique(moduleColors)[modColor]]
  # retrieve expression data
  normData_module <- pos_normList[rownames(pos_normList) %in% geneID_list,]
  if (nrow(normData_module) >= 2) {
    # create heatmap of normalized log expression
    as.ggplot(pheatmap(normData_module, annotation_col = exp_factor, annotation_colors = exp_colour,
                       #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
    )
    )
    # save plot to a file
    ggsave(paste(unique(moduleColors)[modColor], "positive_heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  }
}

# create heatmap of positive genes normalized log expression
as.ggplot(pheatmap(pos_treatment_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("treatment_positive_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmap of positive genes normalized log expression
as.ggplot(pheatmap(pos_tolerance_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("tolerance_positive_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmap of positive genes normalized log expression
as.ggplot(pheatmap(pos_interaction_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("interaction_positive_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmap of positive DE genes normalized log expression
as.ggplot(pheatmap(pos_DEGs_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("positive_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")

# create heatmap of positive interaction and treatment genes normalized log expression
as.ggplot(pheatmap(pos_treatment_interaction_normList, annotation_col = exp_factor, annotation_colors = exp_colour,
                   #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = colorRampPalette(c("#F0E442", plotColors[7], plotColors[5]))(100)
)
)
# save plot to a file
ggsave("positive_interaction_treatment_DEGs_heatmap.png", bg = "white", width = 30, height = 20, units = "cm")
