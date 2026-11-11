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

#Set working directory
#workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/OLYM_30_modules"
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA_9Nov2026/Genotypes/OLYM_30_modules"
dir.create(workingDir)
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

# import positively selected gene IDs
inFile <- paste(posDir, "Pulex_Olympics_kaksResults_dNdS_positive_geneIDs.csv", sep="/")
pos_list <- readLines(inFile)
pos_list <- pos_list[-1]
pos_normList <- normList[rownames(normList) %in% pos_list,]

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
  # define the data range
  data_min <- min(normData_module)
  data_max <- max(normData_module)
  # create a symmetric break vector centered at 0
  n_colors <- 100
  data_colors <- colorRampPalette(c(plotColors[16], "white", plotColors[2]))(n_colors)
  # ensure that 0 is white
  break_length <- n_colors + 1
  data_breaks <- c(
    seq(data_min, 0, length.out = ceiling(n_colors / 2) + 1),
    seq(data_max / n_colors, data_max, length.out = floor(n_colors / 2))
  )
  # create heatmap of normalized log expression
  as.ggplot(pheatmap(normData_module, 
                     annotation_col = exp_factor, 
                     annotation_colors = exp_colour,
                     #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     color = data_colors,
                     breaks = data_breaks
                     )
            )
  # save plot to a file
  ggsave(paste(unique(moduleColors)[modColor], "pheatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  # create scaled heatmap of normalized log expression
  as.ggplot(pheatmap(normData_module, 
                     scale="row",
                     annotation_col = exp_factor, 
                     annotation_colors = exp_colour,
                     #main = paste(unique(moduleColors)[modColor], "Heatmap of Samples"), 
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     color = data_colors
    )
  )
  # save plot to a file
  ggsave(paste(unique(moduleColors)[modColor], "scaled_pheatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  
}

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
                       color = colorRampPalette(c(plotColors[16], "white", plotColors[2]))(100)
    )
    )
    # save plot to a file
    ggsave(paste(unique(moduleColors)[modColor], "positive_heatmap.png", sep = "_"), bg = "white", width = 30, height = 20, units = "cm")
  }
}
