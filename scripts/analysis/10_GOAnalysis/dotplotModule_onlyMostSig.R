#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[5], plotColors[6])

# turn off scientific notation
options(scipen = 999)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
#workingDir <- args[1]
#workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_ks_OLYM_30"
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA_9Nov2026/Genotypes/GOAnalysis_ks_OLYM_30"

# set working directory
setwd(workingDir)

# retrieve subset tag
#set <- args[2]
set <- "OLYM"

# set the minimum module size
#minModSize <- args[3]
minModSize <- "30"

# retrieve WGCNA directory
#inDir <- args[4]
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# import positively selected gene set p-values
#positiveTable <- read.csv(file=args[5], row.names="color")
positiveTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/Genotypes/fisherTest_positiveSelection_modules.csv", row.names="color")

# import effect ANOVA p-values
#effectTable <- read.csv(file=args[6], row.names="module")
effectTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/DEGsANOVA_OLYM_30/aov_summary_pValues.csv", row.names="module")

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
num_mods <- length(unique(moduleColors)) + 1
color_list <- c(unique(moduleColors), "None")

# create data frame to hold plotting data for each module
module_BP_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_BP_results) <- c("Color","Term","Significant","weightFisher")
module_MF_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_MF_results) <- c("Color","Term","Significant","weightFisher")
module_CC_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_CC_results) <- c("Color","Term","Significant","weightFisher")

# loop through each module color
for(j in 1:num_mods){
  # read in data on significant GO terms (BP, MF, and CC) 
  BP_GO_sig <- read.csv(paste(color_list[j], 'BP_sigGO_terms.csv', sep="_"), row.names = NULL)[1,]
  MF_GO_sig <- read.csv(paste(color_list[j], 'MF_sigGO_terms.csv', sep="_"), row.names = NULL)[1,]
  CC_GO_sig <- read.csv(paste(color_list[j], 'CC_sigGO_terms.csv', sep="_"), row.names = NULL)[1,]
  
  # add a column labeling the color to each GO term set
  BP_GO_sig <- cbind("Color" = color_list[j], BP_GO_sig[,2], BP_GO_sig[,4], BP_GO_sig[,6])
  MF_GO_sig <- cbind("Color" = color_list[j], MF_GO_sig[,2], MF_GO_sig[,4], MF_GO_sig[,6])
  CC_GO_sig <- cbind("Color" = color_list[j], CC_GO_sig[,2], CC_GO_sig[,4], CC_GO_sig[,6])
  
  # create tables of all GO terms for each level
  #list_all_BP_GO_included <- unique(c(BP_GO_sig$GO.ID))
  #list_all_MF_GO_included <- unique(c(MF_GO_sig$GO.ID))
  #list_all_CC_GO_included <- unique(c(CC_GO_sig$GO.ID))
  
  # combine all tables
  module_BP_results <- rbind(module_BP_results, BP_GO_sig)
  module_MF_results <- rbind(module_MF_results, MF_GO_sig)
  module_CC_results <- rbind(module_CC_results, CC_GO_sig)
}

# re-name columns
colnames(module_BP_results) <- c("Color","Term","Significant","weightFisher")
colnames(module_MF_results) <- c("Color","Term","Significant","weightFisher")
colnames(module_CC_results) <- c("Color","Term","Significant","weightFisher")

# add GO level tags
module_BP_results <- cbind('Level' = 'BP', module_BP_results)
module_MF_results <- cbind('Level' = 'MF', module_MF_results)
module_CC_results <- cbind('Level' = 'CC', module_CC_results)

# combine GO level results for plotting
plotTable <- rbind(module_BP_results, module_MF_results, module_CC_results)

# remove NAs
plotTable <- na.omit(plotTable)

# remove < signs
plotTable$weightFisher <- gsub("<","",as.character(plotTable$weightFisher))

# loop through each module color
#plotTable$enrichment <- "NA"
#plotTable$depletion <- "NA"
#for(k in 1:num_mods){
  # add enrichment and depletion values
#  plotTable[plotTable$Color == row.names(positiveTable)[k],]$enrichment <- positiveTable$enrichment[k]
#  plotTable[plotTable$Color == row.names(positiveTable)[k],]$depletion <- positiveTable$depletion[k]
#}

# add signifigance tags
#plotTable$selection <- "#000000"
#plotTable[plotTable$enrichment <= 0.05,]$selection <- plotColors[9]
#plotTable[plotTable$depletion <= 0.05,]$selection <- plotColors[4]

# clean up row names
row.names(effectTable) <- gsub("ME","",as.character(row.names(effectTable)))

# loop through each module color
plotTable$treatment <- "NA"
plotTable$tolerance <- "NA"
#plotTable$interaction <- "NA"
for(k in 1:num_mods-1){
  # add enrichment and depletion values
  plotTable[plotTable$Color == row.names(effectTable)[k],]$treatment <- effectTable$treatment[k]
  plotTable[plotTable$Color == row.names(effectTable)[k],]$tolerance <- effectTable$tolerance[k]
  #plotTable[plotTable$Color == row.names(effectTable)[k],]$interaction <- effectTable$interaction[k]
}

# add signifigance tags
plotTable$effect <- "None"
plotTable[plotTable$treatment <= 0.05,]$effect <- "Treatment"
plotTable[plotTable$tolerance <= 0.05,]$effect <- "Tolerance"
#plotTable[plotTable$interaction <= 0.05,]$effect <- "Interaction"
plotTable[plotTable$treatment <= 0.05 & plotTable$tolerance <= 0.05,]$effect <- "Both"

# setup facet groups
x_axis_order <- factor(plotTable$Color, levels = color_list)
facetLevel <- factor(plotTable$Level, levels = c('BP', 'CC', 'MF'))
#facetEffect <- factor(plotTable$effect, levels = c('Treatment', 'Tolerance', 'Interaction'))

# create dot plot of significant GO terms
dotplot <- ggplot(data = plotTable, aes(x = x_axis_order, y = Term, size = as.numeric(Significant), color = as.numeric(weightFisher))) + 
  facet_grid(rows = facetLevel, space = 'free_y', scales = 'free') +
  geom_point() +
  #scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  scale_color_gradientn(colors = plotColorSubset) +
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, color = plotTable$selection)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_text(position = position_dodge(width = 1), aes(x=effect, y=0)) +
  xlab('Color') +
  ylab('GO Term') + 
  #scale_x_discrete(labels=c("sienna3"=expression(bold("sienna3")), "plum1"=expression(bold("plum1")), "darkgrey"=expression(italic("darkgrey")), "darkorange2"=expression(italic("darkorange2")), "floralwhite"=expression(italic("floralwhite")), "royalblue"=expression(italic("royalblue")), "salmon4"=expression(italic("salmon4")), parse=TRUE)) +
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplotModule_mostSigGO.pdf', plot = dotplot, device = 'pdf')


# remove rows not associated with an effect
plotSubset <- plotTable[plotTable$effect != "None",]

# setup facet groups
x_axis_order <- factor(plotSubset$Color, levels = color_list)
facetLevel <- factor(plotSubset$Level, levels = c('BP', 'CC', 'MF'))
#facetEffect <- factor(plotSubset$effect, levels = c('Treatment', 'Tolerance', 'Interaction'))

# create dot plot of significant GO terms
dotplot <- ggplot(data = plotSubset, aes(x = x_axis_order, y = Term, size = as.numeric(Significant), color = as.numeric(weightFisher))) + 
  #facet_grid(rows = facetLevel, space = 'free_y', scales = 'free') +
  #facet_grid(. ~ Level + effect, scales = "free_x", space = "free") +
  #facet_grid(Level + effect ~ ., scales = "free_y", space = "free") +
  facet_grid(Level ~ effect, scales = "free", space = "free") +
  geom_point() +
  #scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  scale_color_gradientn(colors = plotColorSubset) +
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, color = plotSubset$selection)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_text(position = position_dodge(width = 1), aes(x=effect, y=0)) +
  xlab('Color') +
  ylab('GO Term') + 
  scale_x_discrete(labels=c("sienna3"=expression(bold("sienna3")), "plum1"=expression(bold("plum1")), "darkgrey"=expression(italic("darkgrey")), "darkorange2"=expression(italic("darkorange2")), "floralwhite"=expression(italic("floralwhite")), "royalblue"=expression(italic("royalblue")), "salmon4"=expression(italic("salmon4")), parse=TRUE)) +
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplotModule_mostSigGO_faceted.pdf', plot = dotplot, device = 'pdf')


# select the salmon4 module
plotModules <- plotTable[plotTable$Color %in% c("lightyellow","salmon4","floralwhite","skyblue"),]

# setup facet groups
x_axis_order <- factor(plotModules$Color, levels = color_list)
facetLevel <- factor(plotModules$Level, levels = c('BP', 'CC', 'MF'))

# create dot plot of significant GO terms
dotplot <- ggplot(data = plotModules, aes(x = x_axis_order, y = Term, size = as.numeric(Significant), color = as.numeric(weightFisher))) + 
  facet_grid(Level ~ effect, scales = "free", space = "free") +
  geom_point() +
  scale_color_gradientn(colors = plotColorSubset) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Color') +
  ylab('GO Term') + 
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplotModule_mostSigGO_faceted.pdf', plot = dotplot, device = 'pdf')
