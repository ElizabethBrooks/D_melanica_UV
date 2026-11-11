#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)
library(rcartocolor)
library(ggVennDiagram)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[5], plotColors[6])

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
#workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_14Oct2026/interaction_nested_random_group_genotype/GOAnalysis"
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_9Nov2026/interaction_nested_random_group_genotype_LFC_0_FDR_0.05/GOAnalysis"

# set working directory
setwd(workingDir)

# read in data on significant GO terms (BP, MF, and CC) for each effect 
dream_BP_GO_terms <- read.csv('treatment_BP_sigGO_terms.csv', row.names = 1)
dream_MF_GO_terms <- read.csv('treatment_MF_sigGO_terms.csv', row.names = 1)
dream_CC_GO_terms <- read.csv('treatment_CC_sigGO_terms.csv', row.names = 1)

# retrieve previous results
edgeR_BP_GO_terms <- read.csv('/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis_fisher/treatment_BP_sigGO_terms.csv', row.names = 1)
edgeR_MF_GO_terms <- read.csv('/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis_fisher/treatment_MF_sigGO_terms.csv', row.names = 1)
edgeR_CC_GO_terms <- read.csv('/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis_fisher/treatment_CC_sigGO_terms.csv', row.names = 1)

# filter tables
dream_BP_GO_terms_filtered <- dream_BP_GO_terms[grepl("GO:", rownames(dream_BP_GO_terms)), ]
dream_MF_GO_terms_filtered <- dream_MF_GO_terms[grepl("GO:", rownames(dream_MF_GO_terms)), ]
dream_CC_GO_terms_filtered <- dream_CC_GO_terms[grepl("GO:", rownames(dream_CC_GO_terms)), ]
edgeR_BP_GO_terms_filtered <- edgeR_BP_GO_terms[grepl("GO:", rownames(edgeR_BP_GO_terms)), ]
edgeR_MF_GO_terms_filtered <- edgeR_MF_GO_terms[grepl("GO:", rownames(edgeR_MF_GO_terms)), ]
edgeR_CC_GO_terms_filtered <- edgeR_CC_GO_terms[grepl("GO:", rownames(edgeR_CC_GO_terms)), ]

# create combined list of DE gene names
BP_list_venn <- list(dream = rownames(dream_BP_GO_terms_filtered), edgeR = rownames(edgeR_BP_GO_terms_filtered))
MF_list_venn <- list(dream = rownames(dream_MF_GO_terms_filtered), edgeR = rownames(edgeR_MF_GO_terms_filtered))
CC_list_venn <- list(dream = rownames(dream_CC_GO_terms_filtered), edgeR = rownames(edgeR_CC_GO_terms_filtered))

# create venn diagrams

jpeg("treatment_BP_results_venn.jpg")
ggVennDiagram(BP_list_venn, label_alpha=0.25, category.names = c("dream","edgeR")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

jpeg("treatment_MF_results_venn.jpg")
ggVennDiagram(MF_list_venn, label_alpha=0.25, category.names = c("dream","edgeR")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

jpeg("treatment_CC_results_venn.jpg")
ggVennDiagram(CC_list_venn, label_alpha=0.25, category.names = c("dream","edgeR")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()
