#!/usr/bin/env Rscript

#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

# load libraries
library(ggpubr)
#library(multcomp)
library(rcartocolor)
library(stringr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c("#D55E00", plotColors[6], plotColors[4], plotColors[5])

#Set working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_23Jan2026/Genotypes"
setwd(workingDir)

# import normalized gene count data for the Olympics
inFile <- "glmQLF_normalizedCounts_logTransformed.csv"
inputList <- read.csv(file=inFile, row.names="gene")

# set input gene ID
gene_ID <- "gene-LOC124188748"
#gene_ID <- "gene-LOC124190117"
#gene_ID <- "gene-LOC124202665"
#gene_ID <- "gene-LOC124193197"

# subset data for PHR (LOC124188748)
normList <- inputList[gene_ID,]

# data frame for combined expression values
exp_data <- data.frame(
  expression = rep(NA, ncol(normList)),
  genotype = rep(NA, ncol(normList)),
  treatment = rep(NA, ncol(normList)),
  tolerance = rep(NA, ncol(normList))
)

# fill in the data frame of expression values
for (index in 1:ncol(normList)) {
  # add the expression data
  exp_data$expression[index] <- normList[,index]
  exp_data$genotype[index] <- paste(
    str_split_i(colnames(normList)[index], "_", 1),
    str_split_i(colnames(normList)[index], "_", 2),
    sep = "_")
  exp_data$treatment[index] <- str_split_i(colnames(normList)[index], "_", 3)
  exp_data$tolerance[index] <- str_sub(exp_data$genotype[index], end = -3)
}

# create a colored box plot for all genes
exp_data$tolerance <- factor(exp_data$tolerance, levels=c("LT", "HT"))
log_plot <- ggboxplot(data=exp_data, x="treatment", y="expression", color="genotype",
          palette = plotColorSubset) +
          facet_grid(~tolerance, scales = "free", space = "free")
exportFile <- paste("logExp_", gene_ID, ".png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(log_plot)
dev.off()

