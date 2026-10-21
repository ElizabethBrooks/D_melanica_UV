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

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis_14Oct2026/interaction_nested_random_group_genotype/GOAnalysis"

# set working directory
setwd(workingDir)

# read in data on significant GO terms (BP, MF, and CC) for each effect 
BP_GO_terms <- read.csv('treatment_BP_sigGO_terms.csv', row.names = 1)
MF_GO_terms <- read.csv('treatment_MF_sigGO_terms.csv', row.names = 1)
CC_GO_terms <- read.csv('treatment_CC_sigGO_terms.csv', row.names = 1)

## treatment effect
# filter for top 5 significant terms
BP_GO_sig <- BP_GO_terms[1:5, ]
MF_GO_sig <- MF_GO_terms[1:5, ]
CC_GO_sig <- CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
BP_plot_table <- cbind("Effect" = 'Treatment', BP_GO_sig)
MF_plot_table <- cbind("Effect" = 'Treatment', MF_GO_sig)
CC_plot_table <- cbind("Effect" = 'Treatment', CC_GO_sig)

# create tables of all GO terms for each level
#list_all_BP_GO_included <- unique(c(interaction_BP_GO_sig$GO.ID, treatment_BP_GO_sig$GO.ID, BP_GO_sig$GO.ID))
#list_all_MF_GO_included <- unique(c(interaction_MF_GO_sig$GO.ID, treatment_MF_GO_sig$GO.ID, MF_GO_sig$GO.ID))
#list_all_CC_GO_included <- unique(c(interaction_CC_GO_sig$GO.ID, treatment_CC_GO_sig$GO.ID, CC_GO_sig$GO.ID))

# BP
all_BP_plot_table <- BP_plot_table
all_BP_plot_table <- cbind('GO_cat' = 'BP', all_BP_plot_table)
# MF
all_MF_plot_table <- MF_plot_table
all_MF_plot_table <- cbind('GO_cat' = 'MF', all_MF_plot_table)
# CC
all_CC_plot_table <- CC_plot_table
all_CC_plot_table <- cbind('GO_cat' = 'CC', all_CC_plot_table)

# combine all tables
all_plot_table <- rbind(all_BP_plot_table, all_MF_plot_table, all_CC_plot_table)

# remove NAs
plot_table <- na.omit(all_plot_table)

# create dot plot of significant GO terms
facet <- factor(plot_table$GO_cat, levels = c('BP', 'CC', 'MF'))

# create dot plot
dotplot <- ggplot(data = plot_table, aes(x = "Treatment", y = Term, size = Significant, color = as.numeric(weightFisher))) + 
  facet_wrap(facet, ncol = 1, nrow = length(facet), scales = 'free', strip.position="right") +
  geom_point() +
  scale_color_gradientn(colors = plotColorSubset) +
  theme_bw() +
  xlab('Score') +
  ylab('Term') + 
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplot_sigGO.pdf', plot = dotplot, device = 'pdf')


