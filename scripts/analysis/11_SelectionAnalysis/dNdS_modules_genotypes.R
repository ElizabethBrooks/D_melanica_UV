#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(rcartocolor)
library(tidyr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Genotypes"
setwd(workingDir)

# retrieve subset tag
set <- "OLYM"

# set the minimum module size
minModSize <- "30"

# retrieve WGCNA directory
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames1 = load(file = importFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
numMods <- length(unique(moduleColors))
colorTable <- data.frame(
  color = unique(moduleColors),
  number = seq(from = 1, to = numMods, by = 1)
)

# initialize module data frame
resultsTable <- data.frame(
  gene = character(),
  color = character(),
  number = numeric()
)

# match gene IDs with module colors
for(i in 1:numMods){
  gene <- names(datExpr)[moduleColors==colorTable[i,1]]
  color <- rep(colorTable[i,1], length(gene))
  number <- rep(colorTable[i,2], length(gene))
  moduleData <- cbind(gene, color, number)
  resultsTable <- rbind(resultsTable, moduleData)
}


# retrieve dN dS values
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove NAs
dNdSSubset <- na.omit(dNdSTable)

# subset dN dS values to remove outliers
#dNdSSubset <- dNdSSubset[dNdSSubset$dNdS < 99,]


# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# add geneID column
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
colnames(resultsTable)[1] ="geneID"

# full outer join data frames
resultsTable <- merge(x = resultsTable, y = dNdSSubset, 
                      by = "geneID", all=TRUE)

# set color and numner tags for NAs
resultsTable$color <- resultsTable$color %>% replace_na('None')
resultsTable$number <- resultsTable$number %>% replace_na('None')

# remove NAs
resultsTable <- na.omit(resultsTable)

# add column for identifying mode of selection
resultsTable$Selection[resultsTable$dNdS == 99] <- "Error"
resultsTable$Selection[resultsTable$dNdS > 1 & resultsTable$dNdS < 99] <- "Positive"
resultsTable$Selection[resultsTable$dNdS < 1] <- "Negative"

# create label set
labelSetNeg <- resultsTable[resultsTable$Selection == "Negative",]
labelSetPos <- resultsTable[resultsTable$Selection == "Positive",]
labelSetError <- resultsTable[resultsTable$Selection == "Error",]

# subset results
resultsSubset <- resultsTable[resultsTable$dNdS < 99,]
resultsSubset <- resultsSubset[resultsSubset$color != "None",]


# plotting

# box plot
jpeg("dNdS_modules_boxPlot.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# box plot colored by selection
jpeg("dNdS_modules_selection_boxPlot_error.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS, color=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_colour_discrete(type = plotColorSubset) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# box plot colored by selection
# without error values
jpeg("dNdS_modules_selection_boxPlot.jpg")
ggplot(resultsSubset, aes(x=color, y=dNdS, color=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_colour_discrete(type = c(plotColorSubset[2],plotColorSubset[3])) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# violin plots
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# violin plot color
jpeg("dNdS_modules_violinPlot.jpg")
ggplot(resultsSubset, aes(x=color, y=dNdS)) +
  theme_minimal() +
  geom_violin(trim=FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# violin plot color by selection
jpeg("dNdS_modules_selection_violinPlot.jpg")
ggplot(resultsSubset, aes(x=color, y=dNdS, color=Selection)) +
  theme_minimal() +
  geom_violin(trim=FALSE) +
  scale_colour_discrete(type = plotColorSubset) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# violin plot color by selection
jpeg("dNdS_modules_selection_violinPlot_error.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS, color=Selection)) +
  theme_minimal() +
  geom_violin(trim=FALSE) +
  scale_colour_discrete(type = plotColorSubset) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# scatter plots
# http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# https://stackoverflow.com/questions/15015356/how-to-do-selective-labeling-with-ggplot-geom-point
# http://www.sthda.com/english/wiki/ggplot2-point-shapes

# subset dN dS values to remove outliers
resultsTable <- resultsTable[resultsTable$dS < 10,]

# scatter plot colored by selection
# shaped by DEG effect set
jpeg("dNdS_modules_selection_scatterPlot_jittered.jpg")
ggplot(resultsTable, aes(x=dS, y=dN, shape=Selection, color=color)) +
  theme_minimal() +
  geom_point() +
  geom_jitter() +
  scale_shape_manual(values=c(2, 1, 3))
dev.off()

# scatter plot colored by selection
# shaped by DEG effect set
# with error genes labeled by geneID
#jpeg("dNdS_modules_selection_scatterPlot_labeled_error.jpg")
#ggplot(resultsTable, aes(x=dS, y=dN, shape=Selection, color=color)) +
#  theme_minimal() +
#  geom_point() +
#  ggrepel::geom_text_repel(data = labelSetError, aes(label = labelSetError$geneID), max.overlaps=100, color="Black") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  scale_shape_manual(values=c(2, 1, 3))
#dev.off()
