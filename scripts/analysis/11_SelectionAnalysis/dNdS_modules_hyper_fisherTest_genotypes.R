#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(ghibli)

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

#GO enrichment
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


# Plotting Palettes
# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

# retrieve dN dS values
positiveTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove NAs
positiveTable <- na.omit(positiveTable)

# subset the positively selected genes
positiveSubset <- positiveTable[positiveTable$dNdS > 1 & positiveTable$dNdS < 99,]

# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# add geneID column
geneID <- row.names(positiveSubset)
positiveSubset <- cbind(geneID,positiveSubset)
colnames(resultsTable)[1] ="geneID"

# full outer join data frames
resultsSubset <- merge(x = resultsTable, y = positiveSubset, 
                      by = "geneID", all=TRUE)

# remove rows with NAs
resultsSubset <- na.omit(resultsSubset)


# hypergeometric distribution (fishers test)
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/

# loop over each module color
pValues <- data.frame(
  color = unique(moduleColors),
  enrichment = rep(99,numMods),
  depletion = rep(99,numMods)
)
total <- nrow(resultsTable)
positive <- nrow(positiveSubset)
for(i in 1:numMods){
  # initialize variables
  de <- nrow(resultsTable[resultsTable$number == i,])
  overlap <- nrow(resultsSubset[resultsSubset$number == i,])
  
  # test for over-representation (enrichment)
  pValues$enrichment[i] <- phyper(overlap-1, de, total-de, positive,lower.tail= FALSE)
  
  # test for under-representation (depletion)
  pValues$depletion[i] <- phyper(overlap, de, total-de, positive,lower.tail= TRUE)
}

# write results to a csv file
write.csv(pValues, "fisherTest_positiveSelection_modules.csv", row.names=FALSE)

