#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(rcartocolor)
library(tidyr)

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Tolerance"
setwd(workingDir)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[11], plotColors[5], plotColors[6])

# retrieve gene lengths
lengthsTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/geneLengths.pep.csv", row.names="geneID")

# retrieve dN dS values
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove NAs
lengthSubset <- na.omit(lengthsTable)
dNdSSubset <- na.omit(dNdSTable)

# retrieve DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")

# keep only sig
interactionSig <- interactionTable[interactionTable$FDR < 0.05,]
treatmentSig <- treatmentTable[treatmentTable$FDR < 0.05,]
toleranceSig <- toleranceTable[toleranceTable$FDR < 0.05,]

# add effect tags
interactionSig$Effect <- "Interaction"
treatmentSig$Effect <- "Treatment"
toleranceSig$Effect <- "Tolerance"

# add geneID column
geneID <- row.names(lengthSubset)
lengthSubset <- cbind(geneID,lengthSubset)
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
geneID <- row.names(interactionSig)
interactionSig <- cbind(geneID,interactionSig)
geneID <- row.names(treatmentSig)
treatmentSig <- cbind(geneID,treatmentSig)
geneID <- row.names(toleranceSig)
toleranceSig <- cbind(geneID,toleranceSig)

# keep necessary columns
lengthSubset <- lengthSubset[,c("geneID","reference")]
dNdSSubset <- dNdSSubset[,c("geneID","dN","dS","dNdS")]
treatmentSubset <- treatmentSig[,c("geneID","Effect")]
toleranceSubset <- toleranceSig[,c("geneID","Effect")]
interactionSubset <- interactionSig[,c("geneID","Effect")]

# rename length column
colnames(lengthSubset)[colnames(lengthSubset) == "reference"] ="Length"

# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# combine all tables
effectTable <- rbind(interactionSubset, treatmentSubset, toleranceSubset)

# full outer join
geneTable <- merge(x = dNdSSubset, y = lengthSubset, 
                   by = "geneID", all=TRUE)

# full outer join
plotTable <- merge(x = geneTable, y = effectTable, 
                   by = "geneID", all=TRUE)

# set tags for NAs
# non-proteon coding
#plotTable$reference <- plotTable$Length %>% replace_na(-1)
#plotTable$dN <- plotTable$dN %>% replace_na(-1)
#plotTable$dS <- plotTable$dS %>% replace_na(-1)
#plotTable$dNdS <- plotTable$dNdS %>% replace_na(-1)
# non-DE
plotTable$Effect <- plotTable$Effect %>% replace_na('None')

# add selection type
plotTable$Selection[plotTable$dNdS == 99] <- "Error"
plotTable$Selection[plotTable$dNdS > 1 & plotTable$dNdS < 99] <- "Positive"
plotTable$Selection[plotTable$dNdS < 1] <- "Negative"

# remove NAs
plotTable <- na.omit(plotTable)

# set row order for plotting
plotTable$Order <- ifelse(plotTable$Effect=="None", 1, 2)
plotTable <- plotTable[order(plotTable$Order),]


# plotting

# dNdS vs Effect

# box plot colored by selection
jpeg("dNdS_DEGs_boxPlot.jpg")
ggplot(plotTable, aes(x=Effect, y=dNdS)) +
  geom_boxplot() +
  theme_minimal()
dev.off()

# box plot colored by selection
jpeg("dNdS_DEGs_selection_boxPlot.jpg")
ggplot(plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# box plot colored by selection with labels
#jpeg("dNdS_DEGs_selection_boxPlot_labeled.jpg")
#ggplot(plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
#  geom_boxplot() +
#  theme_minimal() +
#  ggrepel::geom_text_repel(data = labelSetAll, aes(label = labelSetAll$geneID), max.overlaps=100) +
#  scale_colour_discrete(type = plotColorSubset)
#dev.off()

# Length vs Effect

# view outliers
plotTable[plotTable$Length >= 5000,]

# subset Length values to remove outliers
plotSubset <- plotTable[plotTable$Length <= 5000,]

# box plot colored by selection
jpeg("length_DEGs_selection_boxPlot.jpg")
ggplot(plotSubset, aes(x=Effect, y=Length, fill=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# Length vs Selection

# box plot colored by selection
jpeg("length_selection_boxPlot.jpg")
ggplot(plotTable, aes(x=Selection, y=Length)) +
  geom_boxplot() +
  theme_minimal()
dev.off()


# scatter plots
# http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# https://stackoverflow.com/questions/15015356/how-to-do-selective-labeling-with-ggplot-geom-point
# http://www.sthda.com/english/wiki/ggplot2-point-shapes
# https://www.statology.org/ggplot-multiple-data-frames/

# view outliers
plotTable[plotTable$dS > 10,]

# subset dN dS values to remove outliers
plotSubset <- plotTable[plotTable$dS < 10,]

# setup facet
facetSelectionSubset <- factor(plotSubset$Selection, levels = c('Error', 'Positive', 'Negative'))

# dN vs dS

# scatter plot colored by effect
# shaped by selection
jpeg("dNdS_DEGs_selection_scatterPlot_jittered.jpg")
ggplot(plotSubset, aes(x=dS, y=dN, shape=Selection, color=Effect)) +
  theme_minimal() +
  geom_point() +
  geom_jitter() +
  scale_shape_manual(values=c(2, 1, 3)) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# scatter plot colored by effect
# faceted by selection
jpeg("dNdS_DEGs_selection_scatterPlot_faceted.jpg")
ggplot(plotSubset, aes(x=dS, y=dN, color=Effect)) +
  theme_minimal() +
  geom_point() +
  scale_shape_manual(values=c(2, 1, 3)) +
  facet_grid(rows = facetSelectionSubset, scales = "free_y") +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# dN vs Length

# scatter plot colored by effect
# shaped by selection
jpeg("dN_DEGs_selection_scatterPlot_jittered.jpg")
ggplot(plotSubset, aes(x=Length, y=dN, shape=Selection, color=Effect)) +
  theme_minimal() +
  geom_point() +
  geom_jitter() +
  scale_shape_manual(values=c(2, 1, 3)) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# scatter plot colored by effect
# faceted by selection
jpeg("dN_DEGs_selection_scatterPlot_faceted.jpg")
ggplot(plotSubset, aes(x=Length, y=dN, color=Effect)) +
  theme_minimal() +
  geom_point() +
  facet_grid(rows = facetSelectionSubset, scales = "free_y") +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# dS vs Length

# scatter plot colored by effect
# shaped by selection
jpeg("dS_DEGs_selection_scatterPlot_jittered.jpg")
ggplot(plotSubset, aes(x=Length, y=dS, shape=Selection, color=Effect)) +
  theme_minimal() +
  geom_point() +
  geom_jitter() +
  scale_shape_manual(values=c(2, 1, 3)) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# scatter plot colored by effect
# faceted by selection
jpeg("dS_DEGs_selection_scatterPlot_faceted.jpg")
ggplot(plotSubset, aes(x=Length, y=dS, color=Effect)) +
  theme_minimal() +
  geom_point() +
  facet_grid(rows = facetSelectionSubset, scales = "free_y") +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# dNdS vs Length

# setup facet
facetSelection <- factor(plotTable$Selection, levels = c('Error', 'Positive', 'Negative'))

# scatter plot colored by effect
# faceted by selection
jpeg("dNdS_length_DEGs_selection_scatterPlot_faceted.jpg")
ggplot(plotTable, aes(x=Length, y=dNdS, color=Effect)) +
  theme_minimal() +
  geom_point() +
  facet_grid(rows = facetSelection, scales = "free_y") +
  scale_colour_discrete(type = plotColorSubset)
dev.off()
