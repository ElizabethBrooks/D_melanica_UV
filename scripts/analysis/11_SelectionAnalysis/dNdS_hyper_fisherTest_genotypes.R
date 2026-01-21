#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(ghibli)

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Genotypes"
setwd(workingDir)

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

# retrieve DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv")
genotypesTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv")

# keep only sig
interactionSig <- subset(interactionTable, interactionTable$FDR < 0.05)
treatmentSig <- subset(treatmentTable, treatmentTable$FDR < 0.05)
genotypesSig <- subset(genotypesTable, genotypesTable$FDR < 0.05)

# add effect tags
interactionSig$Effect <- "Interaction"
treatmentSig$Effect <- "Treatment"
genotypesSig$Effect <- "Tolerance"

# add geneID column
geneID <- row.names(positiveSubset)
positiveSubset <- cbind(geneID,positiveSubset)
geneID <- row.names(interactionSig)
interactionSig <- cbind(geneID,interactionSig)
geneID <- row.names(treatmentSig)
treatmentSig <- cbind(geneID,treatmentSig)
geneID <- row.names(genotypesSig)
genotypesSig <- cbind(geneID,genotypesSig)

# full outer join data frames
interactionSubset <- merge(x = interactionSig, y = positiveSubset, 
                          by = "geneID", all=TRUE)
treatmentSubset <- merge(x = treatmentSig, y = positiveSubset, 
                        by = "geneID", all=TRUE)
genotypesSubset <- merge(x = genotypesSig, y = positiveSubset, 
                          by = "geneID", all=TRUE)

# remove rows with NAs
positiveSubset <- na.omit(positiveSubset)
interactionSubset <- na.omit(interactionSubset)
treatmentSubset <- na.omit(treatmentSubset)
genotypesSubset <- na.omit(genotypesSubset)


# hypergeometric distribution (fishers test)
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/

# results data frame
pValues <- data.frame(
  effect = c("Interaction", "Treatment", "Tolerance"),
  enrichment = rep(99,3),
  depletion = rep(99,3)
)

# interaction
# initialize variables
positive <- nrow(positiveSubset)
de <- nrow(interactionSig)
overlap <- nrow(interactionSubset)
total <- nrow(interactionTable)

# test for over-representation (enrichment)
pValues$enrichment[1] <- phyper(overlap-1, de, total-de, positive,lower.tail= FALSE)

# test for under-representation (depletion)
pValues$depletion[1] <-  phyper(overlap, de, total-de, positive,lower.tail= TRUE)

# treatment
# initialize variables
positive <- nrow(positiveSubset)
de <- nrow(treatmentSig)
overlap <- nrow(treatmentSubset)
total <- nrow(treatmentTable)

# test for over-representation (enrichment)
pValues$enrichment[2] <- phyper(overlap-1, de, total-de, positive,lower.tail= FALSE)

# test for under-representation (depletion)
pValues$depletion[2] <-  phyper(overlap, de, total-de, positive,lower.tail= TRUE)

# genotypes
# initialize variables
positive <- nrow(positiveSubset)
de <- nrow(genotypesSig)
overlap <- nrow(genotypesSubset)
total <- nrow(genotypesTable)

# test for over-representation (enrichment)
pValues$enrichment[3] <- phyper(overlap-1, de, total-de, positive,lower.tail= FALSE)

# test for under-representation (depletion)
pValues$depletion[3] <-  phyper(overlap, de, total-de, positive,lower.tail= TRUE)

# write results to a csv file
write.csv(pValues, "fisherTest_positiveSelection_DEGs.csv", row.names=FALSE)
