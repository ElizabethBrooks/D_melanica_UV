#!/usr/bin/env Rscript

# load libraries
library(tidyr)

# turn off scientific notation
#options(scipen = 999)

# turn on scientific notation
#options(scipen = 0)

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


# retrieve DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
genotypesTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")

# keep only sig
interactionSig <- interactionTable[interactionTable$FDR < 0.05,]
treatmentSig <- treatmentTable[treatmentTable$FDR < 0.05,]
genotypesSig <- genotypesTable[genotypesTable$FDR < 0.05,]

# add effect tags
interactionSig$Effect <- "Interaction"
treatmentSig$Effect <- "Treatment"
genotypesSig$Effect <- "Tolerance"

# add geneID column
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
geneID <- row.names(interactionSig)
interactionSig <- cbind(geneID,interactionSig)
geneID <- row.names(treatmentSig)
treatmentSig <- cbind(geneID,treatmentSig)
geneID <- row.names(genotypesSig)
genotypesSig <- cbind(geneID,genotypesSig)

# keep necessary columns
dNdSSubset <- dNdSSubset[,c("geneID","dN","dS","dNdS")]
treatmentSubset <- treatmentSig[,c("geneID","logFC","FDR","Effect")]
genotypesSubset <- genotypesSig[,c("geneID","logFC","FDR","Effect")]
interactionSubset <- interactionSig[,c("geneID","logFC","FDR","Effect")]


# import annotation data (emailed from ZQ)
#annotationTable <- read.delim(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/KAP4_NCBI_functional_annotation.txt", row.names="Gene_locus")

# add geneID column
#geneID <- row.names(annotationTable)
#annotationTable <- cbind(geneID,annotationTable)

# add gene- tag to each ID
#annotationTable$geneID <- paste("gene", annotationTable$geneID, sep="-")

# https://ncbiinsights.ncbi.nlm.nih.gov/2022/02/24/gene-ontology-refseq/
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/#:~:text=Schema%3A%20Genome%20GFF3-,Genome%20GBFF,features%20annotated%20on%20each%20genome.
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gbff/
# https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
# import annotation data (NCBI gbff)
annotationTable <- read.delim(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/KAP4_NCBI_functional_annotation.txt", row.names="Gene_locus")


# import potentially locally adapted UV responsive gene pathway associates
pathwayTable <- read.delim(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/dMelUV_UVResponseGenes_pathways_16Aug2023.tsv", row.names="Gene")

# add geneID column
geneID <- row.names(pathwayTable)
pathwayTable <- cbind(geneID,pathwayTable)

# add gene- tag to each ID
pathwayTable$geneID <- paste("gene", pathwayTable$geneID, sep="-")


# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# combine all tables
effectTable <- rbind(interactionSubset, treatmentSubset, genotypesSubset)

# full outer join
fullTable <- merge(x = resultsTable, y = effectTable, 
                   by = "geneID", all=TRUE)

# full outer join
fullTable <- merge(x = fullTable, y = annotationTable, 
                   by = "geneID", all=TRUE)

# full outer join
fullTable <- merge(x = fullTable, y = pathwayTable, 
                   by = "geneID", all=TRUE)


# force scientific notation
fullTable$logFC <- formatC(fullTable$logFC, format = "e", digits = 2) 
fullTable$FDR <- formatC(fullTable$FDR, format = "e", digits = 2) 

# convert columns to character type
fullTable$start <- as.character(fullTable$start)
fullTable$end <- as.character(fullTable$end)

# set tags for NAs
fullTable$Effect <- fullTable$Effect %>% replace_na('None')
fullTable$AssociatedPathways <- fullTable$AssociatedPathways %>% replace_na('None')
fullTable$FDR[fullTable$FDR == ' NA'] <- 'None'
fullTable$logFC[fullTable$logFC == ' NA'] <- 'None'
fullTable$chr. <- fullTable$chr. %>% replace_na('None')
fullTable$strand <- fullTable$strand %>% replace_na('None')
fullTable$description1 <- fullTable$description1 %>% replace_na('None')
fullTable$description2 <- fullTable$description2 %>% replace_na('None')
fullTable$GO.terms <- fullTable$GO.terms %>% replace_na('None')
fullTable$X <- fullTable$X %>% replace_na('None')
fullTable$start <- fullTable$start %>% replace_na('None')
fullTable$end <- fullTable$end %>% replace_na('None')

# remove NAs
cleanTable <- na.omit(fullTable)

# fill empty rows
cleanTable$description2 <- sub("^$", "None", cleanTable$description2)
cleanTable$GO.terms <- sub("^$", "None", cleanTable$GO.terms)
cleanTable$X <- sub("^$", "None", cleanTable$X)


# subset genes associated with interesting pathways or under positive selection
subsetTable <- cleanTable[cleanTable$AssociatedPathways != "None" | cleanTable$Selection == "Positive",]

# write table to tsv file
write.table(subsetTable, file = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/dMelUV_UVResponseGenes_networkModules_16Aug2023.tsv", sep = "\t")


# select specific columns
subsetColumns <- subsetTable[,c("geneID","AssociatedPathways","color","number","logFC","FDR","Effect","dN","dS","dNdS","Selection","description1","description2","GO.terms","X")]

# write table to tsv file
write.table(subsetColumns, file = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/dMelUV_UVResponseGenes_networkModules_16Aug2023_subset.tsv", sep = "\t")

