#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("")

#Load the libraries
#library(filesstrings)
library(topGO)
library(edgeR)
#library(GO.db)
#library(reshape2)
library(ggplot2)
library(Rgraphviz)
#library(statmod)
library(tidyr)

# turn off scientific notation
options(scipen = 999)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
#workingDir <- args[1]
#workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance/GOAnalysis_ks_OLYM_30"
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA_9Nov2026/Genotypes/GOAnalysis_ks_OLYM_30"

# set working directory
dir.create(workingDir)
setwd(workingDir)

# retrieve subset tag
#set <- args[2]
set <- "OLYM"

# set the minimum module size
#minModSize <- args[3]
minModSize <- "30"

# retrieve WGCNA directory
#inDir <- args[4]
#inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance"

# retrieve gene to GO map
#GOmaps <- readMappings(file = args[5])
GOmaps <- readMappings(file = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ZQ/geneToGO_tagged_map.txt")

# retrieve dN dS values
#dNdSTable <- read.csv(file = args[6], row.names="geneID")
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# retrieve test type
#testType <- args[7]
testType <- "ks"

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
resultsTable$number <- resultsTable$number %>% replace_na('0')

# remove NAs
resultsTable <- na.omit(resultsTable)

# update the number of modules
numMods <- numMods + 1

# update the table of module colors
colorTable[nrow(colorTable) + 1,] <- c("None","0")


# create named list of all genes (gene universe) and p-values. The gene universe is set to be
# the list of all genes contained in the gene2GO list of annotated genes.
list_genes <- as.numeric(resultsTable$number)
list_genes <- setNames(list_genes, resultsTable$gene)
list_genes_filtered <- list_genes[names(list_genes) %in% names(GOmaps)]
  
# create data frame to hold the FDR, module color, total genes, significant genes, 
# and top significant BP, MF, and CC GO terms
moduleBPResults = data.frame(matrix(ncol = 10, nrow = numMods))
colnames(moduleBPResults) <- c("color","number","total","sigTotal","topID","topTerm","annotated","topSig","topExpected","topWeight")
moduleMFResults = data.frame(matrix(ncol = 10, nrow = numMods))
colnames(moduleMFResults) <- c("color","number","total","sigTotal","topID","topTerm","annotated","topSig","topExpected","topWeight")
moduleCCResults = data.frame(matrix(ncol = 10, nrow = numMods))
colnames(moduleCCResults) <- c("color","number","total","sigTotal","topID","topTerm","annotated","topSig","topExpected","topWeight")

#Loop through each module
for(j in 1:numMods){
  #Create function to return list of interesting DE genes (0 == not significant, 1 == significant)
  get_interesting_DE_genes <- function(geneUniverse){
    interesting_DE_genes <- rep(0, length(geneUniverse))
    for(i in 1:length(geneUniverse)){
      if(geneUniverse[i] == colorTable[j,2]){
        interesting_DE_genes[i] = 1
      }
    }
    interesting_DE_genes <- setNames(interesting_DE_genes, names(geneUniverse))
    return(interesting_DE_genes)
  }

  #Create topGOdata objects for enrichment analysis (1 for each ontology)
  BP_GO_data <- new('topGOdata', ontology = 'BP', allGenes = list_genes_filtered, 
                    geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                    gene2GO = GOmaps)
  MF_GO_data <- new('topGOdata', ontology = 'MF', allGenes = list_genes_filtered, 
                    geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                    gene2GO = GOmaps)
  CC_GO_data <- new('topGOdata', ontology = 'CC', allGenes = list_genes_filtered, 
                    geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                    gene2GO = GOmaps)
  
  #PerformGO enrichment on topGOdata object
  BP_GO_results <- runTest(BP_GO_data, statistic = 'Fisher')
  MF_GO_results <- runTest(MF_GO_data, statistic = 'Fisher')
  CC_GO_results <- runTest(CC_GO_data, statistic = 'Fisher')
  
  #GenTable to get statistics on GO terms
  list_BP_GO_terms <- usedGO(BP_GO_data)
  list_MF_GO_terms <- usedGO(MF_GO_data)
  list_CC_GO_terms <- usedGO(CC_GO_data)
  
  #Create summary data table
  BP_GO_results_table <- GenTable(BP_GO_data, weightFisher = BP_GO_results, orderBy = 'weightFisher', 
                                  topNodes = length(list_BP_GO_terms))
  MF_GO_results_table <- GenTable(MF_GO_data, weightFisher = MF_GO_results, orderBy = 'weightFisher', 
                                  topNodes = length(list_MF_GO_terms))
  CC_GO_results_table <- GenTable(CC_GO_data, weightFisher = CC_GO_results, orderBy = 'weightFisher', 
                                  topNodes = length(list_CC_GO_terms))

  # write table of GO terms to a CSV file
  write.table(BP_GO_results_table, file=paste(colorTable[j,1], "BP_GO_terms.csv", sep="_"), sep=",", row.names=FALSE, quote=FALSE)
  write.table(MF_GO_results_table, file=paste(colorTable[j,1], "MF_GO_terms.csv", sep="_"), sep=",", row.names=FALSE, quote=FALSE)
  write.table(CC_GO_results_table, file=paste(colorTable[j,1], "CC_GO_terms.csv", sep="_"), sep=",", row.names=FALSE, quote=FALSE)
  
  # create table of significant GO terms
  BP_sigGO_results_table <- BP_GO_results_table[BP_GO_results_table$weightFisher <= 0.05, ]
  MF_sigGO_results_table <- MF_GO_results_table[MF_GO_results_table$weightFisher <= 0.05, ]
  CC_sigGO_results_table <- CC_GO_results_table[CC_GO_results_table$weightFisher <= 0.05, ]

  # write table of significant GO terms to a CSV file
  write.table(BP_sigGO_results_table, file=paste(colorTable[j,1], "BP_sigGO_terms.csv", sep="_"), sep=",", row.names=FALSE, quote=FALSE)
  write.table(MF_sigGO_results_table, file=paste(colorTable[j,1], "MF_sigGO_terms.csv", sep="_"), sep=",", row.names=FALSE, quote=FALSE)
  write.table(CC_sigGO_results_table, file=paste(colorTable[j,1], "CC_sigGO_terms.csv", sep="_"), sep=",", row.names=FALSE, quote=FALSE)

  #Summary BP functions
  moduleBPResults$color[j] <- colorTable[j,1]
  moduleBPResults$number[j] <- colorTable[j,2]
  moduleBPResults$total[j] <- numGenes(BP_GO_data)
  moduleBPResults$sigTotal[j] <- length(sigGenes(BP_GO_data))
  moduleBPResults$topID[j] <- BP_GO_results_table[1,1]
  moduleBPResults$topTerm[j] <- BP_GO_results_table[1,2]
  moduleBPResults$annotated[j] <- BP_GO_results_table[1,3]
  moduleBPResults$topSig[j] <- BP_GO_results_table[1,4]
  moduleBPResults$topExpected[j] <- BP_GO_results_table[1,5]
  moduleBPResults$topWeight[j] <- BP_GO_results_table[1,6]
  
  #Summary MF functions
  moduleMFResults$color[j] <- colorTable[j,1]
  moduleMFResults$number[j] <- colorTable[j,2]
  moduleMFResults$total[j] <- numGenes(MF_GO_data)
  moduleMFResults$sigTotal[j] <- length(sigGenes(MF_GO_data))
  moduleMFResults$topID[j] <- MF_GO_results_table[1,1]
  moduleMFResults$topTerm[j] <- MF_GO_results_table[1,2]
  moduleMFResults$annotated[j] <- MF_GO_results_table[1,3]
  moduleMFResults$topSig[j] <- MF_GO_results_table[1,4]
  moduleMFResults$topExpected[j] <- MF_GO_results_table[1,5]
  moduleMFResults$topWeight[j] <- MF_GO_results_table[1,6]
  
  #Summary CC functions
  moduleCCResults$color[j] <- colorTable[j,1]
  moduleCCResults$number[j] <- colorTable[j,2]
  moduleCCResults$total[j] <- numGenes(CC_GO_data)
  moduleCCResults$sigTotal[j] <- length(sigGenes(CC_GO_data))
  moduleCCResults$topID[j] <- CC_GO_results_table[1,1]
  moduleCCResults$topTerm[j] <- CC_GO_results_table[1,2]
  moduleCCResults$annotated[j] <- CC_GO_results_table[1,3]
  moduleCCResults$topSig[j] <- CC_GO_results_table[1,4]
  moduleCCResults$topExpected[j] <- CC_GO_results_table[1,5]
  moduleCCResults$topWeight[j] <- CC_GO_results_table[1,6]
  
  #showGroupDensity example for most significant GO term
  BP_topSigGO_ID <- BP_GO_results_table[1, 'GO.ID']
  MF_topSigGO_ID <- MF_GO_results_table[1, 'GO.ID']
  CC_topSigGO_ID <- CC_GO_results_table[1, 'GO.ID']
  
  # create density plots and export to PDF
  #exportFile <- paste(colorTable[j,1], "TopSigGO_Density.pdf", sep="_")
  #pdf(file = exportFile)
  #showGroupDensity(BP_GO_data, BP_topSigGO_ID, ranks = TRUE)
  #showGroupDensity(MF_GO_data, MF_topSigGO_ID, ranks = TRUE)
  #showGroupDensity(CC_GO_data, CC_topSigGO_ID, ranks = TRUE)
  #dev.off()
  
  #printGraph to plot subgraphs induced by the most significant GO terms...saves to a file
  printGraph(BP_GO_data, BP_GO_results, firstSigNodes = 5, 
             fn.prefix = paste(colorTable[j,1],'BP_sigGO_subgraphs', sep = '_'), useInfo = 'all', pdfSW = TRUE)
  printGraph(MF_GO_data, MF_GO_results, firstSigNodes = 5, 
             fn.prefix = paste(colorTable[j,1], 'MF_sigGO_subgraphs', sep = '_'), useInfo = 'all', pdfSW = TRUE)
  printGraph(CC_GO_data, CC_GO_results, firstSigNodes = 5, 
             fn.prefix = paste(colorTable[j,1], 'CC_sigGO_subgraphs', sep = '_'), useInfo = 'all', pdfSW = TRUE)
  
  #showGroupDensity for UV tolerance assocaited GO terms
  #CC_DNA_repair_complex <- "GO:1990391"
  #BP_DNA_integrity_checkpoint <- "GO:0031570"
  #BP_Response_UV <- "GO:0009411"
  #BP_Mitotic_cell_cycle_checkpoint <- "GO:0007093"
  #BP_Cellular_response_DNA_damage_stimulus <- "GO:0006974"
  #MF_Single_stranded_DNA_binding <- "GO:0003697"
  #MF_Damaged_DNA_binding <- "GO:0003684"
  
  # create density plots for selected GO terms
  #exportFile <- paste(colorTable[j,1], "UVTGO_Density.pdf", sep="_")
  #pdf(file=exportFile)
  #showGroupDensity(CC_GO_data, CC_DNA_repair_complex, ranks = TRUE)
  #showGroupDensity(BP_GO_data, BP_DNA_integrity_checkpoint, ranks = TRUE)
  #showGroupDensity(BP_GO_data, BP_Response_UV, ranks = TRUE)
  #showGroupDensity(BP_GO_data, BP_Mitotic_cell_cycle_checkpoint, ranks = TRUE)
  #showGroupDensity(BP_GO_data, BP_Cellular_response_DNA_damage_stimulus, ranks = TRUE)
  #showGroupDensity(MF_GO_data, MF_Single_stranded_DNA_binding, ranks = TRUE)
  #showGroupDensity(MF_GO_data, MF_Damaged_DNA_binding, ranks = TRUE)
  #dev.off()
}


#Write the resulting tables to files
write.table(moduleBPResults, file="moduleTopGO_BPResults.csv", sep=",", row.names=FALSE)
write.table(moduleMFResults, file="moduleTopGO_MFResults.csv", sep=",", row.names=FALSE)
write.table(moduleCCResults, file="moduleTopGO_CCResults.csv", sep=",", row.names=FALSE)
