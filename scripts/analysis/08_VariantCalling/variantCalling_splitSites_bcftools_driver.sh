#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N VC_splitSites_pipeline_jobOutput

# script to perform bam read quaity filtering
# usage: qsub variantCalling_splitSites_bcftools_driver.sh
# usage Ex: qsub variantCalling_splitSites_bcftools_driver.sh

#Required modules for ND CRC servers
module load bio

# set sorted folder name
sortedFolder="sortedCoordinate_samtoolsHisat2_run1"

# retrieve aligned reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")

# set sorting outputs absolute path
inputsDir=$inputsPath"/"$sortedFolder

# set filter type
filterType="filteredMapQ"

# set and create merged outputs directory
outputPath=$inputsPath"/variantsCalled_splitSites_samtoolsBcftools"
mkdir $outputPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# run script to filter bam files by mapq
bash filterByMapQ_samtools.sh

# merge alignments for each genotype
bash sortingMerge_samtools.sh

# run variant calling script
bash variantCallingMerged_bcftools.sh

# run variant filtering script
#bash variantFiltering_bcftools.sh
bash variantFiltering_splitSites_bcftools.sh

# run script to generate consensus sequences
bash generateConsensusMerged_bcftools.sh
#bash generateConsensusGenotype_bcftools.sh "E05"
#bash generateConsensusGenotype_bcftools.sh "R2"
#bash generateConsensusGenotype_bcftools.sh "Y05"
#bash generateConsensusGenotype_bcftools.sh "Y023"

# incorporate variants and update the sample gff
bash incorporateVariants_liftOver.sh

# clean up
#rm -r $inputsDir
