#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantFiltering_splitSites_jobOutput

# script to perform variant filtering after calling
# usage: qsub variantFiltering_splitSites_bcftools.sh

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_splitSites_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsDir=$inputsDir"/variantsMerged"

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# set output folder
outFolder=$inputsDir

# name output file of inputs
inputOutFile=$outFolder"/variantFiltering_summary.txt"
# name output file of filtering info
outputsFile=$outFolder"/variantFiltering_stats.txt"

# add version to output file of inputs
bcftools --version > $inputOutFile

# check total variants
echo "Total variants from filtered reads with MQ > 60: " > $outputsFile
bcftools view --threads 4 $inputsDir"/"$type"_calls.norm.bcf" | grep -v "#" | wc -l >> $outputsFile

# include sites with quality > 20 
bcftools filter --threads 4 -i '%QUAL>20' $inputsDir"/"$type"_calls.norm.bcf" -Ob -o $outFolder"/"$type"_calls.flt-qual.bcf"
echo "bcftools filter --threads 4 -i '%QUAL>20' "$inputsDir"/"$type"_calls.norm.bcf -Ob -o "$outFolder"/"$type"_calls.flt-qual.bcf" >> $inputOutFile
echo "& including sites with quality > 20: " >> $outputsFile
bcftools view --threads 4 $inputsDir"/"$type"_calls.flt-qual.bcf" | grep -v "#" | wc -l >> $outputsFile

# include sites with average read depth > 10
bcftools filter --threads 4 -i 'INFO/DP>10' $outFolder"/"$type"_calls.flt-qual.bcf" -Ob -o $outFolder"/"$type"_calls.flt-DP.bcf"
echo "bcftools filter --threads 4 -i 'INFO/DP>10' "$outFolder"/"$type"_calls.flt-qual.bcf -Ob -o "$outFolder"/"$type"_calls.flt-DP.bcf" >> $inputOutFile
echo "& including sites with average read depth > 10: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-DP.bcf" | grep -v "#" | wc -l >> $outputsFile

# include sites that have a MQSB value
# MQSB is a Mann-Whitney U test of Mapping Quality vs Strand Bias
# see MQB in https://www.reneshbedre.com/blog/vcf-fields.html
bcftools filter --threads 4 -i 'INFO/MQSB>0.05' $outFolder"/"$type"_calls.flt-DP.bcf" -Ob -o $outFolder"/"$type"_calls.flt-MQSB.bcf"
echo "bcftools filter --threads 4 -i 'INFO/MQSB>0.05' "$outFolder"/"$type"_calls.flt-DP.bcf -Ob -o "$outFolder"/"$type"_calls.flt-MQSB.bcf" >> $inputOutFile
echo "& including variants with MQSB > 0.05: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-MQSB.bcf" | grep -v "#" | wc -l >> $outputsFile

# include site with 8 alleles, which indicates no uncalled genotypes
# integer total number of alleles in called genotypes
bcftools filter --threads 4 -i 'INFO/AN=8' $outFolder"/"$type"_calls.flt-MQSB.bcf" -Ob -o $outFolder"/"$type"_calls.flt-AN.bcf"
echo "bcftools filter --threads 4 -i 'INFO/AN=8' "$outFolder"/"$type"_calls.flt-MQSB.bcf -Ob -o "$outFolder"/"$type"_calls.flt-AN.bcf" >> $inputOutFile
echo "& including sites with 8 alleles: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-AN.bcf" | grep -v "#" | wc -l >> $outputsFile

# include variants with 8 alleles, which indicates no multiallelic sites
# integer allele count in genotypes for each ALT allele, in the same order as listed
bcftools filter --threads 4 -i 'INFO/AC=8' $outFolder"/"$type"_calls.flt-AN.bcf" -Ob -o $outFolder"/"$type"_calls.flt-AC.bcf"
echo "bcftools filter --threads 4 -i 'INFO/AC=8' "$outFolder"/"$type"_calls.flt-AN.bcf -Ob -o "$outFolder"/"$type"_calls.flt-AC.bcf" >> $inputOutFile
echo "& including sites with 8 ALT alleles: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-AC.bcf" | grep -v "#" | wc -l >> $outputsFile

# exclude heterozygous sites
bcftools filter --threads 4 -e 'GT="het"' $outFolder"/"$type"_calls.flt-AC.bcf" -Ob -o $outFolder"/"$type"_calls.flt-homo.bcf"
echo "bcftools filter --threads 4 -e 'GT=\"het\"' "$outFolder"/"$type"_calls.flt-AC.bcf -Ob -o "$outFolder"/"$type"_calls.flt-homo.bcf" >> $inputOutFile
echo "& excluding heterozygous sites: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-homo.bcf" | grep -v "#" | wc -l >> $outputsFile

# exclude sites with an IMF < 1
# Maximum fraction of reads supporting an indel
# https://github.com/samtools/bcftools/issues/911
bcftools filter --threads 4 -e 'INFO/IMF<1' $outFolder"/"$type"_calls.flt-homo.bcf" -Ob -o $outFolder"/"$type"_calls.flt-IDV.bcf"
echo "bcftools filter --threads 4 -e 'INFO/IMF<1' "$outFolder"/"$type"_calls.flt-homo.bcf -Ob -o "$outFolder"/"$type"_calls.flt-IDV.bcf" >> $inputOutFile
echo "& excluding sites with an IMF < 1: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-IDV.bcf" | grep -v "#" | wc -l >> $outputsFile

# turn on left alignment, normalize indels, and collapse multi allelic sites
bcftools norm --threads 4 -m -any -f $genomeFile $outFolder"/"$type"_calls.flt-IDV.bcf" -Ob -o $outFolder"/"$type"_calls.flt-norm.bcf"
echo "bcftools norm --threads 4 --write-index -m +any -f "$genomeFile" "$outFolder"/"$type"_calls.flt-IDV.bcf -Ob -o "$outFolder"/"$type"_calls.flt-norm.bcf" >> $inputOutFile
echo "& with left alignment, normalized indels, and collapsed multi allelic sites: " >> $outputsFile
bcftools view --threads 4 $outFolder"/"$type"_calls.flt-norm.bcf" | grep -v "#" | wc -l >> $outputsFile

# index the vcf
bcftools index $outputsPath"/"$runNum"_calls.flt-norm.bcf"

# clean up
#rm $outFolder"/"$type"_calls.flt-qual.bcf"
#rm $outFolder"/"$type"_calls.flt-DP.bcf"
#rm $outFolder"/"$type"_calls.flt-MQSB.bcf"
#rm $outFolder"/"$type"_calls.flt-AN.bcf"
#rm $outFolder"/"$type"_calls.flt-homo.bcf"b
#rm $outFolder"/"$type"_calls.flt-IDV.bcf"
