#!/bin/bash -l

task6_dir=/path/to/task6/home/directory
task7_dir=/path/to/task7/home/directory/
refGenome=/path/to/reference/genome/fasta/
GATK_path=/path/to/GenomeAnalysisTK.jar
picard_path=/path/to/picard.jar

# First we need to move the alignment files (SAM) of three samples to Task7 directory
mv $task6_dir/Sample*/*.sam $task7_dir

# assumption: here we assume that alignment files are named as sample*.sam (where * is 1, 2, 3) and exist in ${task7_dir}

cd ${task7_dir}

# convert SAM to BAM
# In order to save space and make it easy for softwares to hand the large alignments files, 
# we convert SAM files to a binary format called BAM.
for i in *.sam; do name=$(echo $i | cut -d "." -f 1); bam='.bam'; samtools view -bS $i > $name$bam; done

# sort the BAM files
for i in *.bam; do name=$(echo $i | cut -d "." -f 1); bam='_sorted.bam'; samtools sort $i -o $name$bam; done

# create the BAM files indexes
for i in *_sorted.bam; do name=$(echo $i | cut -d "." -f 1); samtools index $i; done

# calculate flagstat for every aligned file in the directory. 
# This will provide you the alignment statistics for each sample.
# for example, how many reads are there in the dataset, and how many of them properly aligned to the reference genome.
mkdir -p flagstat_results
for i in *_sorted.bam; do name=$(echo $i | cut -d "_" -f 1); stat='_stats.txt'; samtools flagstat $i > ./flagstat_results/$name$stat; done

# check the mapping quality
# This tool provides the basic statistics of the alignment (number of reads, coverage, GC-content, etc.), 
# also a number of useful graphs are produced.
mkdir -p qualimap_results
for i in *_sorted.bam; do name=$(echo $i | cut -d "_" -f 1); result='_result.pdf'; qualimap bamqc -nt 8 -bam $i -outfile ./qualimap_results/$name$result; done


# mark duplicates
# Potential PCR duplicates are marked with Picard Tools
for i in *_sorted.bam
 do \ 
 	name=$(echo $i | cut -d "." -f 1) \ 
 	markDupl='_markDupl.bam'; metric='_markDupl_metrics.txt' \ 
 	java -Xmx8g -jar $picard_path MarkDuplicates \ 
  	I=$i \ 
  	O=$name$markDupl \ 
  	METRICS_FILE=$name$metric \ 
 done


# add or replace read groups
# The GATK pipeline requires read group information in BAM files. This information is used to differentiate samples 
# and to detect artifacts associated with sequencing techniques.
for i in *_sorted_markDupl.bam
do \ 
	name=$(echo $i | cut -d "." -f 1) \ 
	readgroup='_readgroup.bam' \ 
	java -Xmx4g -jar $picard_path AddOrReplaceReadGroups \ 
	I=$i \ 
	O=$name$readgroup \ 
	RGLB=library \ 
	RGPL=illumina \ 
	RGPU=barcode \ 
	RGSM=sample \ 
done


# index the resulting BAM files
for i in *_sorted_markDupl_readgroup.bam; do name=$(echo $i | cut -d "." -f 1); samtools index $i; done


# realign Indels
# The local realignment process is designed to locally realign reads such that the 
# number of mismatching bases is minimized across all the reads.
# The first step is to determine suspicious intervals which are likely in need of realignment by using RealignerTargetCreator.
for i in *_sorted_markDupl_readgroup.bam
do \ 
	name=$(echo $i | cut -d "." -f 1) \ 
	intervals='.intervals'; log='.intervals.log' \ 
	java -Xmx4g -jar $GATK_path \ 
	-T RealignerTargetCreator \ 
	-R $refGenome \ 
	-I $i \ 
	-O $name$intervals \ 
	-log $name$log \ 
done


# perform realignment
# The second step is to run the realigner over the intervals determined in previous step.
for i in *_sorted_markDupl_readgroup.bam
do \ 
	name=$(echo $i | cut -d "." -f 1) \ 
	realigned='_realigned.bam' \ 
	java -Xmx4g -jar $GATK_path \ 
	-T IndelRealigner \ 
	-R $refGenome \ 
	-I $i \ 
	-targetIntervals $name$intervals \ 
	-O $name$realigned \ 
done


# generate GVCF files
# The BAM files are ready now to be processed with the GATK to obtain the genotype data.
for i in *_sorted_markDupl_readgroup_realigned.bam
do \ 
	name=$(echo $i | cut -d "." -f 1) \ 
	gvcf='.g.vcf' \ 
	java -Xmx4g -jar $GATK_path \ 
	-T HaplotypeCaller \ 
	-R $refGenome \ 
	-I $i \ 
	-O $name$gvcf \ 
	-ERC GVCF \ 
done


# genotype all GVCF files
# Use the GenotypeGVCFs in the GATK to genotype all GVCF files. 
# It is necessary because estimation of some population genetics statistics require scaling by a total number of sites.
java -Xmx8g -jar $GATK_path -T GenotypeGVCFs -R $refGenome \ 
-V sample1_sorted_markDupl_readgroup_realigned.g.vcf \ 
-V sample2_sorted_markDupl_readgroup_realigned.g.vcf \ 
-V sample3_sorted_markDupl_readgroup_realigned.g.vcf \ 
-O GVCFall.vcf


# check quality of GVCF file
# Given a variant callset, the VariantEval calculates various quality control metrics.
java -Xmx8g -jar $GATK_path VariantEval -R $refGenome -O GVCFall.eval.grp --eval GVCFall.vcf

echo "GATK pipeline executed successfully"

