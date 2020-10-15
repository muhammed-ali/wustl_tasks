Task 7

Given the samples on task 6 have been successfully processed, define the necessary steps to generate a final joint vcf for these three samples.

Input files:
SAM alignment files that are named as sample*.sam (where * is 1, 2, 3)

Tools required:
GATK (GenomeAnalysisTK.jar)
Picard (picard.jar)
samtools
qualimap


Output files:
convert SAM to BAM: sample*.bam
sorted BAM files: sample*_sorted.bam
indexed BAM files:: sample*_sorted.bam.bai
flagstat_results: sample*_stats.txt
qualimap_results: 

Execution:
sh ./PLINK_script.bash


# http://zzz.bwh.harvard.edu/plink/download.shtml#download # installation
# https://www.cog-genomics.org/plink/1.9/assoc#linear
# http://zzz.bwh.harvard.edu/plink/anal.shtml#glm
# https://www.cog-genomics.org/plink/1.9/input#covar # how to define covariates
# https://www.biostars.org/p/292370/


# ./plink --bfile CHR1_Task_Renamed --linear --covar Phenotype.txt --covar-name age, gender, PC1, PC2 --parameters 1-4, 7

# ./plink --bfile CHR1_Task_Renamed --pheno Phenotype.txt --covar Phenotype.txt --covar-name age, gender, PC1, PC2 --out gwas_results --linear --adjust --ci 0.95
