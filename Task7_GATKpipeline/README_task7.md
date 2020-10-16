## Task 7  

Given the samples on task 6 have been successfully processed, define the necessary steps to generate a final joint vcf for these three samples. Provide the pipeline you would use to run these samples with as much detail as possible, including necessary reference files. No need to actually generate the final files for this task.  

### Input files:  
SAM alignment files from the task 6 output that are named as sample*.sam (where * is 1, 2, 3)  


### Tools required:  
GATK (GenomeAnalysisTK.jar)  
Picard (picard.jar)  
samtools  
qualimap  


### Output files:  
convert SAM to BAM: sample*.bam  
sorted BAM files: sample*_sorted.bam  
indexed BAM files: sample*_sorted.bam.bai  
flagstat results: sample*_stats.txt  
qualimap results: sample*_result.pdf  
mark duplicates: sample*_sorted_markDupl.bam  
add or replace read groups: sample*_sorted_markDupl_readgroup.bam  
realign indels: sample*_sorted_markDupl_readgroup_.intervals and sample*_sorted_markDupl_readgroup_.intervals.log  
indel realigner: sample*_sorted_markDupl_readgroup_realigned.bam  
generate gvcf files: sample*_sorted_markDupl_readgroup_realigned.g.vcf  
genotype all gvcf files: GVCFall.vcf  
check quality of gvcf files: GVCFall.eval.grp


### Execution:  
`sh ./GATKpipeline_script.sh`  

### Brief summary of the pipeline:  

#### Convert .SAM to .BAM files:  
In order to save space and make it easy for softwares to hand the large alignments files, we convert SAM files to a binary format called BAM.  

#### Samtools flagstat:  
The samtool `flagstat` will provide you the alignment statistics for each sample. For example, statistics about how many reads are there in the dataset, and how many of them properly aligned to the reference genome.  

#### Mapping quality check:  
The `qualimap` tool provides the basic statistics of the alignment (number of reads, coverage, GC-content, etc.), 
It also produces a number of useful graphs for checking the quality of mapping.  

#### Mark duplicates:  
Potential PCR duplicates need to be marked with Picard tools. Marking duplicates is advised even if you have used a PCR-free library preparation procedure because reads identified as duplicates are not removed and can be included in the subsequent analyses if needed.  

#### Add or replace read groups:  
The GATK pipeline requires read group information in BAM files. This information is used to differentiate samples and to detect artifacts associated with sequencing techniques.  

#### Realign Indels:  
The local realignment process is designed to locally realign reads such that the number of mismatching bases is minimized across all the reads. It consists of two steps:  
+ Determine suspicious intervals which are likely in need of realignment by using RealignerTargetCreator.  
+ Run the realigner over the intervals determined in previous step by using IndelRealigner.  

#### Generate GVCF files:  
Use `HaplotypeCaller` in the GATK to obtain the genotype data.  

#### Genotype all GVCF files:  
Use the `GenotypeGVCFs` in the GATK to genotype all GVCF files. It is necessary because estimation of some population genetics statistics require scaling by a total number of sites.  

#### check quality of GVCF file:  
Given a variant callset, the `VariantEval` calculates various quality control metrics.  
