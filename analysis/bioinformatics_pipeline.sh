### Bioinformatics pipeline starting with fastq files of raw sequence reads to dataset of called variants for high coverage samples
### All for loops in this file were used to produce a job script that was submitted to the Texas Advanced Computing Center (TACC) with slurm
### Programs were conda installed and run in environments created for each program

## Trim raw reads with trim-galore! Run this job in the same directory with fastq files
for i in *_R2_001.fastq.gz; do 
echo "trim_galore -j 8 -q 20 --phred33 --paired ${i/_R2_001.fastq.gz/}_R1_001.fastq.gz $i" >> trimpe; 
done


## Map reads to reference genome with bowtie2. Run this job in same directory as trimmed files and indexed reference genome. Process mapped alignments (sam files) into s

#Create index files for reference genome
bowtie2-build PanGut_3.0_genomic.fna PanGut_3.0_genomic.fna

#map trimmed files to reference
for i in *_R1_001_val_1.fq.gz; do
echo "bowtie2 --no-unal --local -p 8 -x PanGut_3.0_genomic.fna -1 $i -2 ${i/_R1_001_val_1.fq.gz/}_R2_001_val_2.fq.gz -S ${i/_R1_001_val_1.fq.gz/}.sam" >> bowtie2; 
done

## Process mapped alignments (sam files) into sorted, duplicate-removed, binary (bam) files with read groups added

#Convert sams to bams
for i in *.sam; do 
echo "samtools view -bS $i > ${i/.sam/}.bam" >> s2b; done
	
#Sort bams
for i in *.bam; do
echo "samtools sort $i -o ${i/.bam/}_sorted.bam" >> sort; done

#Remove duplicates with picard. Used -Xmx flag to increase max memory (10gb) 
for i in *_sorted.bam; do
echo "picard MarkDuplicates -Xmx10g -INPUT $i -OUTPUT ${i/_sorted.bam/}_sorted_dedup.bam -METRICS_FILE ${i/_sorted.bam/}_dupMetrics.txt -REMOVE_DUPLICATES true" >> dedup;
done

#Add read groups (add sample name to reads)
for i in *_sorted_dedup.bam; do 
echo "picard AddOrReplaceReadGroups -INPUT $i -OUTPUT ${i/.bam/}_RG.bam -RGID ${i/_sorted_dedup.bam/} -RGPL ILLUMINA -RGLB lib1 -RGPU unit1 -RGSM ${i/_sorted_dedup.bam/}" >> addRG; 
done

#Index bam files; can also be indexed using "samtools index $file"
for i in *RG.bam; do 
echo "picard BuildBamIndex -I $i" >> index;
done


## Call variants with FreeBayes (for High Coverage & Outgroup Samples only). Run in new directory with these samples only and indexed reference genome

# merge bam files of the four high coverage samples and outgroup
samtools merge -o HCSamples_PanObs.bam *RG.bam

# create freebayes job file to run using slurm
echo "freebayes-parallel <(fasta_generate_regions.py PanGut_3.0_genomic.fna.fai 100000) 56 \ -f PanGut_3.0_genomic.fna --gvcf -g 500 HCSamples_PanObs.bam > HCSamples_PanObs_fb.vcf" >> freebayes


## filter vcf

#remove indels
vcftools --vcf HCSamples_PanObs_fb.vcf --remove-indels --recode --recode-INFO-all --out HCSamples_PanObs_fb_snps

#bgzip new vcf to subsample and look at vcfstats (see vcfstats.txt for more info)
bgzip HCSamples_PanObs_fb_snps.recode.vcf
	
#filter for quality, depth of coverage, and missing data
vcftools --gzvcf HCSamples_PanObs_fb_snps.recode.vcf.gz --minQ 30 --max-missing 0.80 --min-meanDP 10 --minDP 10 --recode --recode-INFO-all --out HCSamples_PanObs_fb_snps_filtered

