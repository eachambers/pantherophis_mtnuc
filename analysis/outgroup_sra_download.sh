##download P. obsoletus fastq files from SRA and map to P. guttatus reference genome

#install SRA toolkit 
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -vxzf sratoolkit.tar.gz

#export PATH=$PATH:/work/05104/thomm/ls6/sratoolkit.3.0.5-ubuntu64/bin
#source ~/.bashrc

#download P. obsoletus fastq files from SRA (this was done with multiple fastqs (SRR10405221-SRR10405240)):
prefetch SRR10405221

fastq-dump --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR10405221/SRR10405221.sra

#put commands for downloading the whole set of fastqs in a job file to slurm
for i in SRR104052*; do
	echo "sratoolkit.3.0.5-ubuntu64/bin/fastq-dump --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $i/$i.sra">>fastqdump;
done

#rename fastq files 
for i in *_pass_1.fastq.gz; do
	mv $i ${i/_pass_1.fastq.gz}_R1.fastq.gz;
done

for i in *_pass_2.fastq.gz; do
	mv $i ${i/_pass_2.fastq.gz}_R2.fastq.gz;
done

#merge all R1 and R2 fastqs 
cat *R1.fastq.gz > PanObs_R1.fastq.gz

cat *R2.fastq.gz > PanObs_R2.fastq.gz

#map fastqs to P. guttatus genome using bowtie2
bowtie2 --no-unal --local -p 64 -x PanGut_3.0_genomic.fna -1 PanObs_R1_val_1.fq.gz -2 PanObs_R2_val_2.fq.gz -S PanObs.sam

#convert sam to bam
samtools view -bS PanObs.sam > PanObs.bam

#sort, remove duplicates, and add readgroups
samtools sort PanObs.bam -o PanObs_sorted.bam

picard MarkDuplicates -Xmx10g -INPUT PanObs_sorted.bam -OUTPUT PanObs_sorted_dedup.bam -METRICS_FILE PanObs_dupMetrics.txt -REMOVE_DUPLICATES true

AddOrReplaceReadGroups -I PanObs_sorted_dedup.bam -O PanObs_sorted_dedup_RG.bam -RGID PanObs -RGPL ILLUMINA -RGLB lib1 -RGPU unit1 -RGSM PanObs

#index final bam file
samtools index PanObs_sorted_dedup_RG.bam