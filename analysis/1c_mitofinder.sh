###Extract mitogenomes from each fastq file using Mitofinder
### All analyses were run on the Lonestar 6 (LS6) high performance computing system at the Texas Advanced Computing Center (TACC), UT Austin; job scripts were submitted to the LS6 queue using slurm

#download MitoFinder using singularity 
singularity pull --arch amd64 library://remiallio/default/mitofinder:v1.4.1 

#Add MitoFinder to path
cdw
p=$(pwd)
echo -e "\n#Path to MitoFinder image \nexport PATH=\$PATH:$p" >> ~/.bashrc 
source ~/.bashrc  

#run Mitofinder on low coverage, trimmed sample fastqs using a slowinskii reference mitogenome
for i in *_R1_001_val_1.fq.gz; do 
echo "mitofinder_v1.4.1.sif -j ${i/_R1_001_val_1.fq.gz} -1 $i -2 ${i/_R1_001_val_1.fq.gz}_R2_001_val_2.fq.gz -r P_slowinskii_mt_genome.gb -o 2 -p 32" >> mito; done

##using Mitofinder on high coverage samples is computationally intensive. Mitogenomes can be pulled from a subset of the reads. Reads can be subsampled using seqtk. This command subsamples 20M reads from each file.

#In separate directory with high coverage and outgroup samples only, use this code to subsample fastqs
for i in *fq.gz; do
seqtk sample -s100 $i 20000000 > Mitofinder/${i/.fq.gz/}_sub.fq.gz;
done

#then extract mitogenomes with the same code used for lower coverage samples
for i in *_R1_001_val_1_sub.fq.gz; do 
echo "mitofinder_v1.4.1.sif -j ${i/_R1_001_val_1_sub.fq.gz} -1 $i -2 ${i/_R1_001_val_1_sub.fq.gz}_R2_001_val_2_sub.fq.gz -r P_slowinskii_mt_genome.gb -o 2 -p 32" >> mito; done
