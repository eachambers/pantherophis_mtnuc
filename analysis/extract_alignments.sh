### Extracting phylogenomic blocks, or alignments (scripts modified from Else K. Mikkelsen's github page, https://github.com/elsemikk/Stercorariidae_Phylogenomics)

# create dataset with variant AND invariant sites from which to extract alignments. Run in directory with merged bam file (including high coverage samples plus outgroup) and indexed reference genome. Run job script using slurm. 
echo "freebayes-parallel <(fasta_generate_regions.py PanGut_3.0_genomic.fna.fai 100000) 256 \ -f PanGut_3.0_genomic.fna --gvcf --report-monomorphic -g 500 HCSamples_PanObs.bam > HCSamples_allsites.vcf" >> freebayes

# filter SNPs and invariant sites in allsites file; this creates a snps file and an invariant sites file which we will then concatenate into one file
# genotype quality is not filtered here because of how it would differentially affect variant and invariant sites
bcftools view -m2 -M3 -v snps -e 'AVG(FMT/DP)>80 || MQM<20 || MQMR<20' HCSamples_allsites.vcf | bcftools filter -S . -e 'FMT/DP<10' -Oz -o HCSamples.filteredSNPs.vcf.gz

bcftools view -m1 -M1 -e 'AVG(FMT/DP)>80 || MQMR<20' HCSamples_allsites.vcf | bcftools filter -S . -e 'FMT/DP<10' -Oz -o HCSamples.invariant.vcf.gz

bcftools concat --threads 24 --allow-overlaps HCSamples.filteredSNPs.vcf.gz HCSamples.invariant.vcf.gz -Oz -o HCSamples.filtered.allsites.vcf.gz

mkdir contig_fastas

## EXTRACT 5-kb ALIGNMENT BLOCKS ACROSS GENOME FOR SPECIES TREE ANALYSIS

# make a list of contigs from reference genome that are at least 5-kb long and save this list in text file called "contigs"
# REGIONS=[horizontal list of contigs (separated by a space) that are at least 5kb long; broke these regions into 3 sections - REGIONS, REGIONS2, REGIONS3 - due to length]

#generate consensus sequences for each sample and contig 
parallel --colsep " " 'samtools faidx PanGut_3.0_genomic.fna {2}: | bcftools consensus --sample {1} --haplotype I --absent N --missing N --include '\''TYPE="snp" || TYPE="ref"'\'' HCSamples.filtered.allsites.vcf.gz > contig_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $REGIONS

parallel --colsep " " 'samtools faidx PanGut_3.0_genomic.fna {2}: | bcftools consensus --sample {1} --haplotype I --absent N --missing N --include '\''TYPE="snp" || TYPE="ref"'\'' HCSamples.filtered.allsites.vcf.gz > contig_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $REGIONS2

parallel --colsep " " 'samtools faidx PanGut_3.0_genomic.fna {2}: | bcftools consensus --sample {1} --haplotype I --absent N --missing N --include '\''TYPE="snp" || TYPE="ref"'\'' HCSamples.filtered.allsites.vcf.gz > contig_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $REGIONS3

#replace headers with sample names
parallel 'sed -i "s/>.*$/>{1}/g" contig_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $REGIONS

parallel 'sed -i "s/>.*$/>{1}/g" contig_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $REGIONS2

parallel 'sed -i "s/>.*$/>{1}/g" contig_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $REGIONS3

#combine samples together
cat contigs | while read scaffold ; do cat contig_fastas/"$scaffold".*.fa |  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | tr "\t" "\n" > contig_fastas/"$scaffold".fasta ; done

cd contig_fastas

#download extract_blocks.rb from https://github.com/mmatschiner/tutorials/find/master
#generate 5000kb blocks with no more than 20% missing data and save as nexus files
cat ../contigs | parallel ruby extract_blocks.rb {1}.fasta blocks 5000 0.2

#count how many blocks were made
ls blocks | wc -l

#We need a stricter filter to produce a more manageable number of alignments
#generate 5000kb blocks with no more than 3% missing data and save as nexus files
cat contigs | parallel ruby /work/05104/thomm/ls6/extract_blocks.rb {1}.fasta blocks 5000 0.03

#This produced 216 5kb alignments. We made an ordered list of the coordinates of each alignment and subtracted the ending coordinate from the starting coordinate of each adjacent alignment in the list. We then removed alignments so that the final list included only alignments that were at least 5kb apart (i.e., not adjacent).
#This resulted in a total of 208 5-kb blocks


## EXTRACT N-MT GENE REGIONS

#return to directory with HCSamples_allsites.vcf and indexed reference genome

mkdir nmt_haplotypes

#NMT = [horizontal list of N-mt gene coordinates, separated by space]

#generate consensus sequences for each sample and gene
parallel --colsep " " 'samtools faidx PanGut_3.0_genomic.fna {2} | bcftools consensus --sample {1} --haplotype I --absent N --missing N --include '\''TYPE="snp" || TYPE="ref"'\'' HCSamples.filtered.allsites.vcf.gz > nmt_haplotypes/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $NMT

#replace headers with sample names
parallel 'sed -i "s/>.*$/>{1}/g" nmt_haplotypes/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $NMT

#combine samples together
cat nmt | while read scaffold ; do cat nmt_haplotypes/"$scaffold".*.fa |  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | tr "\t" "\n" > nmt_haplotypes/"$scaffold".fasta ; done


## EXTRACT CONTROL GENE REGIONS

#return to directory with HCSamples_allsites.vcf and indexed reference genome

mkdir control_fastas

#CONT = [horizontal list of control gene coordinates, separated by space]

#generate consensus sequences for each sample and gene
parallel --colsep " " 'samtools faidx PanGut_3.0_genomic.fna {2} | bcftools consensus --sample {1} --haplotype I --absent N --missing N --include '\''TYPE="snp" || TYPE="ref"'\'' HCSamples.filtered.allsites.vcf.gz > control_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $CONT

#replace headers with sample names
parallel 'sed -i "s/>.*$/>{1}/g" control_fastas/{2}.{1}.diploid.fa' ::: H16057 H21189 O42736 PanObs TJH3395 ::: $CONT

#combine samples together
cat controls | while read scaffold ; do cat control_fastas/"$scaffold".*.fa |  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | tr "\t" "\n" > control_fastas/"$scaffold".fasta ; done
