mamba activate mtnuc

# The following script performs linkage disequilibrium pruning of N-mt and control gene datasets in contact zone samples

##          vcfs for NMTs and control loci (cznmtsnps.vcf & czcontsnps.vcf)
##          fastas for NMTs and control loci (cznmtsnps.min4.fasta & czcontsnps.min4.fasta) <- vcfs converted to fasta using vcf2phylip.py
##          mitotype assignments (cz_mitotypes.txt)
##          vcfs for LD-pruned NMTs and control loci (cznmtsnps_ldp.vcf & czcontsnps_ldp.vcf) <- generated using `ld_pruning.sh` script
##          fastas for LD-pruned NMTs and control loci (cznmtsnps_ldp.min1.fasta & czcontsnps_ldp.min1.fasta) <- vcfs converted to fasta using vcf2phylip.py

# Index the vcf such that it can be used with bcftools
bcftools index ../../data/cz_snps.vcf.gz # or tabix cz_snps.vcf.gz

# Extract relevant control gene regions from vcf
bcftools view -R ../../data/control_coords2.txt -Oz -o czcontsnps2.vcf.gz ../../data/cz_snps.vcf.gz

# Control genes vcf should have 73 chroms
bcftools query -f '%CHROM\n' czcontsnps2.vcf.gz | sort | uniq | wc -l

# Convert vcf to fasta file
python vcf2phylip.py -i czcontsnps2.vcf.gz --fasta

### LD pruning CZ NMT and control gene sets

# Specify location of plugins directory for bcftools +prune function
export BCFTOOLS_PLUGINS=/Users/eac/Documents/GitHub/pantherophis_mtnuc/software/bcftools-1.21/plugins

# Run LD pruning
bcftools +prune -m 0.6 -w 50 cznmtsnps.vcf -O v -o cznmtsnps_ldp.vcf
bcftools +prune -m 0.6 -w 50 czcontsnps2.vcf.gz -O v -o czcontsnps2_ldp.vcf

# Verify number of sites remaining
bcftools query -f '%POS\n' cznmtsnps.vcf | wc -l # 38,551 sites in NMT vcf originally
bcftools query -f '%POS\n' cznmtsnps_ldp.vcf | wc -l # 14,330 sites remain in LD-pruned NMT vcf

bcftools query -f '%POS\n' czcontsnps2.vcf.gz | wc -l # 157,026 sites in control vcf originally
bcftools query -f '%POS\n' czcontsnps2_ldp.vcf | wc -l # 59,208 sites remain in LD-pruned control vcf

mv czcontsnps2* ../../data/

# Convert vcf to fasta file
mamba deactivate
python vcf2phylip.py -i ../data/czcontsnps2_ldp.vcf --fasta -m 1

mv czcontsnps2* ../data/

