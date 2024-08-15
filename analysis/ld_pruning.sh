### LD pruning CZ NMT and control gene sets

# Run within data directory of GitHub

./bcftools +prune -m 0.6 -w 50 cznmtsnps.vcf -O v -o cznmtsnps_ldp.vcf
./bcftools +prune -m 0.6 -w 50 czcontsnps.vcf -O v -o czcontsnps_ldp.vcf

# Verify number of sites remaining
./bcftools query -f '%POS\n' cznmtsnps.vcf | wc -l # 38,551 sites in NMT vcf originally
./bcftools query -f '%POS\n' cznmtsnps_ldp.vcf | wc -l # 14,330 sites remain in LD-pruned NMT vcf

./bcftools query -f '%POS\n' czcontsnps.vcf | wc -l # 163,958 sites in control vcf originally
./bcftools query -f '%POS\n' czcontsnps_ldp.vcf | wc -l # 61,697 sites remain in LD-pruned control vcf

# Convert vcf to fasta file

python vcf2phylip.py -i cznmtsnps_ldp.vcf --fasta -m 1
python vcf2phylip.py -i czcontsnps_ldp.vcf --fasta -m 1