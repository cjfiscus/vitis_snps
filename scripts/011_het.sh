#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --output=std/vt%j.stdout
#SBATCH --error=std/vt%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="het"
#SBATCH -p gcluster


## software
PATH=/gpool/cfiscus/bin:$PATH

# SET VARIABLES
# xpclr
I="/gpool/cfiscus/vitis_snps/results/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz"
O="../results/stat"

vcftools --gzvcf "$I" --het --out "$O"/stathet
