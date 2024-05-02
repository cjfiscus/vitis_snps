#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="trios"
#SBATCH -p gcluster

# software dependencies
## Dsuite

## software
PATH=/gpool/cfiscus/bin:$PATH
PATH=/gpool/cfiscus/bin/Dsuite/Build:$PATH
VCF=/gpool/cfiscus/vitis_snps/results/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz
SET=/gpool/cfiscus/vitis_snps/data/vitis_set.txt
RESULTS=/gpool/cfiscus/vitis_snps/results/dsuite
##########
cd "$RESULTS"

Dsuite Dtrios -g -o alldata "$VCF" "$SET" 
