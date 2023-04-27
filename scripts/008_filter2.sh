#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=64G
#SBATCH --output=std/filt%j.stdout
#SBATCH --error=std/filt%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="filter2"
#SBATCH -p gcluster

# software dependencies
## bcftools 1.17

## software
PATH=/gpool/cfiscus/bin:$PATH
PATH=/gpool/bin/samtools-1.10/bin:$PATH

# SET VARIABLES
THREADS=2
LOG=/gpool/cfiscus/vitis_snps/results/logs/008.log
RESULTS=/gpool/cfiscus/vitis_snps/results

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
# top cut depth and filter very rare alleles
bcftools view -e'INFO/DP > (MEAN(INFO/DP) + STDEV(INFO/DP)) | MAF < 0.01' "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1.vcf.gz \
	| bgzip > "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_filter3.vcf.gz

# plink missingness report
plink2 --vcf "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_filter3.vcf.gz --missing --memory 64000 --allow-extra-chr --out filter3
