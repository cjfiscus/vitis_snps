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
LOG=/gpool/cfiscus/vitis_snps/results/logs/009.log
RESULTS=/gpool/cfiscus/vitis_snps/results

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
# filter samps with > 70% missing and sites with > 5% missing
bcftools view -S ^id_rm2.txt -e 'F_MISSING > 0.05' \
	"$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_full_filter3.vcf.gz | bgzip > "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz

# index
tabix "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz

# produce sample list
bcftools query -l "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz > "$RESULTS"/filtered_samples.txt

# plink missingness report
plink2 --vcf "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz --missing --memory 64000 --allow-extra-chr --out filtered
