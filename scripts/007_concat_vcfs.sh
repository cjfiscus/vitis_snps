#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --output=std/filt%j.stdout
#SBATCH --error=std/filt%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="concat"
#SBATCH -p gcluster

# software dependencies
## bcftools 1.17

## software
PATH=/gpool/cfiscus/bin:$PATH
PATH=/gpool/bin/samtools-1.10/bin:$PATH

# SET VARIABLES
THREADS=2
LOG=/gpool/cfiscus/vitis_snps/results/logs/007.log
RESULTS=/gpool/cfiscus/vitis_snps/results

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
bcftools concat --threads 7 -f concat_lst.txt -o "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1.vcf.gz
