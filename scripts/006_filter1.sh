#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --output=std/filt%j.stdout
#SBATCH --error=std/filt%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="filter"
#SBATCH -p gcluster
#SBATCH --array=1-19

# software dependencies
## bcftools 1.17

## software
PATH=/gpool/cfiscus/bin:$PATH
PATH=/gpool/bin/samtools-1.10/bin:$PATH

# SET VARIABLES
THREADS=2
LOG=/gpool/cfiscus/vitis_snps/results/logs/006_"$SLURM_ARRAY_TASK_ID".log
RESULTS=/gpool/cfiscus/vitis_snps/results
CHR=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2d\n", $1}')

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
cd "$RESULTS"/vcf
FILE=VITVarB40-14_v2.0_hap1_chr"$CHR".vcf.gz

## biallelic snps, GATK hard filters, QUAL
bcftools view -e'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0 | QUAL < 20' -m2 -M2 -v snps "$FILE" \
        | bgzip > "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_chr"$CHR"_filter1.vcf.gz

## set calls with ind. depth to missing
bcftools filter -e 'FMT/DP<1' -S . "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_chr"$CHR"_filter1.vcf.gz \
	| bgzip > "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_chr"$CHR"_filter2.vcf.gz
