#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/stat%j.stdout
#SBATCH --error=std/stat%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="stat"
#SBATCH -p gcluster
#SBATCH --array=1-19

# software dependencies
## bcftools 1.17

## software
PATH=/gpool/cfiscus/bin:$PATH

# SET VARIABLES
THREADS=2
LOG=/gpool/cfiscus/vitis_snps/results/logs/005_"$SLURM_ARRAY_TASK_ID".log
RESULTS=/gpool/cfiscus/vitis_snps/results/
CHR=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2d\n", $1}')

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
cd "$RESULTS"/vcf
FILE=VITVarB40-14_v2.0_hap1_chr"$CHR".vcf.gz

## pull site info
bcftools query -f '%CHROM %POS %REF %ALT %AF %QD %FS %SOR %MQ %MQRankSum %ReadPosRankSum\n' $FILE | gzip > "$RESULTS"/vcf/"$CHR".vinfo.gz

## missingness report
plink2 --vcf "$FILE" --missing --memory 32000 --allow-extra-chr --out "$CHR"
