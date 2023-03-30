#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --output=std/com%j.stdout
#SBATCH --error=std/com%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=4-00:00:00
#SBATCH --job-name="combine"
#SBATCH -p gcluster
#SBATCH --array=70-73

# software dependencies
## bcftools

## software
PATH=/gpool/bin/samtools-1.10/bin:$PATH
PATH=/gpool/cfiscus/bin:$PATH
PATH=/gpool/cfiscus/bin/gatk-4.2.6.1:$PATH
PATH=/gpool/bin/bcftools-1.10.2/bin:$PATH

## set vars
ID=$(head -n "$SLURM_ARRAY_TASK_ID" combine.txt | tail -n 1 | cut -f1)
FILES=$(head -n "$SLURM_ARRAY_TASK_ID" combine.txt | tail -n 1| cut -f3)
expression=($FILES)
OUT=$(echo "/gpool/cfiscus/vitis_snps/results/gvcf/""$ID"".g.vcf.gz")

##########
bcftools concat -a --no-version --rm-dups all "${expression[@]}" -O z -o "$OUT"
tabix -p vcf "$OUT"
