#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --output=std/xpclr%j.stdout
#SBATCH --error=std/xpclr%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="sel"
#SBATCH -p gcluster
#SBATCH --array=1-19

# software dependencies
## xpclr

## software
PATH=/gpool/cfiscus/bin:$PATH
PATH=/gpool/cfiscus/bin/miniconda3/bin:$PATH
source activate xpclr

# SET VARIABLES
# xpclr
A=../data/cluster1.txt
B=../data/cluster2.txt
I="/gpool/cfiscus/vitis_snps/results/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz"
F="vcf"
O="../results/xpclr"
C_LST=../data/chr_lst.txt
C=$(head -n "$SLURM_ARRAY_TASK_ID" "$C_LST" | tail -n 1)

xpclr -O "$O"/"$C"_xpclr.txt -F "$F" -I "$I" -Sa "$A" -Sb "$B" -C "$C" --verbose 10 --maxsnps 100 --ld 0.90 --step 5000
