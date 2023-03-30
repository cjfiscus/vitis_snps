#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=std/dl%j.stdout
#SBATCH --error=std/dl%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="geno"
#SBATCH -p gcluster
#SBATCH --array=1-19

# software dependencies
## GATK 4.2.6.1
## JAVA 8

## software
PATH=/gpool/bin/samtools-1.10/bin:$PATH
PATH=/gpool/bin/jre1.8.0_221/bin/:$PATH
PATH=/gpool/cfiscus/bin/gatk-4.2.6.1:$PATH
PATH=/gpool/cfiscus/bin:$PATH

# SET VARIABLES
THREADS=4
REFERENCE=/gpool/cfiscus/vitis_snps/data/VITVarB40-14_v2.0.pseudomolecules.hap1.fasta
LOG=/gpool/cfiscus/vitis_snps/results/logs/003_"$SLURM_ARRAY_TASK_ID".log
TEMP_DIR=/gpool/cfiscus/temp
SAMPLE_MAP=/gpool/cfiscus/vitis_snps/data/sample_map
RESULTS=/gpool/cfiscus/vitis_snps/results/
CHR=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2d\n", $1}')

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
# get filenames from list
TEMP_DIR="$TEMP_DIR"/003_"$NAME"
mkdir -pv "$TEMP_DIR"
echo "$TEMP_DIR"
cd "$TEMP_DIR"

##########
gatk --java-options "-Xmx32g" GenotypeGVCFs \
   -R "$REFERENCE" \
   -V gendb://"$RESULTS"/db/db_"$CHR" \
   -O "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_"$CHR".vcf.gz \
   --tmp-dir "$TEMP_DIR"
