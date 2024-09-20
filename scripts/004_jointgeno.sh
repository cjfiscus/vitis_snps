#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100G
#SBATCH --output=std/gt%j.stdout
#SBATCH --error=std/gt%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="geno"
#SBATCH -p gcluster
#SBATCH --array=19

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
SAMPLE_MAP=/gpool/cfiscus/vitis_snps/data/sample_map2
RESULTS=/gpool/cfiscus/vitis_snps/results/
CHR=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2d\n", $1}')
INT=/gpool/cfiscus/vitis_snps/results/vcf/VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
# get filenames from list
TEMP_DIR="$TEMP_DIR"/003_"$CHR"
mkdir -pv "$TEMP_DIR"
echo "$TEMP_DIR"
cd "$TEMP_DIR"

##########
export TILEDB_DISABLE_FILE_LOCKING=1
gatk --java-options "-Xmx100g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
   -R "$REFERENCE" \
   -V gendb://"$RESULTS"/db/db_"$CHR" \
   -O "$RESULTS"/vcf/VITVarB40-14_v2.0_hap1_chr"$CHR"_allsites.vcf.gz \
   -all-sites \
   -L chr${CHR}.filtered.interval_list \
   -ip 100 \
   --genomicsdb-shared-posixfs-optimizations true \
   --tmp-dir "$TEMP_DIR"
