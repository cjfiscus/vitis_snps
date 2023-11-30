#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=std/db%j.stdout
#SBATCH --error=std/db%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --job-name="db_imp"
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
TRIMMOMATIC=/gpool/cfiscus/bin/Trimmomatic-0.39/trimmomatic-0.39.jar

# SET VARIABLES
THREADS=4
LOG=/gpool/cfiscus/vitis_snps/results/logs/002_"$SLURM_ARRAY_TASK_ID".log
TEMP_DIR=/gpool/cfiscus/temp
SAMPLE_MAP=/gpool/cfiscus/vitis_snps/data/sample_map2
RESULTS=/gpool/cfiscus/vitis_snps/results/db
CHR=$(echo "$SLURM_ARRAY_TASK_ID" | awk '{printf "%.2d\n", $1}')

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
# get filenames from list
echo "$CHR"

TEMP_DIR="$TEMP_DIR"/002_"$CHR"
mkdir -pv "$TEMP_DIR"
echo "$TEMP_DIR"
cd "$TEMP_DIR"

##########
gatk --java-options "-Xmx64g -Xms4g" GenomicsDBImport \
	--sample-name-map "$SAMPLE_MAP" \
	--genomicsdb-workspace-path "$RESULTS"/db_"$CHR" \
	--tmp-dir "$TEMP_DIR" \
	--batch-size 50 \
	--reader-threads 3 \
	-L VITVarB40-14_v2.0.hap1.chr${CHR}
