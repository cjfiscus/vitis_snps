#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=33G
#SBATCH --output=std/dl%j.stdout
#SBATCH --error=std/dl%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=6-00:00:00
#SBATCH --job-name="map"
#SBATCH -p gcluster
#SBATCH --array=270

# 142-300 will be array
# software dependencies
## GATK 4.2.6.1
## bwa 0.7.12-r1039
## samtools/1.10
## JAVA 8

## software
PATH=/gpool/bin/samtools-1.10/bin:$PATH
PATH=/gpool/bin/jre1.8.0_221/bin/:$PATH
PATH=/gpool/cfiscus/bin/gatk-4.2.6.1:$PATH
PATH=/gpool/cfiscus/bin:$PATH
TRIMMOMATIC=/gpool/cfiscus/bin/Trimmomatic-0.39/trimmomatic-0.39.jar

# SET VARIABLES
THREADS=4
REFERENCE=/gpool/cfiscus/vitis_snps/data/VITVarB40-14_v2.0.pseudomolecules.hap1.fasta
ADAPTERSPE=/gpool/cfiscus/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa
ADAPTERSSE=/gpool/cfiscus/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa
LOG=/gpool/cfiscus/vitis_snps/results/logs/"$SLURM_ARRAY_TASK_ID".log
SEQLIST=/gpool/cfiscus/vitis_snps/data/reads_2_gvcf_lst.txt
TEMP_DIR=/gpool/cfiscus/temp

# SET LOGS
exec 1>>${LOG}
exec 2>>${LOG}
###########
# get filenames from list
NUM=$(echo $SLURM_ARRAY_TASK_ID - 140 | bc) 
FILE=$(head -n $NUM $SEQLIST | tail -n 1 | cut -f8)

# obtain sample id
NAME=$(head -n $NUM $SEQLIST | tail -n 1 | cut -f1)
echo "$NAME"
SAMPLE=$(head -n $NUM $SEQLIST | tail -n 1 | cut -f2)
echo "$SAMPLE"
FILE=$(head -n "$NUM" "$SEQLIST" | tail -n 1 | cut -f3)

# make temp directory and go there
TEMP_DIR="$TEMP_DIR"/"$NAME"
mkdir -pv "$TEMP_DIR"
echo "$TEMP_DIR"
cd "$TEMP_DIR"

# Quality/Adapter trimming
echo "trimming with trimmomatic..."
java -jar "$TRIMMOMATIC" PE -threads "$THREADS" \
    "$FILE"-READ1-Sequences.txt.gz \
    "$FILE"-READ2-Sequences.txt.gz \
    "$NAME"_1_trimmed_paired.fq.gz "$NAME"_1_unpaired.fq.gz \
    "$NAME"_2_trimmed_paired.fq.gz "$NAME"_2_unpaired.fq.gz \
    ILLUMINACLIP:"$ADAPTERSPE":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:60

# map to reference genome
echo "mapping with bwa..."
bwa mem -t "$THREADS" -M $REFERENCE "$NAME"_1_trimmed_paired.fq.gz \
	"$NAME"_2_trimmed_paired.fq.gz > "$NAME".sam

# sam to sorted bam
echo "samtools sam to bam"
samtools view -bS "$NAME".sam | samtools sort -T temp - -o "$NAME".bam
##########
# mark dups
gatk --java-options "-Xmx32G" MarkDuplicates \
    -I "$NAME".bam \
    -O "$NAME".nodups.bam \
    -M "$LOG".md.metrics.txt 

## calc stats
samtools index "$NAME".nodups.bam
mosdepth -t "$THREADS" /gpool/cfiscus/vitis_snps/results/stats/"$NAME" "$NAME".nodups.bam
samtools flagstat "$NAME".nodups.bam > /gpool/cfiscus/vitis_snps/results/stats/"$NAME".flagstat

# add read groups
gatk --java-options "-Xmx32G" AddOrReplaceReadGroups \
    -I "$NAME".nodups.bam \
    -O "$NAME".nodups.rg.bam \
    --CREATE_INDEX true \
    -RGID "$NAME" \
    -RGLB "$SAMPLE" \
    -RGSM "$SAMPLE" \
    -RGPL Illumina_HiSeq \
    -RGPU Illumina_HiSeq  

# haplotype caller
gatk --java-options "-Xmx32G" HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$NAME".nodups.rg.bam \
    -O /gpool/cfiscus/vitis_snps/results/gvcf/"$NAME".g.vcf.gz \
    --base-quality-score-threshold 20 \
    --sample-ploidy 2 \
    --native-pair-hmm-threads 2 \
    --tmp-dir "$TEMP_DIR" \
    -ERC GVCF

# clean up temp
#rm -r "$TEMP_DIR"
