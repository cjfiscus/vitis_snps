#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --output=std/dict%j.stdout
#SBATCH --error=std/dict%j.stderr
#SBATCH --mail-user=cfiscus@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --job-name="dict"
#SBATCH -p gcluster

# software dependencies
## GATK 4.2.6.1
## bwa 0.7.12-r1039
## samtools/1.10
## JAVA 8

## software
PATH=/gpool/bin/samtools-1.10/bin:$PATH
PATH=/gpool/bin/jre1.8.0_221/bin/:$PATH
PATH=/gpool/cfiscus/bin/gatk-4.2.6.1:$PATH

# SET VARIABLES
REFERENCE=/gpool/cfiscus/vitis_snps/data/VITVarB40-14_v2.0.pseudomolecules.hap1.fasta
OUT=/gpool/cfiscus/vitis_snps/data/VITVarB40-14_v2.0.pseudomolecules.hap1.dict

# create .dict
gatk --java-options "-Xmx8G" CreateSequenceDictionary \
            -R "$REFERENCE" \
            -O "$OUT"

# index with bwa-mem2
bwa index "$REFERENCE"

# index with samtools
samtools faidx "$REFERENCE"
