#!/bin/bash

for i in /gpool/cfiscus/vitis_snps/results/gvcf/*.vcf.gz
do
## parse number
ID=$(basename "$i" | sed 's/.g.vcf.gz//g')

## parse sample name
SAMP=$(zcat $i | head -n 200 | grep "#CHROM" | cut -f10)

echo "$ID""\t""$SAMP""\t""/gpool/cfiscus/vitis_snps/results/gvcf/""$ID".vcf.gz
done
