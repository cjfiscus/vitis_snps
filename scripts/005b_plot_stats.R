#!/usr/bin/env Rscript
# vcf stats plotting script
# cjfiscus
# 2020-09-01
# 
# usage
# Rscript plot_stats.R raw_stats.txt.gz 0.10 out
# arg2[1] is stats file
# args[2] is proportion of variants to sample
# args[3] is output prefix
##########
# parse args and libs
args = commandArgs(trailingOnly=TRUE)

library(pacman)
p_load(data.table, ggplot2)

# read in data
df<-fread(args[1])
names(df)<-c("CHROM", "POS", "REF", "ALT", "AF", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")

# samp rand X% dataset
size<-round(nrow(df)*as.numeric(args[2]), 0)
sub<-df[sample(nrow(df), size),]
rm(df)

# plot stats
## QD
sub$QD<-as.numeric(as.character(sub$QD))
p1<-ggplot(sub, aes(x=QD)) + geom_density() + theme_classic()
out<-paste0(args[3], "_QD.jpeg")
ggsave(out, p1, height=4, width=5)

## FS
sub$FS<-as.numeric(as.character(sub$FS))
p1<-ggplot(sub, aes(x=log(FS))) + geom_density() + theme_classic()
out<-paste0(args[3], "_FS.jpeg")
ggsave(out, p1, height=4, width=5)

## SOR
sub$SOR<-as.numeric(as.character(sub$SOR))
p1<-ggplot(sub, aes(x=SOR)) + geom_density() + theme_classic()
out<-paste0(args[3], "_SOR.jpeg")
ggsave(out, p1, height=4, width=5)

## MQ
sub$MQ<-as.numeric(as.character(sub$MQ))
p1<-ggplot(sub, aes(x=MQ)) + geom_density() + theme_classic()
out<-paste0(args[3], "_MQ.jpeg")
ggsave(out, p1, height=4, width=5)


## MQRankSum
sub$MQRankSum<-as.numeric(as.character(sub$MQRankSum))
p1<-ggplot(sub, aes(x=MQRankSum)) + geom_density() + theme_classic()
out<-paste0(args[3], "_MQRankSum.jpeg")
ggsave(out, p1, height=4, width=5)

## ReadPosRankSum
sub$ReadPosRankSum<-as.numeric(as.character(sub$ReadPosRankSum))
p1<-ggplot(sub, aes(x=ReadPosRankSum)) + geom_density() + theme_classic()
out<-paste0(args[3], "_ReadPosRankSum.jpeg")
ggsave(out, p1, height=4, width=5)
