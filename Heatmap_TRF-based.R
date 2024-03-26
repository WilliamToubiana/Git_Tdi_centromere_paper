#!/usr/bin/env Rscript

###Heatmap

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/levenstein")
family_all<-read.table("Network_database.txt", header=TRUE, sep = "\t")

#Remove overlaps between repeats of the same family
family_all$family_chr<-paste(family_all$family, family_all$chr, sep="-")
colnames(family_all)[1]<-"sequence"

library(GenomicRanges)
myranges<-GRanges(seqnames=family_all$family_chr,ranges=IRanges(start=family_all$start,end=family_all$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

#Sum array lengths within each family per chromosome
nonoverlapWind_agg<-aggregate(nonoverlapWind$width ~ nonoverlapWind$seqnames, FUN=sum)
nonoverlapWind_agg[c('group', 'chr')] <- str_split_fixed(nonoverlapWind_agg$`nonoverlapWind$seqnames`, '-', 2)
nonoverlapWind_agg[c('Tdi', 'LR', 'scf')] <- str_split_fixed(nonoverlapWind_agg$chr, '_', 3)
nonoverlapWind_agg<-nonoverlapWind_agg[c(3,7,2)]
colnames(nonoverlapWind_agg)[3]<-"array_length"

#Build matrix 
library(tidyr)
matrix<-spread(nonoverlapWind_agg, scf, array_length)
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
matrix[is.na(matrix)] <- 0

#Color palette
library(grDevices)
paletteLength <- 50
myColor <- colorRampPalette(rev(c( rgb(255,42,20,maxColorValue = 255), 
                                   "#F7C5CC" ,"white", "#CCCCCC" ,"#999999" ,"#666666")))(paletteLength)
                                   
#Heatmap figure
library(pheatmap)
pheatmap(as.matrix(matrix), scale = "column", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", fontsize = 9, color= myColor)
