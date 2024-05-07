#!/usr/bin/env Rscript

###Heatmap

#Blast table
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/blast")

blast<-read.table("Enriched_motifs_enriched_10kb_regions_blastn.txt",
                  header=FALSE, col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), sep = "\t")

blast$perc_cov<-blast$length/blast$qlen

subset_blast<-subset(blast, pident>=80 & perc_cov>=0.8) 


#Network database
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/levenstein")
family_all<-read.table("Network_database.txt", header=TRUE, sep = "\t")


#Associate repeat families with blast hits (including motifs with no sequence similarity with other motifs) 

merge<-merge(subset_blast, family_all[c(1,2,10)], by.x="qseqid", by.y="sequence", all.x = T, all.y = F) #adding the motif and repeat family it belongs to
merge <- merge %>% distinct()
merge$array<-abs(merge$send-merge$sstart)


#Convert all blast hits to the same strand
merge2 <- merge %>% mutate(i_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))

merge2 <- merge2 %>% mutate(i_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))


#Remove overlaps between repeats of the same family
library(GenomicRanges)
merge2$family_wind<-paste(merge2$family, merge2$sseqid, sep="-")
myranges<-GRanges(seqnames=merge2$family_wind,ranges=IRanges(start=merge2$i_start,end=merge2$i_end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind[c('family_scf', 'position')] <- str_split_fixed(nonoverlapWind$seqnames, ':', 2)


#Sum hit lengths within each family per chromosome
nonoverlapWind_agg<-aggregate(nonoverlapWind$width ~ nonoverlapWind$family_scf, FUN=sum)
nonoverlapWind_agg[c('family', 'chr')] <- str_split_fixed(nonoverlapWind_agg$`nonoverlapWind$family_scf`, '-', 2)
nonoverlapWind_agg[c('Tdi', 'LR', 'scf')] <- str_split_fixed(nonoverlapWind_agg$chr, '_', 3)
nonoverlapWind_agg<-nonoverlapWind_agg[c(3,7,2)]
colnames(nonoverlapWind_agg)[3]<-"hit_length"


#Build matrix
library(tidyr)
matrix<-spread(nonoverlapWind_agg, scf, hit_length)
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


##Generate bed file
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/blast")
merge2$seq_family<-paste(merge2$qseqid, merge2$family, sep="_")
write.table(merge2[c(2,18,19,21)], file = "Repeat_families_10kb_regions_blastn.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

##Generate Supplementary table 4
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/blast")
merge2[c('scaffold', 'window')] <- str_split_fixed(merge2$sseqid, ':', 2)
write.table(merge2[c(1,15,16,18,19,22,23)], file = "Supplementary_table_4.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
