#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(spgs)

###Coverage ChIP and input from TR annotation within enriched 10kb regions
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/coverage/categories_enriched")

TR_cenh3<-read.table("Tdi_cenh3_testes_1_TR_10kb_coverage.txt",
                     header=FALSE)
TR_input<-read.table("Tdi_input_testes_1_TR_10kb_coverage.txt",
                     header=FALSE)

TR<-cbind(TR_cenh3,TR_input[4])
colnames(TR)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

TR$coverage_cenh3_norm<-TR$coverage_cenh3/69361547 #normalized by the number of mapped reads
TR$coverage_input_norm<-TR$coverage_input/75167800 #normalized by the number of mapped reads
TR$ratiocenH3_norm<-TR$coverage_cenh3_norm/TR$coverage_input_norm
TR$log2cenH3_norm<-log2(TR$ratiocenH3_norm)
TR$genomic_feature<-"TR"

TR$array<-paste(TR$Scaffold_name, TR$start, sep="_")
TR$array<-paste(TR$array, TR$stop, sep="_")


###TR annotation within enriched 10kb regions (minimal rotation on motif sequences)
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_10kb_regions/")
TR2<-read.table("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation_enriched_10kb_regions.txt",
                header=FALSE, col.names = c("chr", "start", "end", "motif"))

TR2$array<-paste(TR2$chr, TR2$start, sep="_")
TR2$array<-paste(TR2$array, TR2$end, sep="_")

merge<-merge(TR2, TR[c(11,9)], by="array", all.x = T, all.y = F)
merge <- merge %>% distinct()
merge<-merge %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


#filter out motifs with log2(ChIP/Input) < 2
TR_enriched<-subset(merge, log2cenH3_norm>=2)

#select enriched motifs with distinct sequences
TR_enriched_uniq<-data.frame(unique(TR_enriched$motif))
TR_enriched_uniq$sequence<-paste("seq", 1:nrow(TR_enriched_uniq), sep = "_")
colnames(TR_enriched_uniq)[1]<-"motif"


#duplicate enriched motifs until reaching 2000bp in length
x = TR_enriched_uniq$motif
n = 2000 #largest motif as reference
TR_enriched_uniq$motif_2kb<-sapply(strsplit(x, ""), function(s){
  paste(rep(s, length.out = n), collapse = "")
})

#generate reverse complements
TR_enriched_uniq$motif_2kb_RC<-reverseComplement(TR_enriched_uniq$motif_2kb, case=c("upper"))

#convert into fasta
TR_enriched_uniq$sequence <- paste0(">",TR_enriched_uniq$sequence)
fasta <- c(rbind(TR_enriched_uniq$sequence, TR_enriched_uniq$motif))

##Enriched motif database
TR_enriched_db<-merge(TR_enriched, TR_enriched_uniq, by="motif", all.x = T, all.y = F)


setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/levenstein")
write.table(TR_enriched_uniq[2:3], file = "Enriched_2kbmotifs_enriched_10kb_regions.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_enriched_uniq[c(2,4)], file = "Enriched_2kbmotifs_reverseComp_enriched_10kb_regions.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_enriched_db, file = "Enriched_motifs_database.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/blast")
write(x = fasta, file = "Enriched_motifs_enriched_10kb_regions.fasta")



