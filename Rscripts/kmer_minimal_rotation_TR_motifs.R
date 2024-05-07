#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(spgs)


setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/mer_centromere/spades_assembly")
TR<-read.table("contigs.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt",
               header=T, col.names = c("chr", "start", "end", "copy_nb", "motif_length", "motif", "motif_minimal_rotation"))


#select enriched motifs with distinct sequences
TR_uniq<-data.frame(unique(TR$motif_minimal_rotation))
colnames(TR_uniq)<-"motif"

#duplicate motifs until reaching 329bp in length
x = TR_uniq$motif
n = 329 #largest motif as reference
TR_uniq$motif_extended<-sapply(strsplit(x, ""), function(s){
  paste(rep(s, length.out = n), collapse = "")
})


#generate reverse complements
TR_uniq$motif_extended_RC<-reverseComplement(TR_uniq$motif_extended, case=c("upper"))


#prepare files for Levenshtein analysis
TR_uniq$TR_name<-paste("seq", 1:nrow(TR_uniq), sep = "_")

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/spades_assembly/levenstein")
write.table(TR_uniq[c(4,3)], file = "Unique_TR_sequences_contigs_reverseComp.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_uniq[c(4,2)], file = "Unique_TR_sequences_contigs.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_uniq, file = "Unique_TR_sequences_contigs_database.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


