#!/usr/bin/env Rscript

library(stringr)

##Parse TR gff3 file
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/TR_annotation_timema/")
TR_db_GW<-read.delim("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")

write.table(TR_db_GW, file = "Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)