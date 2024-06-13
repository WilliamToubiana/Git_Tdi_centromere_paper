#!/usr/bin/env Rscript

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/spades_assembly/levenstein")

table1<-read.table("output_distances_FF_contigs.txt", header=F, sep="\t", quote="")
table2<-read.table("output_distances_FR_contigs.txt", header=F, sep="\t", quote="")
colnames(table1)<-c("V1_FF", "V2_FF", "V3_FF")
colnames(table2)<-c("V1_FR", "V2_FR", "V3_FR")

library(dplyr)
table2_sub<-cbind(table1, table2)
table2_sub<-table2_sub %>% 
  rowwise() %>%
  mutate(min = min(V3_FF, V3_FR))

table2_sub1<-subset(table2_sub[c(1,2,7)], min<=65.8)
table2_sub1 <- table2_sub1[table2_sub1$V1_FF != table2_sub1$V2_FF,]

###Extract non-connected motifs

unique_all<-unique(table2$V1)
unique_connected<-unique(table2_sub1$V1_FF)
unique_unconnected<-setdiff(unique_all,unique_connected)
motifs_unconnected<-data.frame(unique_unconnected)


library(igraph)
library(GGally)
detach("package:purrr", unload=TRUE)#to reload "simplify"
detach("package:dyplr", unload=TRUE)#to reload "simplify"

net1 <- graph.data.frame(table2_sub1, directed=FALSE)
net1<-simplify(net1)


ggnet2(net1, node.size = 1, node.color = "black", edge.size = 0.1, size = "degree", label=FALSE, label.color = "red", label.size=2)+
  guides(size = FALSE)



###Extract connected nodes
components <- igraph::clusters(net1, mode="weak")
components[["csize"]]


split<-split(names(V(net1)), components[["membership"]])
library(purrr)
#split<-sapply(split, function(x) paste0("seq_", x), USE.NAMES=FALSE)

for (i in 1:length(split)) {
  assign(paste0("split", i), as.data.frame(split[[i]]))
}



###Merge the different families (included unconnected motifs) with Unique_TR_sequences_contigs_database

#load database
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/spades_assembly/levenstein")
TR_uniq<-read.table("Unique_TR_sequences_contigs_database.txt", header=T, sep="\t", quote="")

#merge
unconnected<-merge(motifs_unconnected, TR_uniq, by.x="unique_unconnected", by.y="TR_name", all.x = T, all.y = F)
unconnected$family<-"no_family"
colnames(unconnected)[1]<-"split[[i]]"

family1<-merge(split1, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family1$family<-"family1"
family2<-merge(split2, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family2$family<-"family2"
family3<-merge(split3, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family3$family<-"family3"
family4<-merge(split4, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family4$family<-"family4"
family5<-merge(split5, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family5$family<-"family5"
family6<-merge(split6, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family6$family<-"family6"
family7<-merge(split7, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family7$family<-"family7"
family8<-merge(split8, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family8$family<-"family8"
family9<-merge(split9, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family9$family<-"family9"
family10<-merge(split10, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family10$family<-"family10"
family11<-merge(split11, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family11$family<-"family11"
family12<-merge(split12, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family12$family<-"family12"
family13<-merge(split13, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family13$family<-"family13"
family14<-merge(split14, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family14$family<-"family14"
family15<-merge(split15, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family15$family<-"family15"
family16<-merge(split16, TR_uniq, by.x="split[[i]]", by.y="TR_name", all.x = T, all.y = F)
family16$family<-"family16"


family_all<-do.call(rbind, list(unconnected,family1, family2, family3, family4, family5, family6, family7, family8, family9, family10, family11, family12, 
                               family13, family14, family15, family16))
family_all$motif_size<-nchar(family_all$motif)
colnames(family_all)[1]<-"sequence"

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/spades_assembly/levenstein")
write.table(family_all, file = "Supplementary_table3.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
