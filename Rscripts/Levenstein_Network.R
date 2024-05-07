#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/levenstein")

table1<-read.table("levenstein_distances_FF.txt", header=F, sep="\t", quote="")
table2<-read.table("levenstein_distances_FR.txt", header=F, sep="\t", quote="")
colnames(table1)<-c("V1_FF", "V2_FF", "V3_FF")
colnames(table2)<-c("V1_FR", "V2_FR", "V3_FR")

library(dplyr)
table2_sub<-cbind(table1, table2)
table2_sub<-table2_sub %>% 
  rowwise() %>%
  mutate(min = min(V3_FF, V3_FR))

table2_sub1<-subset(table2_sub[c(1,2,7)], min<=400)
table2_sub1 <- table2_sub1[table2_sub1$V1_FF != table2_sub1$V2_FF,]


###Network

library(igraph)
library(GGally)
detach("package:purrr", unload=TRUE)#to reload "simplify"
detach("package:dyplr", unload=TRUE)#to reload "simplify"

net1 <- graph.data.frame(table2_sub1, directed=FALSE)
net1<-simplify(net1)


ggnet2(net1, node.size = 1, node.color = "black", edge.size = 0.1, size = "degree", label=FALSE, label.color = "red", label.size=2)+
  guides(size = "none")


###Extract non-connected motifs

unique_all<-unique(table2$V1)
unique_connected<-unique(table2_sub1$V1_FF)
unique_unconnected<-setdiff(unique_all,unique_connected)
motifs_unconnected<-data.frame(unique_unconnected)


#Extract connected motifs (i.e., repeat families)

components <- igraph::clusters(net1, mode="weak")
components[["csize"]]

split<-split(names(V(net1)), components[["membership"]])

library(purrr)
for (i in 1:length(split)) {
  assign(paste0("split", i), as.data.frame(split[[i]]))
}


###Associate repeat families with enriched motif database (including motifs with no sequence similarity with other motifs) 

#Enriched motif database
TR_enriched_db<-read.table("Enriched_motifs_database.txt", header=TRUE, sep = "\t")

#Merge datasets
unconnected<-merge(motifs_unconnected, TR_enriched_db, by.x="unique_unconnected", by.y="sequence", all.x = T, all.y = F)
unconnected$family<-"no_family"
colnames(unconnected)[1]<-"split[[i]]"

family1<-merge(split1, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family1$family<-"family1"
family2<-merge(split2, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family2$family<-"family2"
family3<-merge(split3, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family3$family<-"family3"
family4<-merge(split4, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family4$family<-"family4"
family5<-merge(split5, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family5$family<-"family5"
family6<-merge(split6, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family6$family<-"family6"
family7<-merge(split7, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family7$family<-"family7"
family8<-merge(split8, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family8$family<-"family8"
family9<-merge(split9, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family9$family<-"family9"
family10<-merge(split10, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family10$family<-"family10"
family11<-merge(split11, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family11$family<-"family11"
family12<-merge(split12, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family12$family<-"family12"
family13<-merge(split13, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family13$family<-"family13"
family14<-merge(split14, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family14$family<-"family14"
family15<-merge(split15, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family15$family<-"family15"
family16<-merge(split16, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family16$family<-"family16"
family17<-merge(split17, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family17$family<-"family17"
family18<-merge(split18, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family18$family<-"family18"
family19<-merge(split19, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family19$family<-"family19"
family20<-merge(split20, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family20$family<-"family20"
family21<-merge(split21, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family21$family<-"family21"
family22<-merge(split22, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family22$family<-"family22"
family23<-merge(split23, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family23$family<-"family23"
family24<-merge(split24, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family24$family<-"family24"
family25<-merge(split25, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family25$family<-"family25"
family26<-merge(split26, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family26$family<-"family26"
family27<-merge(split27, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family27$family<-"family27"
family28<-merge(split28, TR_enriched_db, by.x="split[[i]]", by.y="sequence", all.x = T, all.y = F)
family28$family<-"family28"


family_all<-do.call(rbind, list(unconnected,family1, family2, family3, family4, family5, family6, family7, family8, family9, family10, family11, family12, 
                               family13, family14, family15, family16, family17, family18, family19, family20, family21, family22,
                               family23, family24, family25, family26, family27, family28))


write.table(family_all, file = "Network_database.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


##Supplementary table1

family_all$motif_size<-nchar(family_all$motif)
supp_table1<-family_all[,c(1, 2, 8, 9, 10, 12, 4, 5, 6)]
colnames(supp_table1)<-c("sequence", "motif", "motif_extended", "motif_extended_RC", "group", "motif_size", "scaffold", "scaffold_start", "scaffold_end")

write.table(supp_table1, file = "Supplementary_table2.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
