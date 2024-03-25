#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/coverage/categories_enriched")

###Coverage sequence categories within enriched 10kb regions

#TRs
TR_cenh3<-read.table("Tdi_cenh3_testes_1_TR_10kb_coverage.txt",
                     header=FALSE)
TR_input<-read.table("Tdi_input_testes_1_TR_10kb_coverage.txt",
                     header=FALSE)

TR<-cbind(TR_cenh3,TR_input[4])
colnames(TR)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
TR$genomic_feature<-"TR"



#exons
exon_cenh3<-read.table("Tdi_cenh3_testes_1_exon_10kb_coverage.txt",
                       header=FALSE)
exon_input<-read.table("Tdi_input_testes_1_exon_10kb_coverage.txt",
                       header=FALSE)

exon<-cbind(exon_cenh3,exon_input[4])
colnames(exon)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
exon$genomic_feature<-"exon"


#introns
intron_cenh3<-read.table("Tdi_cenh3_testes_1_intron_10kb_coverage.txt",
                         header=FALSE)
intron_input<-read.table("Tdi_input_testes_1_intron_10kb_coverage.txt",
                         header=FALSE)

intron<-cbind(intron_cenh3,intron_input[4])
colnames(intron)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
intron$genomic_feature<-"intron"


#5-UTRs
UTR5_cenh3<-read.table("Tdi_cenh3_testes_1_5-UTR_10kb_coverage.txt",
                       header=FALSE)
UTR5_input<-read.table("Tdi_input_testes_1_5-UTR_10kb_coverage.txt",
                       header=FALSE)

UTR5<-cbind(UTR5_cenh3,UTR5_input[4])
colnames(UTR5)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
UTR5$genomic_feature<-"UTR5"



#3-UTRs
UTR3_cenh3<-read.table("Tdi_cenh3_testes_1_3-UTR_10kb_coverage.txt",
                       header=FALSE)
UTR3_input<-read.table("Tdi_input_testes_1_3-UTR_10kb_coverage.txt",
                       header=FALSE)

UTR3<-cbind(UTR3_cenh3,UTR3_input[4])
colnames(UTR3)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
UTR3$genomic_feature<-"UTR3"



#LINEs
LINE_cenh3<-read.table("Tdi_cenh3_testes_1_TE_LINE_10kb_coverage.txt",
                       header=FALSE)
LINE_input<-read.table("Tdi_input_testes_1_TE_LINE_10kb_coverage.txt",
                       header=FALSE)

LINE<-cbind(LINE_cenh3,LINE_input[4])
colnames(LINE)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
LINE$genomic_feature<-"LINE"



#SINEs
SINE_cenh3<-read.table("Tdi_cenh3_testes_1_TE_SINE_10kb_coverage.txt",
                       header=FALSE)
SINE_input<-read.table("Tdi_input_testes_1_TE_SINE_10kb_coverage.txt",
                       header=FALSE)

SINE<-cbind(SINE_cenh3,SINE_input[4])
colnames(SINE)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
SINE$genomic_feature<-"SINE"



#LTRs
LTR_cenh3<-read.table("Tdi_cenh3_testes_1_TE_LTR_10kb_coverage.txt",
                      header=FALSE)
LTR_input<-read.table("Tdi_input_testes_1_TE_LTR_10kb_coverage.txt",
                      header=FALSE)

LTR<-cbind(LTR_cenh3,LTR_input[4])
colnames(LTR)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
LTR$genomic_feature<-"LTR"



#RCs
RC_cenh3<-read.table("Tdi_cenh3_testes_1_TE_RC_10kb_coverage.txt",
                     header=FALSE)
RC_input<-read.table("Tdi_input_testes_1_TE_RC_10kb_coverage.txt",
                     header=FALSE)

RC<-cbind(RC_cenh3,RC_input[4])
colnames(RC)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
RC$genomic_feature<-"RC"



#DNAs
DNA_cenh3<-read.table("Tdi_cenh3_testes_1_TE_DNA_10kb_coverage.txt",
                      header=FALSE)
DNA_input<-read.table("Tdi_input_testes_1_TE_DNA_10kb_coverage.txt",
                      header=FALSE)

DNA<-cbind(DNA_cenh3,DNA_input[4])
colnames(DNA)<-c("Scaffold_name", "start", "end","coverage_cenh3", "coverage_input")
DNA$genomic_feature<-"DNA"



###Merge databases

data<-do.call("rbind", list(exon, intron, UTR5, UTR3, TR, LINE, SINE, LTR, RC, DNA))
data$coverage_cenh3_norm<-data$coverage_cenh3/69361547 #normalized by the number of mapped reads
data$coverage_input_norm<-data$coverage_input/75167800 #normalized by the number of mapped reads
data$ratiocenH3_norm<-data$coverage_cenh3_norm/data$coverage_input_norm
data$log2cenH3_norm<-log2(data$ratiocenH3_norm)

###Enriched genomic features
data<-data %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


enriched<-subset(data, log2cenH3_norm >= 2) #only select sequence categories with log2(ChIP/input) >= 2)
enriched$region<-"enriched"


###Annotation sequence categories Genome-Wide

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/gene_annotation_timema_v2")

exon_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_exons.bed",
                     header=FALSE)
colnames(exon_ann)<-c("Scaffold_name","start","end")
exon_ann$region<-"genome"
exon_ann$genomic_feature<-"exon"


intron_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_introns.bed",
                       header=FALSE)
colnames(intron_ann)<-c("Scaffold_name","start","end")
intron_ann$region<-"genome"
intron_ann$genomic_feature<-"intron"


UTR5_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_5UTRs.bed",
                     header=FALSE)
colnames(UTR5_ann)<-c("Scaffold_name","start","end")
UTR5_ann$region<-"genome"
UTR5_ann$genomic_feature<-"UTR5"


UTR3_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_3UTRs.bed",
                     header=FALSE)
colnames(UTR3_ann)<-c("Scaffold_name","start","end")
UTR3_ann$region<-"genome"
UTR3_ann$genomic_feature<-"UTR3"


ncRNA_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_ncRNA.bed",
                      header=FALSE)
colnames(ncRNA_ann)<-c("Scaffold_name","start","end")
ncRNA_ann$region<-"genome"
ncRNA_ann$genomic_feature<-"ncRNA"


rRNA_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_rRNA.bed",
                     header=FALSE)
colnames(rRNA_ann)<-c("Scaffold_name","start","end")
rRNA_ann$region<-"genome"
rRNA_ann$genomic_feature<-"rRNA"


tRNA_ann<-read.table("Tdi_LRv5a_mtDNAv350_v2.1_tRNA.bed",
                     header=FALSE)
colnames(tRNA_ann)<-c("Scaffold_name","start","end")
tRNA_ann$region<-"genome"
tRNA_ann$genomic_feature<-"tRNA"

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/TR_annotation_timema")
TR_ann<-read.table("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed",
                   header=FALSE)
colnames(TR_ann)<-c("Scaffold_name","start","end")
TR_ann$region<-"genome"
TR_ann$genomic_feature<-"TR"

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/TE_annotation_timema")
LTR_ann<-read.table("Tdi_LTRs.bed",
                    header=FALSE)
colnames(LTR_ann)<-c("Scaffold_name","start","end")
LTR_ann$region<-"genome"
LTR_ann$genomic_feature<-"LTR"


DNA_ann<-read.table("Tdi_DNAs.bed",
                    header=FALSE)
colnames(DNA_ann)<-c("Scaffold_name","start","end")
DNA_ann$region<-"genome"
DNA_ann$genomic_feature<-"DNA"


LINE_ann<-read.table("Tdi_LINEs.bed",
                     header=FALSE)
colnames(LINE_ann)<-c("Scaffold_name","start","end")
LINE_ann$region<-"genome"
LINE_ann$genomic_feature<-"LINE"


SINE_ann<-read.table("Tdi_SINEs.bed",
                     header=FALSE)
colnames(SINE_ann)<-c("Scaffold_name","start","end")
SINE_ann$region<-"genome"
SINE_ann$genomic_feature<-"SINE"


RC_ann<-read.table("Tdi_RCs.bed",
                   header=FALSE)
colnames(RC_ann)<-c("Scaffold_name","start","end")
RC_ann$region<-"genome"
RC_ann$genomic_feature<-"RC"

data2<-rbind(do.call("rbind", list(rRNA_ann, tRNA_ann, ncRNA_ann, exon_ann, intron_ann, UTR5_ann, UTR3_ann, TR_ann[,c(1,2,3,5,6)], LINE_ann, SINE_ann, LTR_ann, RC_ann, DNA_ann, enriched[,c(1,2,3,6,11)])))
data2$array<-data2$end-data2$start
data2$region_genomic_feature<-paste(data2$region, data2$genomic_feature, sep="-")
data3<-aggregate(data2$array~data2$region_genomic_feature, FUN=sum)
data3[c('region', 'genomic_feature')] <- str_split_fixed(data3$`data2$region_genomic_feature`, '-', 2)


data4 <- data3 %>%
  group_by(region) %>%                     # grouped by rep
  mutate(sum_rep=sum(`data2$array`)) %>%          # sum of each rep
  group_by(region,genomic_feature) %>%            # grouped by DB
  summarise(perc=sum(`data2$array`)/unique(sum_rep))


ggplot(data4, aes(x = region, y = perc, fill = genomic_feature)) +
  geom_col(color="black", size=0.1) +scale_x_discrete(limits = c("genome",  "enriched"))+
  scale_fill_manual(values = c("#E2D1F9", "#2F3C7E", "#FCF6F5", "#8A307F", "#79A7D3", "white","#6883BC", "#EA738D", "#DDC3A5" , "#F7C5CC", "#990011", "#101820", "#E0A96D"))+
  theme_classic()



##only for autosomes
data2<-rbind(do.call("rbind", list(rRNA_ann, tRNA_ann, ncRNA_ann, exon_ann, intron_ann, UTR5_ann, UTR3_ann, TR_ann[,c(1,2,3,5,6)], LINE_ann, SINE_ann, LTR_ann, RC_ann, DNA_ann, enriched[,c(1,2,3,6,11)])))
data2$array<-data2$end-data2$start
data2A<-subset(data2, Scaffold_name=="Tdi_LRv5a_scf1" | Scaffold_name=="Tdi_LRv5a_scf2" | Scaffold_name=="Tdi_LRv5a_scf4" | Scaffold_name=="Tdi_LRv5a_scf5" |
                 Scaffold_name=="Tdi_LRv5a_scf6" | Scaffold_name=="Tdi_LRv5a_scf7" | Scaffold_name=="Tdi_LRv5a_scf8" | Scaffold_name=="Tdi_LRv5a_scf9" | Scaffold_name=="Tdi_LRv5a_scf10" |
                 Scaffold_name=="Tdi_LRv5a_scf11" | Scaffold_name=="Tdi_LRv5a_scf12")
data2A$region_genomic_feature<-paste(data2A$region, data2A$genomic_feature, sep="-")
data3<-aggregate(data2A$array~data2A$region_genomic_feature, FUN=sum)
data3[c('region', 'genomic_feature')] <- str_split_fixed(data3$`data2A$region_genomic_feature`, '-', 2)


data4 <- data3 %>%
  group_by(region) %>%                     # grouped by rep
  mutate(sum_rep=sum(`data2A$array`)) %>%          # sum of each rep
  group_by(region,genomic_feature) %>%            # grouped by DB
  summarise(perc=sum(`data2A$array`)/unique(sum_rep))


ggplot(data4, aes(x = region, y = perc, fill = genomic_feature)) +
  geom_col(color="black", size=0.1) +scale_x_discrete(limits = c("genome",  "enriched"))+
  scale_fill_manual(values = c("#E2D1F9", "#2F3C7E", "#FCF6F5", "#8A307F", "#79A7D3", "white","#6883BC", "#EA738D", "#DDC3A5" , "#F7C5CC", "#990011", "#101820", "#E0A96D"))+
  theme_classic()


##Statistics

table<-xtabs(`data2A$array` ~ region+genomic_feature, data=data3)
table<-xtabs(`data2$array` ~ region+genomic_feature, data=data3)

chisq.test(table)
#chisq.test(table[c(1,3,5,7,9,11,13,15,17,19,21,23,25)], p = data4$perc[data4$region=="genome"])
#The p-value of the test is < 2.2e-16, which is lower than the significance level alpha = 0.05.
#We can conclude that the observed proportions are significantly different.
