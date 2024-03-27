#!/usr/bin/env Rscript

library(ggplot2)

###Select enriched 10kb regions
setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/coverage/")

cenh3<-read.table("Tdi_cenh3_testes_1_GW_coverage_DR_w10000.txt",
                  header=FALSE)
input1<-read.table("Tdi_input_testes_1_GW_coverage_DR_w10000.txt",
                   header=FALSE)

data1<-cbind(cenh3,input1[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1$coverage_cenh3_norm<-data1$coverage_cenh3/69361547 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/75167800 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratiocenH3_norm<-data1$coverage_cenh3_norm/data1$coverage_input_norm
data1$log2cenH3_norm<-log2(data1$ratiocenH3_norm)


subset<-subset(data1, log2cenH3_norm >= 2)
subset<-subset %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values

setwd("/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_10kb_regions")
write.table(subset[,1:3], file = "tdi_cenh3_testes_GW_coverage_DR_w10000_logRatio2.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)



###Plots

chms<-subset(data1, Scaffold_name=="Tdi_LRv5a_scf1" | Scaffold_name=="Tdi_LRv5a_scf2" | Scaffold_name=="Tdi_LRv5a_scf3" | Scaffold_name=="Tdi_LRv5a_scf4" | Scaffold_name=="Tdi_LRv5a_scf5" 
             | Scaffold_name=="Tdi_LRv5a_scf6" | Scaffold_name=="Tdi_LRv5a_scf7" | Scaffold_name=="Tdi_LRv5a_scf8" | Scaffold_name=="Tdi_LRv5a_scf9" | Scaffold_name=="Tdi_LRv5a_scf10"
             | Scaffold_name=="Tdi_LRv5a_scf11" | Scaffold_name=="Tdi_LRv5a_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$Scaffold_name)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms2<-chms[order(chms$scaff_number, chms$start),]
chms2$cumul_start<-seq(0, by = 10000, length.out = nrow(chms2))
chms2$chrom<-ifelse(chms2$scaff_number=="3", "sex-chrom", "autosome")
chms3<-chms2[chms2$scaff_number=="3",]

plot_chip<-ggplot(chms2, aes(x=cumul_start, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=1) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms2$scaff_number)))/2))+
  scale_y_continuous(limits=c(-3,6)) +
  geom_point(data=chms3, aes(x=chms3$cumul_start, y=chms3$log2cenH3_norm), stat="identity", colour="maroon", size=1) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept=c(2,0,0), linetype=c("dashed","solid", "dashed"), color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )

