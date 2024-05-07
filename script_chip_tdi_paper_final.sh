#!/bin/bash


#Slurm options:

#SBATCH -p cpu
#SBATCH --time=0-1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB


#Commands start here:

echo "Starting job $SLURM_JOB_NAME with ID $SLURM_JOB_ID".

module load gcc
module load bwa
module load perl
module load bedtools2
module load samtools
module load python


#######################################################
## Genome annotation of specific sequence categories ##
#######################################################

  ## Genic features
#Refer to ...

## Tandem Repeats (TRs)

#Refer to script "script_tandem-repeat_annotation.sh"

  ## Transposable Elements (TEs)
#Refer to ...


########################
## Prepare chip reads ##
########################

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/

cat *cenh3*R1*.gz > Tdi_testes_cenh3_1_R1.fq.gz
cat *cenh3*R3*.gz > Tdi_testes_cenh3_1_R3.fq.gz

cat *input*R1*.gz > Tdi_testes_input_1_R1.fq.gz
cat *input*R3*.gz > Tdi_testes_input_1_R3.fq.gz


################
## Trim reads ##
################

module load trimmomatic

for i in *_R1.fq.gz ; do
        foo1=`echo $i`
                basename=`echo $foo1 | sed 's/_R1.fq.gz*//' | sed 's/.*\///'`
        infileR1=`echo $foo1`
        infileR2=`echo $foo1 | sed 's/_R1.fq.gz/_R3.fq.gz/'`
        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq"`
        outfileR2=`echo "./"$basename"_R3_qtrimmed.fq"`
        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq"`
        outfileR2_UP=`echo "./"$basename"_R3_qtrimmed_UNPAIRED.fq"`

#        echo $infileR1
#        echo $infileR2
#        echo $outfileR1
#        echo $outfileR1_UP
#        echo $outfileR2
#        echo $outfileR2_UP

        trimmomatic PE -threads 16 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
done


#############################
## Prepare genome assembly ##
#############################

module load bedtools2
module load samtools

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/genomes

samtools faidx Tdi_LRv5a_mtDNAv350.fasta

cut -f1,2 Tdi_LRv5a_mtDNAv350.fasta.fai > Tdi_chm_size_mtDNAv350.txt

bedtools makewindows -g Tdi_chm_size_mtDNAv350.txt -w 10000  > Tdi_chm_size_mtDNAv350_w10000.bed


###################
## Mapping reads ##
###################

	## bwa index

module load bwa
module load samtools

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip

bwa index genomes/Tdi_LRv5a_mtDNAv350.fasta


	## bwa map
mkdir bwa

bwa mem -t 16 -c 1000000000 -T 30 genomes/Tdi_LRv5a_mtDNAv350.fasta Tdi_testes_cenh3_1_R1_qtrimmed.fq Tdi_testes_cenh3_1_R3_qtrimmed.fq > bwa/chip_cenh3_tdi_testes_1_bwa.sam
bwa mem -t 16 -c 1000000000 -T 30 genomes/Tdi_LRv5a_mtDNAv350.fasta chip/Tdi_testes_input_1_R1_qtrimmed.fq chip/Tdi_testes_input_1_R3_qtrimmed.fq > bwa/chip_input_tdi_testes_1_bwa.sam


	## flagstat

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/bwa

samtools flagstat chip_cenh3_tdi_testes_1_bwa.sam > chip_cenh3_tdi_testes_1_bwa_flagstat.txt
samtools flagstat chip_input_tdi_testes_1_bwa.sam > chip_input_tdi_testes_1_bwa_flagstat.txt


	## sort bam

samtools view -u chip_cenh3_tdi_testes_1_bwa.sam | samtools sort -o chip_cenh3_tdi_testes_1_bwa.bam
samtools view -u chip_input_tdi_testes_1_bwa.sam | samtools sort -o chip_input_tdi_testes_1_bwa.bam

#rm chip_cenh3_tdi_testes_1_bwa.sam
#rm chip_input_tdi_testes_1_bwa.sam


	## remove supp (chimeric) alignements

samtools view chip_cenh3_tdi_testes_1_bwa.bam | fgrep SA:Z: | cut -f 1 > chip_cenh3_tdi_testes_1_bwa_badnames.txt
samtools view -h chip_cenh3_tdi_testes_1_bwa.bam | fgrep -vf chip_cenh3_tdi_testes_1_bwa_badnames.txt | samtools view -b > chip_cenh3_tdi_testes_1_bwa_final.bam
#rm chip_cenh3_tdi_testes_1_bwa_badnames.txt
samtools flagstat chip_cenh3_tdi_testes_1_bwa_final.bam > chip_cenh3_tdi_testes_1_bwa_final_flagstat.txt


samtools view chip_input_tdi_testes_1_bwa.bam | fgrep SA:Z: | cut -f 1 > chip_input_tdi_testes_1_bwa_badnames.txt
samtools view -h chip_input_tdi_testes_1_bwa.bam | fgrep -vf chip_input_tdi_testes_1_bwa_badnames.txt | samtools view -b > chip_input_tdi_testes_1_bwa_final.bam
#rm chip_input_tdi_testes_1_bwa_badnames.txt
samtools flagstat chip_input_tdi_testes_1_bwa_final.bam > chip_input_tdi_testes_1_bwa_final_flagstat.txt


	##remove PCR duplicates
module load picard
module load samtools

for i in  *tdi_testes_1*.bam; do
    outbam=`echo $i | sed 's/_bwa_final.bam/_bwa_final_DR.bam/'`
       flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
       metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

#       echo $i
#       echo $outbam
#       echo $metric_file
#       echo $flagstat_out_bam

       picard MarkDuplicates REMOVE_DUPLICATES=true \
       INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file

       samtools flagstat $outbam > $flagstat_out_bam
       mv $flagstat_out_bam /scratch/wtoubian/bwa_timema_genomes/mapping_genomes/BWA_out/flagstat_out_paired

done


########################################
## Calculate coverage GW (10 kb bins) ##
########################################

module load bedtools2

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip
mkdir coverage

bedtools coverage -a genomes/Tdi_chm_size_mtDNAv350_w10000.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/Tdi_cenh3_testes_1_GW_coverage_DR.txt
bedtools coverage -a genomes/Tdi_chm_size_mtDNAv350_w10000.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/Tdi_input_testes_1_GW_coverage_DR.txt


###########################################################
## Prepare bed files for specific sequence categories GW ##
###########################################################

module load bedtools2

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip

	## Tandem Repeats
sortBed -faidx genomes/Tdi_chm_size_mtDNAv350.txt -i TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed > TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted.bed

awk '{print $1 "\t" $2 "\t" $3;}' TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted.bed > TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_noSeq.bed

	## Genic features

#Add introns to annotation of genic features
#We used a script from AGAT toolkit: Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717.

perl agat_sp_add_introns.pl --gff gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1.gff --out gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns.gff #agat_sp_add_introns.pl can be find here: https://github.com/NBISweden/AGAT/blob/master/bin/agat_sp_add_introns.pl
awk '{print $1 "\t" $4 "\t" $5 "\t" $3;}' gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns.gff > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns.bed
sortBed -faidx genomes/Tdi_chm_size_mtDNAv350.txt -i gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns.bed > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed

#exons
grep -w "exon" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_exons.bed

#UTRs
grep -w "5'-UTR" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_5UTRs.bed
grep -w "3'-UTR" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_3UTRs.bed

#introns
grep -w "intron" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_introns.bed

#ncRNA
grep -w "ncRNA" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_ncRNA.bed

#rRNA
grep -w "rRNA" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_rRNA.bed

#tRNA
grep -w "tRNA" gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_add_introns_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > gene_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_tRNA.bed

	## TEs
sortBed -faidx genomes/Tdi_chm_size_mtDNAv350.txt -i TE_annotation_timema/Tdi_AllRepeats.classi.bed > TE_annotation_timema/Tdi_AllRepeats.classi_sorted.bed
grep -w "LINE" TE_annotation_timema/Tdi_AllRepeats.classi_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > TE_annotation_timema/Tdi_LINEs.bed
grep -w "SINE" TE_annotation_timema/Tdi_AllRepeats.classi_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > TE_annotation_timema/Tdi_SINEs.bed
grep -w "DNA" TE_annotation_timema/Tdi_AllRepeats.classi_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > TE_annotation_timema/Tdi_DNAs.bed
grep -w "LTR" TE_annotation_timema/Tdi_AllRepeats.classi_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > TE_annotation_timema/Tdi_LTRs.bed
grep -w "RC" TE_annotation_timema/Tdi_AllRepeats.classi_sorted.bed | awk '{print $1 "\t" $2 "\t" $3;}' > TE_annotation_timema/Tdi_RCs.bed


############################################################
## Calculate coverage for specific sequence categories GW ##
############################################################

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip
mkdir coverage/GW

## TR coverage
bedtools coverage -a TR_annotation_timema/Tdi_LRv5a/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_noSeq.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam  -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_TR_coverage_GW.txt
bedtools coverage -a TR_annotation_timema/Tdi_LRv5a/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_noSeq.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam  -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_TR_coverage_GW.txt

## TE coverage
#LINEs
bedtools coverage -a TE_annotation_timema/Tdi_LINEs.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_TE_LINE_coverage_GW.txt
bedtools coverage -a TE_annotation_timema/Tdi_LINEs.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_TE_LINE_coverage_GW.txt

#SINEs
bedtools coverage -a TE_annotation_timema/Tdi_SINEs.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_TE_SINE_coverage_GW.txt
bedtools coverage -a TE_annotation_timema/Tdi_SINEs.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_TE_SINE_coverage_GW.txt

#LTRs
bedtools coverage -a TE_annotation_timema/Tdi_LTRs.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_TE_LTR_coverage_GW.txt
bedtools coverage -a TE_annotation_timema/Tdi_LTRs.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_TE_LTR_coverage_GW.txt

#RCs
bedtools coverage -a TE_annotation_timema/Tdi_RCs.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_TE_RC_coverage_GW.txt
bedtools coverage -a TE_annotation_timema/Tdi_RCs.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_TE_RC_coverage_GW.txt

#DNAs
bedtools coverage -a TE_annotation_timema/Tdi_DNAs.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_TE_DNA_coverage_GW.txt
bedtools coverage -a TE_annotation_timema/Tdi_DNAs.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_TE_DNA_coverage_GW.txt

## exons coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_exons.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_exon_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_exons.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_exon_coverage_GW.txt

## introns coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_introns.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_intron_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_introns.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_intron_coverage_GW.txt

        ## 5-UTR coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_5UTRs.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_5-UTR_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_5UTRs.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_5-UTR_coverage_GW.txt

        ## 3-UTR coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_3UTRs.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_3-UTR_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_3UTRs.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_3-UTR_coverage_GW.txt

        ## ncRNA coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_ncRNA.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_ncRNA_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_ncRNA.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_ncRNA_coverage_GW.txt

        ## rRNA coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_rRNA.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_rRNA_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_rRNA.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_rRNA_coverage_GW.txt

        ## tRNA coverage
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_tRNA.bed -b bwa/chip_input_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_input_testes_1_tRNA_coverage_GW.txt
bedtools coverage -a genome_annotation_timema_v2/Tdi_LRv5a_mtDNAv350_v2.1_tRNA.bed -b bwa/chip_cenh3_tdi_testes_bwa_final_DR.bam -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/GW/Tdi_cenh3_testes_1_tRNA_coverage_GW.txt


##################################################################
## Estimate proportion and enrichment of sequence categories GW ##
##################################################################

module load r

./Proportion_categories.R


#################################################################
## Select enriched 10kb windows based on log2(ChIP/Input) >= 2 ##
#################################################################

module load r
mkdir enriched_10kb_regions

./Enriched_windows.R
#tdi_cenh3_testes_GW_coverage_DR_w10000_logRatio2.bed file was subsequently generated, comprising 10kb windows with log2(ChIP/Input) >= 2


#######################################################################
## Prepare bed files for Tandem Repeats within enriched 10kb windows ##
#######################################################################

mkdir enriched_10kb_regions/categories

	#TRs
bedtools intersect -a enriched_10kb_regions/tdi_cenh3_testes_GW_coverage_DR_w10000_logRatio2.bed -b TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_noSeq.bed > enriched_10kb_regions/categories/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_10kb.bed

sortBed -faidx genomes/Tdi_chm_size_mtDNAv350.txt -i enriched_10kb_regions/categories/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_10kb.bed > enriched_10kb_regions/categories/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_10kb_sorted.bed


#######################################################################
## Calculate coverage on Tandem Repeats within enriched 10kb windows ##
#######################################################################

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip
mkdir coverage/enriched_windows

        ## TR coverage
bedtools coverage -a enriched_10kb_regions/categories/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_10kb_sorted.bed -b bwa/chip_input_tdi_testes_1_bwa_final_DR.bam  -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/enriched_windows/Tdi_input_testes_1_TR_10kb_coverage.txt
bedtools coverage -a enriched_10kb_regions/categories/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_10kb_sorted.bed -b bwa/chip_cenh3_tdi_testes_1_bwa_final_DR.bam  -sorted -g genomes/Tdi_LRv5a_mtDNAv350.fasta.fai -mean > coverage/enriched_windows/Tdi_cenh3_testes_1_TR_10kb_coverage.txt


#############################################################################
## Extract TR motifs (with minimal rotation) within enriched 10 kb windows ##
#############################################################################

##TRF parsed .txt file was further processed to compute minimal rotations on motif sequences from parse TRF output
perl script_minimal_rotation_parse.pl TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt > TR_annotation_timema/Tdi_LRv5a/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt

awk '{print $1 "\t" $2 "\t" $3 "\t" $7;}' TR_annotation_timema/Tdi_LRv5a/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt | awk 'NR>1' > TR_annotation_timema/Tdi_LRv5a/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.bed


##Intersect .bed file with minimal rotation sequences and .bed file of enriched 10kb regions
#bedtools version (2.30) on the cluster gives truncated TR sequence (below 1000bp) so I download an updated version (2.31) on my local
source ~/.bashrc
conda activate /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/softwares/bedtools_env
bedtools intersect -a TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.bed -b enriched_10kb_regions/tdi_cenh3_testes_GW_coverage_DR_w10000_logRatio2.bed > enriched_10kb_regions/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation_enriched_10kb_regions.txt


##Extract motifs with minimal rotation sequences and log2(ChIP/Input) >= 2
mkdir enriched_motifs
mkdir enriched_motifs/levenstein
mkdir enriched_motifs/blast

./Enriched_minimal_rotation_TR_motifs.R
#This script 1) filters out motifs with log2(ChIP/Input) < 2, 2) selects only motifs with distinct sequences (based on minimal rotation), 3) duplicates motif sequences until reaching 2000bp in length, 4) generates reverse complements
#Enriched_2kbmotifs_enriched_10kb_regions.txt file was subsequently generated for Levenstein distance analysis
#Enriched_2kbmotifs_reverseComp_enriched_10kb_regions.txt file was subsequently generated for Levenstein distance analysis
#Enriched_motifs_database.txt file was subsequently generated for Levenstein distance network analysis
#Enriched_motifs_enriched_10kb_regions.fasta file was subsequently generated for Blast analysis


####################################
## Calculate Levenstein distances ##
####################################

module load python

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/levenstein

python script_levenshtein_rotations_F-F_inputFiles.py #computes pair-wise levenstein distances between duplicated motif sequences. Outputs levenstein_distances_FF.txt file
python script_levenshtein_rotations_F-R_inputFiles.py #computes pair-wise levenstein distances between duplicated motif sequences and reverse complements. Outputs levenstein_distances_FR.txt file
#levenstein_distances_FF.txt file was subsequently generated for Levenstein distance network analysis
#levenstein_distances_FR.txt file was subsequently generated for Levenstein distance network analysis


#################################
## Levenstein distance network ##
#################################

module load r

./Levenstein_network.R
#Network_database.txt file was subsequently generated for Heatmap analyses
#Supplementary_table2.txt was subsequently generated


##########################
## Heatmap based on TRF ##
##########################

module load r

./Heatmap_TRF-based.R


########################################################################
## blast enriched motifs with distinct sequences onto genome assembly ##
########################################################################

module load blast-plus

##Generate fasta of enriched 10kb regions
bedtools getfasta -fi genomes/Tdi_LRv5a_mtDNAv350.fasta -bed enriched_10kb_regions/tdi_cenh3_testes_GW_coverage_DR_w10000_logRatio2.bed > enriched_motifs/blast/tdi_cenh3_testes_GW_bwa_coverage_RS_logRatio2.fasta

##Blast database of  enriched 10kb regions
makeblastdb -in enriched_motifs/blast/tdi_cenh3_testes_GW_bwa_coverage_RS_logRatio2.fasta -dbtype nucl -parse_seqids -max_file_sz '3.9GB' -input_type fasta

##blastn
blastn -task blastn-short -query enriched_motifs/blast/Enriched_motifs_enriched_10kb_regions.fasta -db enriched_motifs/blast/tdi_cenh3_testes_GW_bwa_coverage_RS_logRatio2.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore" -num_threads 16 -out enriched_motifs/blast/Enriched_motifs_enriched_10kb_regions_blastn.txt


############################
## Heatmap based on blast ##
############################

module load r

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/enriched_motifs/blast

./Heatmap_blast-based.R
#Repeat_families_10kb_regions_blastn.bed was subsequently generated for visualizing blast coordinates within Geneious Prime





######################################################################################
## k-mer analysis (adapted from https://github.com/straightlab/xenla-cen-dna-paper) ##
######################################################################################

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/
mkdir kmer_centromere

source ~/.bashrc
conda activate /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/kmc_env/
#snakemake was conda intalled in the kmc_env beforehand. this avoid issues with pandas module import in python.

snakemake -s kmer_centromere/snakefile_v9_genomev10PRE_SE.py --configfile kmer_centromere/config_xla_merge_final_SE_v9_yf_genome10.2.yaml -pr --cores 8


##################################################################
## De novo assembly CenH3-ChIP reads containing enriched k-mers ##
##################################################################

module load gcc spades

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/
mkdir kmer_centromere/spades_assembly

spades.py --careful --only-assembler -s kmer_centromere/result_tdi_cenh3_testes_GW_kmer25_CIVAL100/yf_seqdat/xla_genome_chunks_v10.2/xla_merge_final/25/100/madx25/CAoINPUT_extract.fa -o kmer_centromere/spades_assembly/


####################################################################
## Extract TR motifs (with minimal rotation) from de novo contigs ##
####################################################################

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/spades_assembly
mkdir levenstein

./kmer_minimal_rotation_TR_motifs.R
#Unique_TR_sequences_contigs.txt file was subsequently generated for Levenstein distance analysis
#Unique_TR_sequences_contigs_reverseComp.txt file was subsequently generated for Levenstein distance analysis
#Unique_TR_sequences_contigs_database.txt file was subsequently generated for Levenstein distance network analysis


####################################
## Calculate Levenstein distances ##
####################################

module load python

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/kmer_centromere/spades_assembly/levenstein

#Here the same scripts as for the genome assembly were used. Only input files were modified and output files were renamed.
python script_levenshtein_rotations_F-F_inputFiles.py #computes pair-wise levenstein distances between duplicated motif sequences. Outputs levenstein_distances_FF_contigs.txt file
python script_levenshtein_rotations_F-R_inputFiles.py #computes pair-wise levenstein distances between duplicated motif sequences and reverse complements. Outputs levenstein_distances_FR_contigs.txt file
#levenstein_distances_FF_contigs.txt file was subsequently generated for Levenstein distance network analysis
#levenstein_distances_FR_contigs.txt file was subsequently generated for Levenstein distance network analysis


#################################
## Levenstein distance network ##
#################################

module load r

./Levenstein_network_contigs.R
#Supplementary_table3.txt was subsequently generated


##END##
