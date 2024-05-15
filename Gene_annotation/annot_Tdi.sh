## 2a_run_annot.sh

DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot'
SDIR='/scratch/dparker/annot'
genome_pref="Tdi_LRv5a"

#############################################################################################
### get genomes

cd $DIR
mkdir genomes
cp /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/fasta/*fasta genomes

#############################################################################################
### Repeat identification with Repeat Modeler through Dfam TE Tools # Container v1.4 (https://github.com/Dfam-consortium/TETools)

## Info
# https://www.repeatmasker.org/RepeatModeler/

## Brůna T, Hoff KJ, Lomsadze A, Stanke M, Borodovsky M (2021). BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database. NAR Genom Bioinform 3: lqaa108.
### https://academic.oup.com/nargab/article/3/1/lqaa108/6066535#supplementary-data
#Repeat masking by RepeatModeler and RepeatMasker with default settings was sufficient to achieve high gene prediction accuracy in all the tested genomes except for X. tropicalis.
#BuildDatabase -engine wublast -name genome genome.fasta
#RepeatModeler -engine wublast -database genome
#RepeatMasker -engine wublast -lib genome-families.fa -xsmall genome.fasta

## maker does one round of a specif lib, then a general one from model taxa.
## https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2 ### example with boa

## calc in the /scratch

cd $SDIR

module load singularity/3.7.4
#singularity pull docker://dfam/tetools

## build database for rep modeler
singularity exec --bind $DIR,$SDIR $DIR/tetools_latest.sif BuildDatabase -engine ncbi -name $genome_pref $DIR/genomes/$genome_pref.fasta


## /opt/RepeatModeler/RepeatModeler - 2.0.2
## run rep modeler to make sp rep lib # 32 cpu, 80GB ram ## took: 2-02:22:38
singularity exec --bind $DIR,$SDIR $DIR/tetools_latest.sif RepeatModeler -engine ncbi -database $genome_pref -pa 30 -LTRStruct


#RepeatMasker version 4.1.2-p1
#No query sequence file indicated
#/opt/RepeatMasker/RepeatMasker - 4.1.2-p1
## run rep masker # 32 cpu, 80GB ram ## took 07:01:29

genome_file=`echo $DIR/genomes/$genome_pref".fasta"`
singularity exec --bind $DIR,$SDIR $DIR/tetools_latest.sif RepeatMasker -engine ncbi -gff -xsmall -pa 30 -lib $SDIR/$genome_pref"-families.fa" $genome_file

## cp stats
cp  $DIR/genomes/*tbl $DIR/Timema_LR_genomic_code/output/repmask/



######################################################################################################################
########## BRAKER 2 run
## Braker is stocastic due to augustus. Gonna rune each 5 times and take the best.
## prots from 2b_run_annot_Tpa.sh
# arth_proteins.fasta 
# arth_and_timemav8_proteins.fasta
# braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif available here: https://doi.org/10.5281/zenodo.11196713

#######################################################################################################################################
### prot runs | arth+Timema |  5 times each - take best

###### 3a | Run time 12:03:04
mkdir $SDIR/Tdi_braker_prot_run_3a
cd    $SDIR/Tdi_braker_prot_run_3a
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tdi sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_prot_run_3a --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_prot_run_3a/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi


###### 3b | Run time 
mkdir $SDIR/Tdi_braker_prot_run_3b
cd    $SDIR/Tdi_braker_prot_run_3b
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tdi sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_prot_run_3b --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_prot_run_3b/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi


###### 3c | Run time 
mkdir $SDIR/Tdi_braker_prot_run_3c
cd    $SDIR/Tdi_braker_prot_run_3c
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tdi sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_prot_run_3c --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_prot_run_3c/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi


###### 3d | Run time 
mkdir $SDIR/Tdi_braker_prot_run_3d
cd    $SDIR/Tdi_braker_prot_run_3d
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tdi sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_prot_run_3d --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_prot_run_3d/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi


###### 3e | Run time 
mkdir $SDIR/Tdi_braker_prot_run_3e
cd    $SDIR/Tdi_braker_prot_run_3e
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tdi sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_prot_run_3e --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_prot_run_3e/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi



## cp back from scratch

cp -r $SDIR/Tdi_braker_prot_run_3a $DIR
cp -r $SDIR/Tdi_braker_prot_run_3b $DIR
cp -r $SDIR/Tdi_braker_prot_run_3c $DIR
cp -r $SDIR/Tdi_braker_prot_run_3d $DIR
cp -r $SDIR/Tdi_braker_prot_run_3e $DIR


awk '{print $3}' $DIR/Tdi_braker_prot_run_3a/braker.gff3 | grep "gene" | wc -l
# 47958
awk '{print $3}' $DIR/Tdi_braker_prot_run_3b/braker.gff3 | grep "gene" | wc -l
# 47596
awk '{print $3}' $DIR/Tdi_braker_prot_run_3c/braker.gff3 | grep "gene" | wc -l
# 47844
awk '{print $3}' $DIR/Tdi_braker_prot_run_3d/braker.gff3 | grep "gene" | wc -l
# 55764
awk '{print $3}' $DIR/Tdi_braker_prot_run_3e/braker.gff3 | grep "gene" | wc -l
# 50736


# tarball (all - I can untar the ones I need. 
cd $DIR


tar -czf Tdi_braker_prot_run_3a.tar.gz Tdi_braker_prot_run_3a
tar -czf Tdi_braker_prot_run_3b.tar.gz Tdi_braker_prot_run_3b
tar -czf Tdi_braker_prot_run_3c.tar.gz Tdi_braker_prot_run_3c
tar -czf Tdi_braker_prot_run_3d.tar.gz Tdi_braker_prot_run_3d
tar -czf Tdi_braker_prot_run_3e.tar.gz Tdi_braker_prot_run_3e


rm -r Tdi_braker_prot_run_3a/
rm -r Tdi_braker_prot_run_3b/
rm -r Tdi_braker_prot_run_3c/
rm -r Tdi_braker_prot_run_3d/
rm -r Tdi_braker_prot_run_3e/



#######################################################################################################
### map - STAR with --twopassMode Basic (to get the splice junction right)
### map invid then merge.
### map all from Tps and Tdi to Tdi

### 117 conditions (NOT using NA dev stages, and single-end GN, LG, HD)
### 114 paired, 3 single end
### should get 364 paired bams, 12 single bams
### See TpsTdi_Annot_RNAseq_samples.xlsx
### Reads trimmed with trimmomatic v. 0.39; options: ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:80


##########################
### Index
### Build index for star
### not using soft masked genome here - STAR ignores soft masking anyway (tested)

module load gcc
module load star/2.7.8a

mkdir $DIR/RNAseq_mapping
cd    $DIR/RNAseq_mapping

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir $genome_pref"_STAR" \
     --genomeFastaFiles $DIR"/genomes/"$genome_pref".fasta" \
     --genomeChrBinNbits 20 --limitGenomeGenerateRAM 79000000000


## map paired end reads


Tps_Tdi_paired_samples=(
Tdi_F_A_Ad
Tdi_F_A_Ha
Tdi_F_A_J2
Tdi_F_A_J3
Tdi_F_A_J4
Tdi_F_A_J5
Tdi_F_B_Ad
Tdi_F_B_Ha
Tdi_F_B_J2
Tdi_F_B_J3
Tdi_F_B_J4
Tdi_F_B_J5
Tdi_F_DG_Ad
Tdi_F_FB_Ad
Tdi_F_Fe_Ad
Tdi_F_Go_Ad
Tdi_F_Go_J3
Tdi_F_Go_J4
Tdi_F_Go_J5
Tdi_F_Gu_Ad
Tdi_F_Gu_Ha
Tdi_F_Gu_J2
Tdi_F_Gu_J3
Tdi_F_Gu_J4
Tdi_F_Gu_J5
Tdi_F_Lg_Ad
Tdi_F_Lg_Ha
Tdi_F_Lg_J2
Tdi_F_Lg_J3
Tdi_F_Lg_J4
Tdi_F_Lg_J5
Tdi_F_Ta_Ad
Tdi_F_WB_Ha
Tdi_M_A_Ad
Tdi_M_AG_Ad
Tdi_M_B_Ad
Tdi_M_DG_Ad
Tdi_M_FB_Ad
Tdi_M_Fe_Ad
Tdi_M_Gu_Ad
Tdi_M_Ta_Ad
Tdi_M_Te_Ad
Tps_F_A_Ad
Tps_F_A_Ha
Tps_F_A_J2
Tps_F_A_J3
Tps_F_A_J4
Tps_F_A_J5
Tps_F_A_J6
Tps_F_B_Ad
Tps_F_B_Ha
Tps_F_B_J2
Tps_F_B_J3
Tps_F_B_J4
Tps_F_B_J5
Tps_F_B_J6
Tps_F_DG_Ad
Tps_F_FB_Ad
Tps_F_Fe_Ad
Tps_F_Go_Ad
Tps_F_Go_J3
Tps_F_Go_J4
Tps_F_Go_J5
Tps_F_Go_J6
Tps_F_Gu_Ad
Tps_F_Gu_Ha
Tps_F_Gu_J2
Tps_F_Gu_J3
Tps_F_Gu_J4
Tps_F_Gu_J5
Tps_F_Gu_J6
Tps_F_Lg_Ad
Tps_F_Lg_Ha
Tps_F_Lg_J2
Tps_F_Lg_J3
Tps_F_Lg_J4
Tps_F_Lg_J5
Tps_F_Lg_J6
Tps_F_Ta_Ad
Tps_M_A_Ad
Tps_M_A_Ha
Tps_M_A_J2
Tps_M_A_J3
Tps_M_A_J4
Tps_M_A_J5
Tps_M_AG_Ad
Tps_M_B_Ad
Tps_M_B_Ha
Tps_M_B_J2
Tps_M_B_J3
Tps_M_B_J4
Tps_M_B_J5
Tps_M_DG_Ad
Tps_M_FB_Ad
Tps_M_Fe_Ad
Tps_M_Go_Ad
Tps_M_Go_J3
Tps_M_Go_J4
Tps_M_Go_J5
Tps_M_Gu_Ad
Tps_M_Gu_Ha
Tps_M_Gu_J2
Tps_M_Gu_J3
Tps_M_Gu_J4
Tps_M_Gu_J5
Tps_M_Lg_Ad
Tps_M_Lg_Ha
Tps_M_Lg_J2
Tps_M_Lg_J3
Tps_M_Lg_J4
Tps_M_Lg_J5
Tps_M_Ta_Ad
Tps_M_Te_Ad
Tps_U_Ca_Ha
)

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmed_*/*_R1_*qtrimmed.fq.gz; do
R2_f=`echo $R1_f | sed 's/_R1_/_R2_/' `
out_prefix=`echo $R1_f | sed 's/.*\///' | sed 's/_R1_.*//'`
sp=`echo $R1_f | sed 's/.*\///' | sed 's/_.*//'`

for s in ${Tps_Tdi_paired_samples[@]}; do

if  [[ $out_prefix == $s* ]];
then
echo $R1_f
echo $R2_f
echo $sp
echo $genome_pref
echo $out_prefix
STAR --twopassMode Basic --genomeDir $genome_pref"_STAR" \
     --readFilesIn $R1_f $R2_f --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN 12 \
     --outFileNamePrefix  $out_prefix"_to_"$genome_pref
fi
done
done 




### single end

mkdir $SDIR/RNAseq_mapping

Tps_Tdi_single_samples=(
Tdi_F_WB_Ad
Tps_F_WB_Ad
Tps_M_WB_Ad
)

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmedSE*/*.fq.gz; do
out_prefix=`echo $R1_f | sed 's/.*\///' | sed 's/_AandQtrimmed.*//' | sed 's/_trimmed.*//'  | sed 's/sampled.*//' | sed 's/_DJP.*//' | sed 's/_Swb.*//'  | sed 's/_Eth.*//' | sed 's/_S2e.*//'  `
sp=`echo $R1_f | sed 's/.*\///' | sed 's/_.*//'`

for s in ${Tps_Tdi_single_samples[@]}; do

if  [[ $out_prefix == $s* ]];
then
echo $R1_f
echo $sp
echo $genome_pref
echo $out_prefix
STAR --twopassMode Basic --genomeDir $SDIR/$genome_pref"_STAR" \
     --readFilesIn $R1_f --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN 12 \
     --outFileNamePrefix  $SDIR"/RNAseq_mapping/"$out_prefix"_to_"$genome_pref
fi
done
done




ls $SDIR"/RNAseq_mapping/"*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" | wc -l
# 376
nohup mv $SDIR"/RNAseq_mapping/"*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping &

### STORE map stats HERE /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping

mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/to_$genome_pref
cp $SDIR/RNAseq_mapping/*$genome_pref"Log.final.out" /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/to_$genome_pref
cd /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/
mkdir to_$genome_pref/Ad
mkdir to_$genome_pref/HaJ
mv  to_$genome_pref/*_Ad* to_$genome_pref/Ad
mv  to_$genome_pref/*out to_$genome_pref/HaJ
module load singularity/3.7.4
singularity exec -e --bind /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq,/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/ \
/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/multiqc_latest.sif multiqc \
to_$genome_pref/Ad -o to_$genome_pref"_Ad_STAR" --interactive
singularity exec -e --bind /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq,/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/ \
/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/multiqc_latest.sif multiqc \
to_$genome_pref/HaJ -o to_$genome_pref"_HaJ_STAR" --interactive
cp to_$genome_pref"_Ad_STAR/multiqc_report.html"  to_$genome_pref"_Ad_STAR_multiqc_report.html" 
cp to_$genome_pref"_HaJ_STAR/multiqc_report.html" to_$genome_pref"_HaJ_STAR_multiqc_report.html" 
rm -r to_$genome_pref"_Ad_STAR"
rm -r to_$genome_pref"_HaJ_STAR"
tar -czvf to_$genome_pref".tar.gz" to_$genome_pref
rm -r to_$genome_pref






############################################################################################################################
#### Merge bams
## by sp / cond

module load gcc
module load samtools/1.12 
mkdir $SDIR"/RNAseq_mapping_merged"

Tps_Tdi_all_samples=(
Tdi_F_A_Ad
Tdi_F_A_Ha
Tdi_F_A_J2
Tdi_F_A_J3
Tdi_F_A_J4
Tdi_F_A_J5
Tdi_F_B_Ad
Tdi_F_B_Ha
Tdi_F_B_J2
Tdi_F_B_J3
Tdi_F_B_J4
Tdi_F_B_J5
Tdi_F_DG_Ad
Tdi_F_FB_Ad
Tdi_F_Fe_Ad
Tdi_F_Go_Ad
Tdi_F_Go_J3
Tdi_F_Go_J4
Tdi_F_Go_J5
Tdi_F_Gu_Ad
Tdi_F_Gu_Ha
Tdi_F_Gu_J2
Tdi_F_Gu_J3
Tdi_F_Gu_J4
Tdi_F_Gu_J5
Tdi_F_Lg_Ad
Tdi_F_Lg_Ha
Tdi_F_Lg_J2
Tdi_F_Lg_J3
Tdi_F_Lg_J4
Tdi_F_Lg_J5
Tdi_F_Ta_Ad
Tdi_F_WB_Ha
Tdi_M_A_Ad
Tdi_M_AG_Ad
Tdi_M_B_Ad
Tdi_M_DG_Ad
Tdi_M_FB_Ad
Tdi_M_Fe_Ad
Tdi_M_Gu_Ad
Tdi_M_Ta_Ad
Tdi_M_Te_Ad
Tps_F_A_Ad
Tps_F_A_Ha
Tps_F_A_J2
Tps_F_A_J3
Tps_F_A_J4
Tps_F_A_J5
Tps_F_A_J6
Tps_F_B_Ad
Tps_F_B_Ha
Tps_F_B_J2
Tps_F_B_J3
Tps_F_B_J4
Tps_F_B_J5
Tps_F_B_J6
Tps_F_DG_Ad
Tps_F_FB_Ad
Tps_F_Fe_Ad
Tps_F_Go_Ad
Tps_F_Go_J3
Tps_F_Go_J4
Tps_F_Go_J5
Tps_F_Go_J6
Tps_F_Gu_Ad
Tps_F_Gu_Ha
Tps_F_Gu_J2
Tps_F_Gu_J3
Tps_F_Gu_J4
Tps_F_Gu_J5
Tps_F_Gu_J6
Tps_F_Lg_Ad
Tps_F_Lg_Ha
Tps_F_Lg_J2
Tps_F_Lg_J3
Tps_F_Lg_J4
Tps_F_Lg_J5
Tps_F_Lg_J6
Tps_F_Ta_Ad
Tps_M_A_Ad
Tps_M_A_Ha
Tps_M_A_J2
Tps_M_A_J3
Tps_M_A_J4
Tps_M_A_J5
Tps_M_AG_Ad
Tps_M_B_Ad
Tps_M_B_Ha
Tps_M_B_J2
Tps_M_B_J3
Tps_M_B_J4
Tps_M_B_J5
Tps_M_DG_Ad
Tps_M_FB_Ad
Tps_M_Fe_Ad
Tps_M_Go_Ad
Tps_M_Go_J3
Tps_M_Go_J4
Tps_M_Go_J5
Tps_M_Gu_Ad
Tps_M_Gu_Ha
Tps_M_Gu_J2
Tps_M_Gu_J3
Tps_M_Gu_J4
Tps_M_Gu_J5
Tps_M_Lg_Ad
Tps_M_Lg_Ha
Tps_M_Lg_J2
Tps_M_Lg_J3
Tps_M_Lg_J4
Tps_M_Lg_J5
Tps_M_Ta_Ad
Tps_M_Te_Ad
Tps_U_Ca_Ha
Tdi_F_WB_Ad
Tps_F_WB_Ad
Tps_M_WB_Ad
)

for s in ${Tps_Tdi_all_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/RNAseq_mapping_merged/"$s"_to_"$genome_pref"_merged.bam" "/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/"$s*"_to_"$genome_pref"Aligned.sortedByCoord.out.bam"
done

for b in $SDIR"/RNAseq_mapping_merged/"*"_to_"$genome_pref"_merged.bam"; do
echo $b
out_b=`echo $b | sed 's/.bam/_sorted.bam/'`
echo $out_b
samtools sort  -@ 40 $b -o $out_b
done

#### rm merged bams
rm $SDIR"/RNAseq_mapping_merged/"*"_to_"$genome_pref"_merged.bam"
ls $SDIR"/RNAseq_mapping_merged/"*"_to_"$genome_pref"_merged_sorted.bam" -l | wc -l
## 117

#### merge for the UTR ## here need to have sep for Ad and non-Ad
module load gcc
module load samtools/1.12
ulimit -n 100000

### ALL ADULT
samtools merge -@ 40 $SDIR"/ALL_Ad_TpsTdisamples_to_"$genome_pref"_merged.bam"    $SDIR"/RNAseq_mapping_merged/"*"_Ad_to_"$genome_pref"_merged_sorted.bam" 
samtools sort  -@ 40 $SDIR"/ALL_Ad_TpsTdisamples_to_"$genome_pref"_merged.bam" -o $SDIR"/ALL_Ad_TpsTdisamples_to_"$genome_pref"_merged_sorted.bam" ### 667G 

### ALL non-ADULT 
samtools merge -@ 40 $SDIR"/ALL_J_TpsTdisamples_to_"$genome_pref"_merged.bam"   $SDIR"/RNAseq_mapping_merged/"*"_J"*"_to_"$genome_pref"_merged_sorted.bam" 
samtools merge -@ 40 $SDIR"/ALL_Ha_TpsTdisamples_to_"$genome_pref"_merged.bam"  $SDIR"/RNAseq_mapping_merged/"*"_Ha"*"_to_"$genome_pref"_merged_sorted.bam" 

samtools merge -@ 40 $SDIR"/ALL_HaJ_TpsTdisamples_to_"$genome_pref"_merged.bam" $SDIR"/ALL_Ha_TpsTdisamples_to_"$genome_pref"_merged.bam" $SDIR"/ALL_J_TpsTdisamples_to_"$genome_pref"_merged.bam" 
samtools sort  -@ 40 $SDIR"/ALL_HaJ_TpsTdisamples_to_"$genome_pref"_merged.bam" -o $SDIR"/ALL_HaJ_TpsTdisamples_to_"$genome_pref"_merged_sorted.bam" ## 532G


#######################################################################################################################################
### RNAseq runs |  5 times - take best to add UTR onto 

###### 4a 
mkdir $SDIR/Tdi_braker_RNAseq_run_4a
cd    $SDIR/Tdi_braker_RNAseq_run_4a
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4a --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4a/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi

###### 4b
mkdir $SDIR/Tdi_braker_RNAseq_run_4b
cd    $SDIR/Tdi_braker_RNAseq_run_4b
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4b --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4b/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi

###### 4c
mkdir $SDIR/Tdi_braker_RNAseq_run_4c
cd    $SDIR/Tdi_braker_RNAseq_run_4c
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4c --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4c/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi

###### 4d
mkdir $SDIR/Tdi_braker_RNAseq_run_4d
cd    $SDIR/Tdi_braker_RNAseq_run_4d
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4d --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4d/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi


###### 4e
mkdir $SDIR/Tdi_braker_RNAseq_run_4e
cd    $SDIR/Tdi_braker_RNAseq_run_4e
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4e --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4e/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi

#10263634                                   braker_rnaseq_4a_Tdi              11:04:42  COMPLETED 
#10265790                                   braker_rnaseq_4b_Tdi              09:45:19  COMPLETED 
#10265791                                   braker_rnaseq_4c_Tdi              09:52:06  COMPLETED 
#10265792                                   braker_rnaseq_4d_Tdi              09:35:31  COMPLETED 
#10265793                                   braker_rnaseq_4e_Tdi              09:59:35  COMPLETED

awk '{print $3}' $DIR/Tdi_braker_RNAseq_run_4a/braker.gff3 | grep "gene" | wc -l
#49260
awk '{print $3}' $DIR/Tdi_braker_RNAseq_run_4b/braker.gff3 | grep "gene" | wc -l
#48200
awk '{print $3}' $DIR/Tdi_braker_RNAseq_run_4c/braker.gff3 | grep "gene" | wc -l
#49327
awk '{print $3}' $DIR/Tdi_braker_RNAseq_run_4d/braker.gff3 | grep "gene" | wc -l
#49546
awk '{print $3}' $DIR/Tdi_braker_RNAseq_run_4e/braker.gff3 | grep "gene" | wc -l
#48358





#########################################################################################################################################################################################################
### UTR

########################################## bigmem 
## add the UTR ## this will delete the bam.
### use modified gushr bigmem ## 230GB mem 48 cores

cd $DIR
git clone https://github.com/Gaius-Augustus/GUSHR.git
mv GUSHR GUSHR_DJP_bigmem
cp Accessory_scripts/GUSHR_DJP_bigmem/gushr.py GUSHR_DJP_bigmem/

## make a config file to point braker to it
gushr_config_DJP_bigmem.txt 
GUSHR_PATH=/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/GUSHR_DJP_bigmem

### cp braker dicts
cp -r $DIR/Tdi_braker_RNAseq_run_4a $SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd
cp -r $DIR/Tdi_braker_RNAseq_run_4b $SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd
cp -r $DIR/Tdi_braker_RNAseq_run_4c $SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd
cp -r $DIR/Tdi_braker_RNAseq_run_4d $SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd
cp -r $DIR/Tdi_braker_RNAseq_run_4e $SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd

cp -r $DIR/Tdi_braker_RNAseq_run_4a $SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ
cp -r $DIR/Tdi_braker_RNAseq_run_4b $SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ
cp -r $DIR/Tdi_braker_RNAseq_run_4c $SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ
cp -r $DIR/Tdi_braker_RNAseq_run_4d $SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ
cp -r $DIR/Tdi_braker_RNAseq_run_4e $SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ

### cp merged bam
nohup cp $SDIR/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd &
nohup cp $SDIR/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd &
nohup cp $SDIR/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd &
nohup cp $SDIR/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd &
nohup cp $SDIR/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd &

### cp merged bam
nohup cp $SDIR/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ &
nohup cp $SDIR/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ &
nohup cp $SDIR/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ &
nohup cp $SDIR/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ &
nohup cp $SDIR/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" $SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ &


## Ad
module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd/ALL_Ad_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd/augustus.hints.gtf --skipAllTraining


## HaJ
module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ/ALL_HaJ_TpsTdisamples_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tdi --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ/augustus.hints.gtf --skipAllTraining



#10285283                           braker_rnaseq_4a_utrBMAd_Tdi            1-02:38:07  COMPLETED 
#10285291                           braker_rnaseq_4b_utrBMAd_Tdi            1-02:12:07  COMPLETED 
#10285297                           braker_rnaseq_4c_utrBMAd_Tdi            1-02:38:32  COMPLETED 
#10285300                           braker_rnaseq_4d_utrBMAd_Tdi            1-02:16:54  COMPLETED 
#10285304                           braker_rnaseq_4e_utrBMAd_Tdi            1-02:26:18  COMPLETED 
#10285323                          braker_rnaseq_4a_utrBMHaJ_Tdi            1-01:08:03  COMPLETED 
#10285354                          braker_rnaseq_4b_utrBMHaJ_Tdi            1-02:00:19  COMPLETED 
#10285358                          braker_rnaseq_4c_utrBMHaJ_Tdi            1-01:06:25  COMPLETED 
#10285364                          braker_rnaseq_4d_utrBMHaJ_Tdi            1-01:59:42  COMPLETED 
#10285367                          braker_rnaseq_4e_utrBMHaJ_Tdi            1-01:25:56  COMPLETED 

cp -r $SDIR/Tdi_braker_RNAseq_run_4a_utrBMHaJ $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4b_utrBMHaJ $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4c_utrBMHaJ $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4d_utrBMHaJ $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4e_utrBMHaJ $DIR

cp -r $SDIR/Tdi_braker_RNAseq_run_4a_utrBMAd  $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4b_utrBMAd  $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4c_utrBMAd  $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4d_utrBMAd  $DIR
cp -r $SDIR/Tdi_braker_RNAseq_run_4e_utrBMAd  $DIR





############################################################
### merge prot + RNASeq + UTR
# Tdi_braker_prot_run_3d, Tdi_braker_RNAseq_run_4d


mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/
mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/orig_merge

tar -zxf Tdi_braker_prot_run_3d.tar.gz
./TSEBRA/bin/fix_gtf_ids.py --gtf Tdi_braker_prot_run_3d/braker.gtf   --out Tdi_braker_prot_run_3d/braker_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tdi_braker_RNAseq_run_4d/braker.gtf --out Tdi_braker_RNAseq_run_4d/braker_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tdi_braker_RNAseq_run_4d_utrBMHaJ/braker_utr.gtf --out Tdi_braker_RNAseq_run_4d_utrBMHaJ/braker_utr_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tdi_braker_RNAseq_run_4d_utrBMAd/braker_utr.gtf  --out Tdi_braker_RNAseq_run_4d_utrBMAd/braker_utr_fixed.gtf

#### no-UTR version pref_braker1.cfg (RNAseq evi more weight)
./TSEBRA/bin/tsebra.py -g Tdi_braker_prot_run_3d/braker_fixed.gtf,Tdi_braker_RNAseq_run_4d/braker_fixed.gtf \
-c TSEBRA/config/pref_braker1.cfg \
-e Tdi_braker_prot_run_3d/hintsfile.gff,Tdi_braker_RNAseq_run_4d/hintsfile.gff \
-o Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4d_combined_2.gtf

#### UTR version pref_braker1.cfg (RNAseq evi more weight)
./TSEBRA/bin/tsebra.py -g Tdi_braker_prot_run_3d/braker_fixed.gtf,Tdi_braker_RNAseq_run_4d/braker_fixed.gtf,Tdi_braker_RNAseq_run_4d_utrBMHaJ/braker_utr_fixed.gtf,Tdi_braker_RNAseq_run_4d_utrBMAd/braker_utr_fixed.gtf \
-c TSEBRA/config/pref_braker1.cfg \
-e Tdi_braker_prot_run_3d/hintsfile.gff,Tdi_braker_RNAseq_run_4d/hintsfile.gff,Tdi_braker_RNAseq_run_4d_utrBMHaJ/hintsfile.gff,Tdi_braker_RNAseq_run_4d_utrBMAd/hintsfile.gff \
-o Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4dUTR_combined_2.gtf

### add to v2 
mv Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4d_combined_2.gtf    $DIR/annotations_v2/orig_merge
mv Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4dUTR_combined_2.gtf $DIR/annotations_v2/orig_merge

python3 Accessory_scripts/fix_braker_gtf.py -i $DIR/annotations_v2/orig_merge/Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4d_combined_2.gtf -d 6 -p Tdi -G 
python3 Accessory_scripts/fix_braker_gtf.py -i $DIR/annotations_v2/orig_merge/Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4dUTR_combined_2.gtf -d 6 -p Tdi -G 
python3 Accessory_scripts/UTR_GUSHR_fix.py -g $DIR/annotations_v2/orig_merge/Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4d_combined_2_FBGgi.gff -u $DIR/annotations_v2/orig_merge/Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4dUTR_combined_2_FBGgi.gff -m 1000 -G
cp $DIR/annotations_v2/orig_merge/Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4d_combined_2_FBGgiUGF1000M0.gff $DIR/annotations_v2
mv $DIR/annotations_v2/Tdi_braker_prot_run_3d_Tdi_braker_RNAseq_run_4d_combined_2_FBGgiUGF1000M0.gff $DIR/annotations_v2/Tdi_braker_prot_run_3d_RNAseq_run_4d_combined_2_FBGgiUGF1000M0.gff 




### run BUSCO on the annotation

module load gcc
module load bedtools2/2.29.2

for i in $DIR/annotations_v2/Tdi*_combined_2_FBGgiUGF1000M0.gff; do
echo $i
out_pref=`echo $i | sed 's/.gff.*//'`
echo $out_pref
awk '$3 == "gene"' $i > $out_pref"_genes.gff"
cut -f 1,4,5,9 $out_pref"_genes.gff"        | sort -k1,3 -u > $out_pref"_genes_coord.temp" ## as some rev orient genes
cut -f 4       $out_pref"_genes_coord.temp" | sed 's/.*ID=//' | sed 's/;.*//' > $out_pref"_genes_names.temp"
cut -f 1,2,3   $out_pref"_genes_coord.temp" > $out_pref"_genes_coord_2.temp" 
paste $out_pref"_genes_coord_2.temp" $out_pref"_genes_names.temp" > $out_pref"_genes.bed"
bedtools getfasta -fi genomes/$genome_pref.fasta -bed  $out_pref"_genes.bed" -name >  $out_pref"_genes.fa"
python3 Accessory_scripts/fasta_len_0.2.py $out_pref"_genes.fa"
rm $out_pref"_genes.gff"
rm $out_pref"_genes_coord.temp"
rm $out_pref"_genes_names.temp"
rm $out_pref"_genes_coord_2.temp"
rm $out_pref"_genes.bed"
done

### run BUSCO
module load singularity/3.7.4
for f in $DIR/annotations_v2/Tdi*_genes.fa; do
out_f=`echo $f | sed 's/.fa/_BUSCO/' `
echo $f
echo $out_f
singularity exec --bind $DIR /work/FAC/FBM/DEE/tschwand/default/dparker/busco_v5.3.2_cv1.sif busco -c 30 -m genome -i $f -o $out_f -l insecta_odb10 --offline --download_path $DIR/busco_downloads/ 
done

mkdir annotations_v2/BUSCO
mv    work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/* annotations_v2/BUSCO/

##### add missing busco

python3 Accessory_scripts/Add_missing_busco.py -m annotations_v2/BUSCO/Tdi_braker_*_combined_2_FBGgiUGF1000M0_genes_BUSCO/run_insecta_odb10/missing_busco_list.tsv -b BUSCO/$genome_pref/ -g annotations_v2/Tdi_*_combined_2_FBGgiUGF1000M0.gff 
#N missing BUSCO genes:
#55
#Number of complete busco found in genome busco run:
#48





##############################################################################
### Infernal
##  https://docs.rfam.org/en/latest/genome-annotation.html

## set up in default
wget eddylab.org/infernal/infernal-1.1.2.tar.gz
tar xf infernal-1.1.2.tar.gz
cd infernal-1.1.2
./configure 
make
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
./src/cmpress Rfam.cm



cd $DIR
mkdir infernal_runs
cd $DIR/infernal_runs

### workout size
for f in  $DIR/genomes/*.fasta; do
echo $f
/work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/easel/miniapps/esl-seqstat $f
done


# Because we want millions of nucleotides on both strands, we multiply 1339820256 by 2, and divide by 1,000,000 = 2679.64051
# needed for accurate eval

##/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/genomes/Tdi_LRv5a.fasta	1292580112	2585.160224



## 42 cpu 30GB RAM # Tps run time = 12:28:31
/work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/src/cmscan -Z 2585.160224 --cut_ga --rfam --nohmmonly --tblout Tdi_LRv5a-genome.tblout --fmt 2 --cpu 40 \
--clanin /work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/Rfam.clanin /work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/Rfam.cm $DIR/genomes/Tdi_LRv5a.fasta > Tdi_LRv5a-genome.cmscan


### filter overlapping output

for f in infernal_runs/T*-genome.tblout; do
echo $f
out_f=`echo $f | sed 's/-genome.tblout/-genome.deoverlapped.tblout/'`
echo $out_f
grep -v " = " $f > $out_f
done


#ncRNA derived pseudogenes pose the biggest problem for eukaryotic genome annotation using Rfam/Infernal. Many genomes contain repeat elements that are derived from a non-coding RNA gene, sometimes in huge copy number. For example, Alu repeats in human are evolutionarily related to SRP RNA, and the active B2 SINE in mouse is recently derived from a tRNA.
#In addition, specific RNA genes appear to have undergone massive pseudogene expansions in certain genomes. For example, searching the human genome using the Rfam U6 family yields over 1000 hits, all with very high score. These are not “false positives” in the sequence analysis sense, because they are closely related by sequence to the real U6 genes, but they completely overwhelm the small number (only 10s) of expected real U6 genes.
#At present we don’t have computational methods to distinguish the real genes from the pseudogenes (of course the standard protein coding gene tricks - in frame stop codons and the like - are useless).
# The sensible and precedented method for ncRNA annotation in large vertebrate genomes is to annotate the easy-to-identify RNAs, such as tRNAs and rRNAs,
# and then trust only hits with very high sequence identity (>95% over >95% of the sequence length) to an experimentally verified real gene. tRNAscan-SE has a very nice method for detecting tRNA pseudogenes.
#

### so don't trust the ncRNA (other than rRNA and tRNA)?

### need to filter output and make gff.
### filter eval.
### store full output.

for f in infernal_runs/*-genome.deoverlapped.tblout; do
prefix=`echo $f | sed 's/-genome.deoverlapped.tblout.*//' | sed 's/.*\///' | sed 's/_LR.*//'`
echo $f
echo $prefix
python3 ~/Gen_BioInf/infernal_to_gff.py -i $f -g annotations_v2/$prefix"_braker_prot_run"*_combined_2_FBGgiUGF1000M0_wB.gff ### default evalue (1e-10)
done


#Tdi: Number of infernal records kept: 465

# add to gff.




