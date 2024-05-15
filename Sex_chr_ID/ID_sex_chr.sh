DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/sex_chr'
SDIR='/scratch/dparker/sex_chr'
mkdir -p $DIR
mkdir -p $SDIR


#########################################################################################################################
#### get reads

### reads stored /nas/FAC/FBM/DEE/tschwand/sex_chromosomes/D2c/raw_reads/

#2 male Tdi - loc 4 Orr spring road. Samples 18-3997, 18-3998. 
#5 female Tdi from Jaron KS, Parker DJ, Anselmetti Y, Van PT, Bast J, Dumas Z, et al. (2022). Convergent consequences of parthenogenesis on stick insect genomes. Science Advances 8: eabg3842.

########## clean and trim

mkdir -p $SDIR/temp_reads_for_cc
cd $SDIR

# rename 
for file in $DIR/READS/raw_reads/*.gz; do

	echo $file
	i=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_L._R._.*//' `
	lane_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_R._.*//' | sed 's/.*_//'`
	R_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R*//' | sed 's/_.*//'`	
	file_number=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._//' | sed 's/_.*//' `
	Group_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._..._//' | sed 's/_//g'  `
	echo $i
	echo $lane_N	
	echo $R_N
	echo $file_number
	echo $Group_N

	sp=""

	elif [ $i = "18-3997" ]; then
		   sp="Tdi_M"
	elif [ $i = "18-3998" ]; then
		   sp="Tdi_M"

	elif [ $i = "ReSeq_Di02" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di04" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di06" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di08" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di10" ]; then
		   sp="Tdi_F"
		   
		   
	else
			foo5="nope"
	fi

	new_name=`echo $sp"_"$i"_R"$R_N"_G"$Group_N$lane_N"_"$file_number".fq.gz"` 
	echo $new_name
	echo ""
	cp $file $SDIR"/temp_reads_for_cc/"$new_name
done


#### cat and clean by read group

cd $SDIR
mkdir $SDIR/cat_clean_by_read_group

for i in temp_reads_for_cc/*_R1_* ; do
	foo_R1=`echo $i`
	foo_R2=`echo $i | sed 's/_R1_/_R2_/g'`
	base_out_name=`echo $i | sed 's/.*\///' | sed 's/R1.*//' `
	read_group=`echo $i | sed 's/.*\///' | sed 's/.*R1_//' | sed 's/_.*//'`	
	
	out_R1=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R1_CC.fq"  `
	out_R2=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R2_CC.fq"  `
	echo $foo_R1
	echo $foo_R2
	echo $base_out_name
	echo $read_group
	echo $out_R1
	echo $out_R2
	echo ""
	zcat $foo_R1  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' >> $out_R1
	zcat $foo_R2  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' >> $out_R2
done


### zip

gzip $SDIR/cat_clean_by_read_group/*fq

rm -r $DIR/READS/raw_reads
rm -r $SDIR/temp_reads_for_cc 


#####################################################################################################
#####################################################################################################
### read trimming


mkdir  $DIR/READS/trimmed_reads_by_RG
cd     $DIR/READS/trimmed_reads_by_RG

cp Accessory_scripts/AllIllumina-PEadapters.fa .

module load gcc
module load trimmomatic/0.39

for i in $SDIR/cat_clean_by_read_group/*_R1_CC.fq.gz ; do
        foo1=`echo $i`
		basename=`echo $foo1 | sed 's/_R1_CC.fq.gz.*//' | sed 's/.*\///'`
        infileR1=`echo $foo1`
        infileR2=`echo $foo1 | sed 's/_R1_CC.fq.gz/_R2_CC.fq.gz/'`
        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq.gz"`
        outfileR2=`echo "./"$basename"_R2_qtrimmed.fq.gz"`
        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq.gz"`
        outfileR2_UP=`echo "./"$basename"_R2_qtrimmed_UNPAIRED.fq.gz"`
        
        echo $infileR1
        echo $infileR2
		echo $basename
        echo $outfileR1
        echo $outfileR1_UP
        echo $outfileR2
        echo $outfileR2_UP
        trimmomatic PE -threads 30 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
done


## tidy up

rm $DIR/READS/trimmed_reads_by_RG/*UNPAIRED.fq.gz
rm $SDIR/cat_clean_by_read_group/*fq.gz


#####################################################################################################
#####################################################################################################
### map reads

### prep refs

module load gcc
module load bwa/0.7.17 
module load samtools/1.15.1 
 
cd   $DIR/REFS


for f in ./*fasta; do
echo $f
bwa index $f
done


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### map reads as paired reads with BWA # gave 40GB RAM
### removes multi + supp reads



mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tdi_LRv5a_mtDNAv350"

for i in $read_dir/Tdi*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt



#########################################################################################################################
##### merge bams by samp


mkdir $SDIR/BWA_out/mapped_as_paired_merged

Tdi_samples=(
Tdi_F_ReSeq_Di02
Tdi_F_ReSeq_Di04
Tdi_F_ReSeq_Di06
Tdi_F_ReSeq_Di08
Tdi_F_ReSeq_Di10
Tdi_M_18-3997
Tdi_M_18-3998
)


genome_pref="Tdi_LRv5a_mtDNAv350"
for s in ${Tdi_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done


#########################################################################################################################
##### remove PCR duplicates  (Not just mark, as I don't think angsD pays attention to this) (actually it should, but I also need to calc cov.)
###  60GB RAM


module load gcc
module load bwa/0.7.17
module load samtools/1.15.1
module load picard/2.26.2

for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30a_sorted.bam/_mapqfilt_30aDR_sorted.bam/'`
        flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
        metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

        echo $i
        echo $outbam
        echo $metric_file
        echo $flagstat_out_bam

        picard MarkDuplicates REMOVE_DUPLICATES=true \
        INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file

        samtools flagstat $outbam > $flagstat_out_bam

done




###############################################################################################################################
### First indel realignment ##

module load gcc
module load bwa/0.7.17
module load samtools/1.15.1
module load picard/2.26.2
module load gatk/3.8.1 


## index and dict of fasta

for f in $DIR/REFS/*_LR*_mtDNAv350.fasta  ; do 
    dict_name=`echo $f | sed 's/.fasta/.dict/'`
    samtools faidx $f
	picard CreateSequenceDictionary R= $f O= $dict_name
done


# indel realignment 

ref_v="Tdi_LRv5a_mtDNAv350"
for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_to_Tdi_*_mapqfilt_30aDR_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30aDR_sorted.bam/_mapqfilt_30aDRra_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	interval_file=`echo $outbam | sed 's/.bam/.intervals/'`
	ref_fa=`echo $DIR"/REFS/"$ref_v".fasta"`
	
	echo $i
	echo $outbam
	echo $flagstat_out_bam
	echo $interval_file
	echo $ref_fa
	echo ""

	# index bam
	samtools index $i
	
	## make target intervals list
	gatk -T RealignerTargetCreator -R $ref_fa -I $i -o $interval_file
	
	## realign
	gatk -T IndelRealigner -R $ref_fa --targetIntervals $interval_file -I $i -o $outbam
	samtools flagstat $outbam > $flagstat_out_bam	

done


##################################################################################################################################################################################3
##### coverage - genome hist 


module load gcc
module load bedtools2 # bedtools v2.30.0

for s in $DIR/BWA_out/mapped_as_paired_merged_recip/*_30aDRra_sorted.bam; do
outfile=`echo $s | sed 's/_sorted.bam/_cov.txt/'`
echo $s
echo $outfile
genomeCoverageBed -ibam $s > $outfile
done


##################################################################################################################################################################################3
##### coverage - per base


module load gcc
module load bedtools2 # bedtools v2.30.0

for s in $DIR/BWA_out/mapped_as_paired_merged/*_30aDRra_sorted.bam; do
outfile=`echo $s | sed 's/_sorted.bam/_covBEDGRAPHbga.txt/'`
echo $s
echo $outfile
genomeCoverageBed -ibam $s -bga > $outfile
done



### per window
## ~30 mins to run per samp

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/dparker/mapping_to_LR_genomes/BWA_out/mapped_as_paired_merged


#### 100000 window
for s in $DIR/BWA_out/mapped_as_paired_merged/*_30aDRra_covBEDGRAPHbga.txt; do
outfile=`echo $s | sed 's/_covBEDGRAPHbga.txt//'`
echo $s
echo $outfile
python3 Accessory_scripts/bedgraph_cov_to_windows.py -i $s -w 100000 -o $outfile
done

## saved in output/

## then use plot Tdi_cov_plot_LR.R


