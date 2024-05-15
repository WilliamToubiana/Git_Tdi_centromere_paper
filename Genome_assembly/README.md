## Pipeline description for genome assembly.
### by Patrick Tran Van 

## Table of contents

1. [Oxford Nanopore reads filtration](#1)
2. [Contig assembly](#2)
3. [Contamination removal](#3)
4. [Haplotype purging](#4)

With *Timema douglasi* (Tdi).

#### <a name="1"></a>1) Oxford Nanopore reads filtration. 

Using [Filtlong](https://github.com/rrwick/Filtlong).

```
module load filtlong/0.2.0 

filtlong --min_length 1000 --keep_percent 90 --target_bases 69050000000 Tdi_nanopore_raw_reads.fastq.gz | gzip > Tdi_nanopore_filtlong50x.fastq.gz
```

#### <a name="2"></a>2) Contig assembly

2.1) Nanopore reads assembly

Using [Flye](https://github.com/fenderglass/Flye).

```
source flye/2.8.1.sh

flye --nano-raw  Tdi_nanopore_filtlong50x.fastq --out-dir Tdi_flye --genome-size 1.381g
```

2.2) First step of polishing using filtered Oxford Nanopore reads and Racon. 

Using [Minimap2](https://github.com/lh3/minimap2) and [Racon](https://github.com/lbcb-sci/racon).

```
# Map the reads against the contig assembly generated previously

module load minimap2/2.19

minimap2 -c -x map-ont Tdi_flye_contig.fasta Tdi_nanopore_filtlong50x.fastq > Tdi_flye_contig_filtlong50x.paf

# Use Racon for polishing

source racon/1.4.3.sh

racon Tdi_nanopore_filtlong50x.fastq Tdi_flye_contig_filtlong50x.paf Tdi_flye_contig.fasta > Tdi_flye_racon.fasta
```

2.3) Second step of polishing using filtered Illumina short reads and three rounds of Pilon. 

Using [BWA](https://github.com/lh3/bwa), [Samtools](https://www.htslib.org/) and [Pilon](https://github.com/broadinstitute/pilon).

```
module load bwa/0.7.17
module load samtools/1.12

# Round 1

bwa mem Tdi_flye_racon Tdi_illumina_R1.fastq.gz Tdi_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tdi_flye_racon_pilon1.sorted.bam

java -jar pilon-1.23.jar --genome Tdi_flye_racon --frags Tdi_flye_racon_pilon1.sorted.bam --output Tdi_flye_racon_pilon1 --diploid

# Round 2

bwa mem Tdi_flye_racon_pilon1 Tdi_illumina_R1.fastq.gz Tdi_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tdi_flye_racon_pilon2.sorted.bam

java -jar pilon-1.23.jar --genome Tdi_flye_racon_pilon1 --frags Tdi_flye_racon_pilon2.sorted.bam --output Tdi_flye_racon_pilon2 --diploid

# Round 3

bwa mem Tdi_flye_racon_pilon2 Tdi_illumina_R1.fastq.gz Tdi_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tdi_flye_racon_pilon3.sorted.bam

java -jar pilon-1.23.jar --genome Tdi_flye_racon_pilon2 --frags Tdi_flye_racon_pilon3.sorted.bam --output Tdi_flye_racon_pilon3 --diploid
```


#### <a name="3"></a>3) Contamination removal

Using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [BlobTools](https://github.com/DRL/blobtools) and [BBMap](https://sourceforge.net/projects/bbmap/).

3.1) Get the coverage of each contigs from the Flye output in `contigs_stats.txt`.

3.2) Blastn against nt:

```
module add Blast/ncbi-blast/2.7.1+
 
blastn -query Tdi_flye_racon_pilon3.fasta -db /nt \
-outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' \
-max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -out Tdi_flye_racon_pilon3.vs.nt.max10.1e25.blastn.out
```

3.3) Run BlobTools:

```
source blobtools/1.1.1.sh

blobtools create -i Tdi_flye_racon_pilon3.fasta -t Tdi_flye_racon_pilon3.vs.nt.max10.1e25.blastn.out \
--nodes data/nodes.dmp --names data/names.dmp \
-c coverage.txt -x bestsumorder -o Tdi_flye_racon_pilon3

blobtools blobplot -i Tdi_flye_racon_pilon3.blobDB.json --sort count --hist count -x bestsumorder

blobtools view -i Tdi_flye_racon_pilon3.blobDB.json --hits --rank all -x bestsumorder
```

3.4) Filter out contigs without hits to metazoans using our custom [script](https://github.com/AsexGenomeEvol/HD_Oppiella/blob/master/assembly/contamination_filtration.py).

```
python contamination_filtration.py -s contamination_identification -i1 Tdi_flye_racon_pilon3.blobDB.table.txt

module add UHTS/Analysis/BBMap/37.82

filterbyname.sh in=Tdi_flye_racon_pilon3.fasta names=contaminant.txt out=Tdi_flye_racon_pilon3_blobtools.fasta include=f 
```

#### <a name="4"></a>4) Haplotype purging

Using [Minimap2](https://github.com/lh3/minimap2), [Samtools](https://www.htslib.org/) and [Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/).

4.1) Map the filtered Oxford Nanopore reads against the decontaminated genome

```
module load minimap2/2.19

minimap2 -ax map-ont Tdi_flye_racon_pilon3_blobtools.fasta Tdi_nanopore_filtlong50x.fastq --secondary=no | samtools sort -o Tdi_map.sorted.bam
```

4.2) Run Purge Haplotigs

```
source purge_haplotigs/1.1.1.sh

purge_haplotigs cov -i Tdi_map.sorted.bam.gencov -l3 -m 17 -h 190 -j 101

purge_haplotigs purge -g Tdi_flye_racon_pilon3_blobtools.fasta -c coverage_stats.csv -o Tdi_flye_racon_pilon3_blobtools_purgehap
```

#### <a name="5"></a>5) Hi-C scaffolding

Using [BWA](https://github.com/lh3/minimap2), [Juicer](https://www.htslib.org/) and [3D-DNA](https://github.com/aidenlab/3d-dna).

5.1) Map Hi-C reads

```
# Create a BWA reference

module add UHTS/Aligner/bwa/0.7.17;

bwa index Tdi_flye_racon_pilon3_blobtools_purgehap.fasta 

# Generate a restriction sites file 

python juicer/generate_site_positions.py Sau3AI draft Tdi_flye_racon_pilon3_blobtools_purgehap.fasta

# Create chrom.sizes file

awk 'BEGIN{OFS="\t"}{print $1, $NF}' draft_Sau3AI.txt > Tdi_flye_racon_pilon3_blobtools_purgehap.chrom.sizes

# Run Juicer

source juicer/1.6.sh

./scripts/juicer.sh -g Tdi -z Tdi_flye_racon_pilon3_blobtools_purgehap.fasta -y draft_Sau3AI.txt -p Tdi_flye_racon_pilon3_blobtools_purgehap.chrom.sizes -D Tdi_juicer
```

5.2) Perform scaffolding

```
source 3d-dna/v180922.sh

3d-dna/run-asm-pipeline.sh --editor-coarse-resolution 25000 Tdi_flye_racon_pilon3_blobtools_purgehap.fasta merged_nodups.txt
```

Final version: `Tdi_LRv5a.fasta`

## <a name="6"></a>6) Evaluation


Using [BUSCO](http://busco.ezlab.org/).

```
source busco/5.1.2.sh

busco -m geno -i Tdi_LRv5a.fasta -o Tdi_LRv5a  -l insecta_odb10 --long
```


SM Table 1:  Parameters for Purge Haplotigs and 3D-DNA

|   **Species**         | Purge Haplotigs          | 3D-DNA                            |
|-----------------------|--------------------------|-----------------------------------|
|   *T. douglasi*       | -l 3 -m 17 -h 190 -j 101 | --editor-coarse-resolution 25000  |


