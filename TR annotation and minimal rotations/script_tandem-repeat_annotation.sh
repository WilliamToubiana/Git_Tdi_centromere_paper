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



###############################################################
## Tandem repeat annotation on Tdi genome assembly using TRF ##
###############################################################

module load gcc
module load trf
module load perl
module load r

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip
mkdir TR_annotation_timema
cp genomes/Tdi_LRv5a_mtDNAv350.fasta TR_annotation_timema/

#run TRF and generate gff3 output file
trf TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h

trf2gff -i TR_annotation_timema/Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)

#Parse TRF
./TRF_parsing.R #this script was used to parse TRF output .gff3 file
#Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt file was subsequently generated
#Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed file was subsequently generated



#########################################################
## Tandem repeat annotation on enriched kmers ChIPtigs ##
#########################################################
module load gcc
module load trf
module load perl

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/mer_centromere/spades_assembly

trf contigs.fasta 2 7 7 80 10 50 2000 -d -m -h

trf2gff -i contigs.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing

#Parse TRF
./TRF_parsing.R #this script was used to parse TRF output contigs.fasta.2.7.7.80.10.50.2000.gff3 file

#convert TR motifs with minimal rotation
perl script_minimal_rotation_parse.pl contigs.fasta.2.7.7.80.10.50.2000_parse.txt > contigs.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt #parse.txt file was generated from gff3 in R script
