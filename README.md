# Centromere sequence identification in Timema douglasi

This repository contains the codes used for analysing centromere sequences in the "Functional monocentricity with holocentric characteristics and chromosome-specific centromeres in a stick insect" article.

Input data files for executing the codes are available upon request.

## Documentation

script_chip_tdi_paper_final.sh describes the general pipeline to analyse the centromere sequences in T. douglasi.

The separated folders include the different scripts used in the overall pipeline:

### TE_annotation
script_transposable_element_annotation.Rmd: annotates transposable elements in the T. douglasi genome assembly.

### TR annotation and minimal rotations
script_tandem-repeat_annotation.sh: annotates tandem repeat sequences in the T. douglasi genome assembly.

script_minimal_rotation_parse.pl: orders every repeated motif sequence alphabetically.

### Levenshtein distances
script_levenshtein_rotations_F-F_inputFiles.py: computes pair-wise levenstein distances between motif sequences.

script_levenshtein_rotations_F-R_inputFiles.py: computes pair-wise levenstein distances between motif sequences and their reverse complements.

### Rscripts
TRF_parsing.R: parses gff3 file obtained from tandem repeat annotation (see script_tandem-repeat_annotation.sh)

Proportion_categories.R: estimates proportion and enrichment of sequence categories annotated in the T. douglasi genome assembly.

Enriched_windows.R: selects 10kb windows based on coverage ratio of CenH3-ChIP to input

Enriched_minimal_rotation_TR_motifs.R: extract and duplicates motif sequences with minimal rotations to compute Levenshtein distances

Levenstein_network.R: builds network of sequence similarities among tandem repeat motifs identified in the genome assembly.

Levenstein_network_contigs.R: builds network of sequence similarities among tandem repeat motifs identified in the de novo contigs (k-mer approach).

Heatmap_TRF-based.R: creates a heatmap with hierarchical clustering based on total array length inferred by summing Tandem Repeat Finder array lengths per repeat family

Heatmap_blast-based.R: creates a heatmap with hierarchical clustering based on total array length inferred by summing the lengths of sequence motif blast hits with 80% sequence similarity and 80% query coverage.
