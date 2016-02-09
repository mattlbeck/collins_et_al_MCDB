#!/bin/bash
set -e # exit after an error from any command
# created: 15/01/2016
# This script should run the preprocessing analysis provided correct starting files and dependencies
#  - Bioperl
#  - PatMan
# Perl modules:
#  - Sort::Naturally
# R modules:
#  - data.table
# tRNAscan v1.3.1

# Input files:
# FASTQ files available from GEO accession ###. Place all fastq files in the fastq directory
# Genome file available from ftp://ftp.hgsc.bcm.edu/Bterrestris/Bter_1.1/LinearScaffolds/ftp://ftp.hgsc.bcm.edu/Bterrestris/Bter_1.1/LinearScaffolds/Bter20110317-genome.fa
# Enter genome file here:
GENOMEFILE="./input_files/bter.fa"
# An Rfam 11 fasta file is also required to map misc. ncRNAs:
RFAM11_FASTA="./input_files/rfam11.fa"

GENOMEFORMATTED="genome.fa"
COMBINEDFASTA="preprocess/combined.fa"
PATMANOUT="mapped/mapped.pat"
BEDFILE="mapped/mapped.bed"
FASTAMAPPED="mapped/mapped.fasta"
COUNTCSV="preprocess/counts.csv"
TRNASCAN_OUT="trnas.out"
TRNASCAN_FASTA="trnas.fasta"
TRNAS_CSV="trnas.csv"
TRNAS_MAPPED="trna_mapped.pat"
TRNAS_BED="trna_mapped.bed"
RFAM_MAPPED="rfam_mapped.pat"
RFAM_CSV="rfam_mapped.csv"

mkdir preprocess || true
mkdir mapped || true
echo "Preprocessing FASTQ -> FASTA"
perl ./scripts/preprocess2.pl -c preprocess.config /scratch/mbeckers/bees/fastq/*fastq 2> preprocess/preprocess.log


# combine fasta files
perl ./scripts/combine_fasta.pl preprocess/*na.fasta > $COMBINEDFASTA

# melt reads into a csv count matrix from the separate fasta files
# this step is dpendent on the ordering of the files in the directory. If the files
# are named using the SRA/SRR id or similar, this shouldn't be a problem
perl ./scripts/meltReads.pl -r "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16" -o $COUNTCSV preprocess/*na.fasta

echo "Aligning reads to genome"
# The genome file used needs to be formatted so that lines are 70 chars long
perl ./scripts/format_genome_fasta.pl $GENOMEFILE > $GENOMEFORMATTED

# Patman was used map sequences to the genome with no mismatches
patman -P $COMBINEDFASTA -D $GENOMEFORMATTED -e 0 -o $PATMANOUT

# convert to bed format
perl ./scripts/pat2bed.pl $PATMANOUT > $BEDFILE

# retrieve mapped sequences in fasta format
perl ./scripts/pat2fasta.pl $PATMANOUT > $FASTAMAPPED

echo "Aligning reads to Rfam"
### ANNOTATION WITH Rfam and tRNAs
patman -P $FASTAMAPPED -D $RFAM11_FASTA  -e 1 -o $RFAM_MAPPED

# format Rfam identifiers from mapped file
perl ./scripts/parse_rfam_mapped.pl $RFAM_MAPPED > $RFAM_CSV

echo "Aligning reads to tRNAs"
# Using tRNAscan-SE 1.3.1
tRNAscan-SE -o $TRNASCAN_OUT $GENOME

# process tRNAs to retrieve post transcriptionally modified transcripts
perl process_trnascan.pl -t $TRNASCAN_OUT -g $GENOMEFORMATTED -f $TRNASCAN_FASTA > $TRNAS_CSV

# map to tRNAs
patman -P $COMBINEDFASTA -D $TRNASCAN_FASTA -e 0 -o $TRNAS_MAPPED
perl ./scripts/pat2bed.pl $TRNAS_MAPPED > $TRNAS_BED

# Quality Checks in R
cd ./r-work
Rscript do_qc.R $BEDFILE $COUNTCSV

# Merge together mapmi, miRCat, and miRDeep results using ~/scripts/miRCat/melt_mircat_csvs.pl
# - name miRCat predictions that overlap mapmi results
# - output only distinct mapmi results, selecting the best result each time.
# - Check miRDeep results against these
# Finish of with some stats on shared predictions and general miRNA characteristics etc. (stuff where no counts are required)
#
# # Input files required:
cd ../miRNAs
MAPMI="./data/mapmi_formatted.csv mapmi results" # formatted mapmi file
MIRCAT="./data/mircat/bees-combined-scaffmapped_s*" # all mircat results
MIRDEEP= "./data/mirdeep_results.csv" # merged mirdeep results
COUNTS="../r-work/all_counts.csv" # produced by do_qc.R
MATUREMIRS_PAT="mature_miRNAs.pat"
MATUREMIR_REGIONS="mature_miRNA_Regions.fasta"
 
perl melt_mircat_csvs2.pl \
    -m $MAPMI -M mapmi2_unique.csv -b bad_hps.csv \
    -D $MIRDEEP \
    -r "((?:s\d+)|(?:[QW]\.(?:InA\.)?(?:OA\.)?[ob][vr]\.\d))" $MIRCAT  \
    > merge.log

perl coordinates2fasta.pl -g $GENOMEFORMATTED summary_predictions.csv > summary_precursors.fasta

patman -P /scratch/mbeckers/mirbase/mature.fa -D summary_precursors.fasta -e 2 > mirbase_to_precursors.pat

perl find_mature_prediction_regions.pl -p predictions.obj -c $COUNTS \ 
    -g $GENOMEFORMATTED -f mature_miRNA_regions.fasta > mature_miRNA_regions.csv 2> mature_miRNA_regions.log

patman -P ../$COMBINEDFASTA -D $MATUREMIR_REGIONS -e 0 > $MATUREMIRS_PAT

# Annotation in R
cd ../r-work
Rscript do_annotation.R ../miRNAs/$MATUREMIRS_PAT ../miRNAs/$MATUREMIR_REGIONS ../$TRNAS_BED ../$RFAM_CSV

# Differential Expression in R
Rscript do_de.R
