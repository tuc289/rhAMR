#!/bin/bash

# Modified Amr++ pipeline for nanopore sequences
# Usage: AMRplusplus_nano.sh <input_file.fastq>
# This script is in-house script designed for the custom environment & downloaded AMRplusplus V2 piepline


# Prerequisite -  export PATH=$PAHT:/gpfs/group/jzk303/default/data/tuc289/rhAMR/amrplusplus_v2/bin
# conda activate thesis (all the dependencies were installed in this environment)

input = $1
sample.id = ${1%.fastq}
amr = /gpfs/group/jzk303/default/data/tuc289/rhAMR/Nanopore_June1st/fastq_pass/combined/megares_modified_database_v2.00.fasta
annotation = /gpfs/group/jzk303/default/data/tuc289/rhAMR/Nanopore_June1st/fastq_pass/combined/megares_modified_annotations_v2.00.csv

#Alginment to MEGAres
#1. indexing
bwa index ${database_path}/megares_modified_database_v2.00.fasta

#2. Alignment
bwa mem ${amr} ${input} -t 20 > ${sample.id}.amr.alignment.sam

#Running ResistomeAnalyzer
resistome -ref_fp ${amr} -annot_fp ${annotation} -sam_fp ${sample.id}.amr.alignment.sam -gene_fp ${sample.id}.gene.tsv -group_fp ${sample.id}.group.tsv -mech_fp ${sample.id}.mechanism.tsv -class_fp ${sample.id}.class.tsv -type_fp ${sample.id}.type.tsv -t 0
