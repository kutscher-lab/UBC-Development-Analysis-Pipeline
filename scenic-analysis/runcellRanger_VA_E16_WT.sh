#!/bin/bash

reference_path="path/to/mouse/ref/dir/refdata-gex-GRCm39-2024-A"
fastq_path="path/to/fastq/sequence"
sample="e16_tbr2_wt" #change here based on the sample name
id_seq="E16_Lmx1a-cre_TdTom_Tbr2_wt" #change here based on the ID


cd /path/to/dir
module load CellRanger/8.0.1
cellranger count --id="$id_seq" \
--transcriptome="$reference_path" \
--fastqs="$fastq_path" \
--sample="$sample" \
--create-bam=false
--chemistry=SC3Pv1
