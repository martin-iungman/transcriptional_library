#!/bin/bash

trimmed1="Trimming/trimmed_fastq_$1_R1.fastq"
trimmed2="Trimming/trimmed_fastq_$1_R2.fastq"
sam_file="Alignment/alignment_$1.sam"
sam_summary_file="Alignment/alignment_summary_$1.sam"
bam_file="Alignment/alignment_$1.bam"
hisat2 --no-spliced-alignment -I 230 -X 270 --no-discordant -x /home/martin/Documents/Labo/splicing_noise/Library_data/data/library_index/library_idx -1 $trimmed1 -2 $trimmed2 -S $sam_file --summary-file $sam_summary_file
tar -czf $trimmed1
tar -czf $trimmed2
samtools view -b -F 256 $sam_file > $bam_file
tar -czf $sam_file
