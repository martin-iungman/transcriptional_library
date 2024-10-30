#!/bin/bash

fastq1="FASTQ/Gate-$1_R1.fastq.gz"
fastq2="FASTQ/Gate-$1_R2.fastq.gz"
trimmed1="Trimming/trimmed_fastq_$1_R1.fastq"
trimmed2="Trimming/trimmed_fastq_$1_R2.fastq"
trimmed_summary="Trimming/cutadapt_summary_$1.fastq"
cutadapt -g  file:pcr_primers.fa -G file:pcr_primers.fa --discard-untrimmed -o $trimmed1 -p $trimmed2 $fastq1 $fastq2 --info-file $trimmed_summary

