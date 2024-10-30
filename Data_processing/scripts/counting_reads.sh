#!/bin/bash

for Gate in 2 3 4 5 6 7
do
	for Lib in 1 3
	do 	
		echo "Counting reads in sample ${Gate}_${Lib}"
		bam_file="Alignment/alignment_${Gate}_${Lib}_R1.bam"
		output="Count_table/counts_ampliconseq2023_${Gate}_${Lib}"	
		for seq in $(cat seq_id.txt); do
		    num_reads=$(samtools view -c -F 4 $bam_file $seq)
		    echo -e "$seq\t$num_reads"
		done > $output
	done
done
