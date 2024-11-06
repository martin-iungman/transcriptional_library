#!/bin/sh

wget -O External_data/genome/GCF_000001405.40_GRCh38.p14_genomic.fna https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED 

awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' External_data/genome/GCF_000001405.40_GRCh38.p14_genomic.fna | grep -E "^>NC" | tr "\t" "\n" > External_data/genome/GRCh38.p14_canonical.fa

wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz    | gunzip -c    | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }'   > External_data/cpgIslandExt.hg38.bed

wget -O External_data/hg38.phyloP100way.bw https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw 

wget -O External_data/Sup1_mouse_human_Young2015.txt https://genome.cshlp.org/content/suppl/2015/08/17/gr.190546.115.DC1/Supplemental_File1.txt

wget  -O External_data/Tissue_cell_ontology https://static-content.springer.com/esm/art%3A10.1038%2Fnature12787/MediaObjects/41586_2014_BFnature12787_MOESM8_ESM.zip | gzip 