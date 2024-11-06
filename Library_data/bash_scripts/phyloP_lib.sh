#!/bin/sh

multiBigwigSummary BED-file -b External_data/hg38.phyloP100way.bw -o hg38.phyloP100way_summary_lib.npz --BED Labo/splicing_noise/Library_data/res/library.bed --outRawCounts ~/Documents/Labo/splicing_noise/External_data/hg38.phyloP100way_summary_lib.bed

multiBigwigSummary BED-file -b External_data/hg38.phyloP100way.bw -o hg38.phyloP100way_summary_lib_1bp.npz --BED Labo/splicing_noise/Library_data/res/1bp_library.bed --outRawCounts ~/Documents/Labo/splicing_noise/External_data/hg38.phyloP100way_summary_lib_1bp.bed
