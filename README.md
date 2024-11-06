# transcriptional_library

## Library_data
### data
- <b>promoters.fa:</b> FASTA file for the whole ordered library as received from Fiszbein lab.
- <b>spipke_in.fa:</b> FASTA for the thre Spike In used in 2023 experiment.
- <b>enhancers_hek.bed: </b> Bed file indicating the position of the enhancers to use.
- <b>promoters_wo_dupl.fa:</b> Same as promoters.fa, after removing duplicated sequences. Output of /scripts/Reports/seq_analysis.qmd and /scripts/Library_data.R
- <b>promoters_wo_dupl_SpikeIns.fa:</b> Same as promoters_wo_dupl.fa, but after the addition of the three sequences used as Spike In in 2023's experiment. Output of /scripts/Reports/seq_analysis.qmd and /scripts/SpikeIn_add_to_fasta.R
- <b>library_index:</b> Indexation of the library FASTA, for alignment by HISAT2. 

### res
- <b>seq_id.txt</b>: List of IDs in Library_data/data/promoters.fa. Output of /scripts/Reports/seq_analysis.qmd and /Library_data/bash_scripts_get_fasta_id.sh
- <b>seq_data.txt</b>: Output of /scripts/Reports/seq_analysis.qmd and /scripts/Library_data.R. Table indicating for each sequence in the original library:
      - seq_id (ID as is in the FASTA, unique for each sequence)
      - name (shared by variants of a sequence, as mutations; removes the EPD prefix)
      - gene_sym (gene associated with the ID, as indicated in the name. None for enhancers)
      - type (class of the sequence: promoter, enhancer, cancer_wt, cancer_mut)
      - dupl_id (indicating the name of the deprecated ID that is dupllicated for each pair)
  
- <b>cancer_ids.txt</b>: List of IDs in cancer WT, disregarding those from EPD database considered class=promoter.
- <b>cancer_wt.fa</b>: FASTA for the sequences in cancer_ids.txt
- <b>cancer_wt_blast.txt</b>: Output of /Library_data/bash_scripts/blastcode.sh. BLASTN for the sequences in cancer_wt.fa in genome GRCh38.p14
- <b>library.bed</b>: Bed file with the genomic position of each sequence in the library. Without duplicates nor cancer_mut sequences (its positions are the same as tthe WT pair). Output of /scripts/Reports/seq_analysis.qmd and /scripts/Library_data.R. 
-  <b>mutation_position.tsv</b>: Output of /scripts/Reports/seq_analysis.qmd and /scripts/mut_pos.R. Table indicating the relative position of the mutation of cancer_mut sequence respect to cancer_wt.
-  <b>1bp_library.bed</b>: Same as library.bed, but with 1bp ranges.
-  <b>upstream_250_library.bed</b>: BED file with only the 250pb upstream of each original sequence.
-  <b>downstream_250_library.bed</b>: BED file with only the 250pb downstream of each original sequence.

### bash_scripts
Some bash script to assist in the creation of Library_data/res files
- <b>get_fasta_id.sh</b>: Input Library_data/data/promoters.fa and outputs Library_data/res/seq_id.txt.
- <b> get_prom_symbol.sh</b>: Input Library_data/res/seq_id.txt. Creates list of name for the IDs corresponding to EPD promoters.
- <b>blastcode.sh</b>: BLASTN for Library_data/res/cancer_wt.fa. Output is Library_data/res/cancer_wt_blast.txt
- <b>phyloP_lib.sh</b>: multiBigwigSummary to convert from .BW file to .BED file and restrict to the desired regions. Assigning phyloP score to the library.bed. 


