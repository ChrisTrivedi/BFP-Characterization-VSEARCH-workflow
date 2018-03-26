#!/bin/bash

### Description ###
# This script will use a combination of QIIME and VSEARCH in order to create an OTU table that is suitable for input into multiple analysis pipelines (QIIME, R, etc.)

# This script assumes that you are in the directory that you are analyzing your data in. Paired end reads are located in Raw/ and subsequent files will be in seq/.

# You'll notice I used 'time' in front of many of the scripts. This is because I want to note how long some of them take to run depending on the data I am running.

# Also note that I use macqiime which is just the Mac version of QIIME. They are interchangable whether you are running on Mac, Windows VirtualBox, or a Linux distro.

### Dependent Software ###
# QIIME 1.9.1 at the time of this analysis
# VSEARCH (v2.3.4.osx_x86_64 used for this analysis), https://github.com/torognes/vsearch

### Workflow ###

# Concatenate raw fastq files
# Forward read (R1)
cat Raw/<raw R1.fastq.gz file> Raw/<raw R1.fastq.gz file> > Raw/<concatenated_R1.fastq.gz file>

# Reverse read (R2)
cat Raw/<raw R2.fastq.gz file> Raw/<raw R2.fastq.gz file> > Raw/<concatenated_R2.fastq.gz file>

# Validate mapping file
# To test for duplicate barcodes. If this is the case (and it may very well be if  samples were sequenced on different Illumina runs), then concatenate files (seqs.fastq) after running 'Split Libraries' and not before.
validate_mapping_file.py -o Mapping_Validate -m Mapping_file.txt

# Starting with QIIME to join paired ends, extract barcodes, and split libraries.

# Join paired ends - overlap of 50
time join_paired_ends.py -f Raw/R1.fastq -r Raw/R2.fastq -o seq/Joined/ -j 50 -v

# Extract barcodes using the mapping file - 10 mins
time extract_barcodes.py -f seq/Joined/fastqjoin.join.fastq -m Mapping_file.txt -l 12 -o seq/Prepped -a -v

# Split libraries (separate out samples from fastq) - 14 mins
time split_libraries_fastq.py --store_demultiplexed_fastq --phred_quality_threshold 20 -i seq/Prepped/reads.fastq -b seq/Prepped/barcodes.fastq -m Mapping_file.txt --barcode_type 12 -o seq/SlOut/ -v

mkdir seq/VsearchOut

# Beginning of VSEARCH portion of the workflow
# get quality stats
vsearch -fastq_stats seq/SlOut/seqs.fastq --log seq/VsearchOut/seqs_stats.log

# Remove low quality reads
time vsearch -fastx_filter seq/SlOut/seqs.fna -fastaout seq/VsearchOut/seqs_filtered.fasta --fastq_maxee 0.5 --threads 4

# Dereplicate seqs
time vsearch -derep_fulllength seq/VsearchOut/seqs_filtered.fasta --output seq/VsearchOut/seqs_filtered_derep.fasta --sizeout --minuniquesize 2 --threads 4

# Reference chimera check
time vsearch -uchime_ref seq/VsearchOut/seqs_filtered_derep.fasta --db /macqiime/DB/gold.fasta --strand plus --nonchimeras seq/VsearchOut/seqs_filtered_derep_nochimeras.fasta --threads 4

# Cluster OTUs @ 97%
time vsearch -cluster_fast seq/VsearchOut/seqs_filtered_derep_nochimeras.fasta --centroids seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta --sizein --xsize --relabel OTU_ --id 0.97 --threads 4

# Make an otus folder
mkdir otus/

# Copy this file to use as RepSet at a later time
cp seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta otus/RepSet.fna

# Map the original quality filtered reads back to OTUs
time vsearch -usearch_global seq/VsearchOut/seqs_filtered.fasta --db seq/VsearchOut/seqs_filtered_derep_nochimeras_repset.fasta --strand plus --id 0.97 -uc  seq/VsearchOut/OTU_map.uc --threads 4

# Modify OTU table for input into QIIME
python /macqiime/bin/uc2otutab.py seq/VsearchOut/OTU_map.uc > seq/VsearchOut/seqs_filtered_derep_nochimeras_repset_OTU-table.txt

# Convert to HDF5 biom type
biom convert --table-type="OTU table" -i seq/VsearchOut/seqs_filtered_derep_nochimeras_repset_OTU-table.txt --to-hdf5 -o otus/<output BIOM filename>.biom

# Summarize BIOM table to check general stats
biom summarize-table -i otus/<output BIOM filename>.biom -o otus/Summary.txt

# Moving back into QIIME
# Assign taxonomy
# You'll need to modify the "-t" flag to point to wherever your reference database is.  
# In this case I am using SILVA123 as my database of choice, but GreenGenes is the default for QIIME.
time assign_taxonomy.py -t /macqiime/SILVA/taxonomy/taxonomy_all/97/taxonomy_7_levels.txt -r /macqiime/SILVA/rep_set/rep_set_all/97/97_otus.fasta -i otus/RepSet.fna -o otus/TaxonomyOut/ -v

# Add taxonomy to BIOM table - instant
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp otus/TaxonomyOut/RepSet_tax_assignments.txt -i otus/<output BIOM filename>.biom -o otus/otuTable.biom

# Computing Summaries of Taxa
biom summarize-table -i otus/otuTable.biom -o otus/otuTable_Summary.txt

# Align seqs (default: pynast - can also do in parallel)
time align_seqs.py -i otus/RepSet.fna -t /macqiime/SILVA/core_alignment/core_alignment_SILVA123.fasta -o otus/RepSet_Aligned/ -v

# Filter alignment to certain quality (0.01 in this case)
time filter_alignment.py -i otus/RepSet_Aligned/RepSet_aligned.fasta -o otus/RepSet_Aligned/ -e 0.01

# Make tree file (.tre) for downstream use
time make_phylogeny.py -i otus/RepSet_Aligned/RepSet_aligned_pfiltered.fasta -o otus/RepSet_Aligned/rep_set.tre -l otus/RepSet_Aligned/tree_log.txt

## At this point you will have an OTU table and tree file, from here there are a number of different things to do. You can continue on in QIIME or bring your data over to R. I will include a few more commands below.

# Sort OTU table by sample type (this is optional, but can definitely make things easier for downstream QIIME analyses)
sort_otu_table.py -i otus/otuTable.biom -o otus/otuTable_sorted.biom -m <sorted mapping file>.txt -s SortOrder

# Summarize Taxa/plots on this new OTU table and then you can filter for further diversity analyses (i.e. heatmap, BetaDiversity, etc.)
time summarize_taxa_through_plots.py -i otus/otuTable_sorted.biom -o TaxaSummary -v

# Filter OTU table and convert for analysis in R
filter_samples_from_otu_table.py -i otus/otuTable_sorted.biom -o otus/otuTable_sorted_filtered.biom --sample_id_fp IDs.txt

# Check that filtering worked
biom summarize-table -i otus/sorted_filtered.biom -o otus/otuTable_filtered_Summary.txt






