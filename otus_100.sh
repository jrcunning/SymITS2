#!/bin/bash

# Cluster entire dataset at 100% similarity
pick_otus.py -i data/fasta/combined_seqs_trimmed.fasta -s 1.0 -o data/otus_100/

# Remove singletons before picking rep set and assigning taxonomy (to save time)
awk '$3 ~ /./ {print}' data/otus_100/combined_seqs_trimmed_otus.txt > data/otus_100/nosingles_otus.txt

# Get rep set of 100% OTUs
pick_rep_set.py -i data/otus_100/nosingles_otus.txt \
-f data/fasta/combined_seqs_trimmed.fasta \
-o data/otus_100/100_otus_rep_set.fasta

# Make OTU table
make_otu_table.py -i data/otus_100/nosingles_otus.txt -o data/otus_100/100_otus.biom
rm -f  data/otus_100/100_otus.tsv  # delete old OTU .tsv if present
biom convert -i data/otus_100/100_otus.biom -o data/otus_100/100_otus.tsv --to-tsv




