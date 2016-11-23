#!/bin/bash

# Cluster entire dataset at 100% similarity
pick_otus.py -i $1 -s 1.0 -o $2

# Remove singletons before picking rep set and assigning taxonomy (to save time)
awk '$3 ~ /./ {print}' $2/$3 > $2/nosingles_otus.txt

# Get rep set of 100% OTUs
pick_rep_set.py -i $2/nosingles_otus.txt \
-f $1 \
-o $2/100_otus_rep_set.fasta

# Make OTU table
make_otu_table.py -i $2/nosingles_otus.txt -o $2/100_otus.biom
rm -f  $2/100_otus.tsv  # delete old OTU .tsv if present
biom convert -i $2/100_otus.biom -o $2/100_otus.tsv --to-tsv




