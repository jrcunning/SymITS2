#!/bin/bash

# Arguments: 1: input fasta file
#            2: output directory

# Delete output directory if it exists to avoid recursion
rm -rf $2

# Cluster entire dataset at 97% similarity
# Cluster entire dataset at 100% similarity
pick_otus.py -i $1 -s 0.97 -o $2

# Remove singletons before picking rep set and assigning taxonomy (to save time)
gawk '$3 ~ /./ {print}' $2/*_otus.txt > $2/nosingles_otus.txt

# Get rep set of 97% OTUs
pick_rep_set.py -i $2/nosingles_otus.txt \
-f $1 \
-o $2/97_otus_rep_set.fasta

# Make OTU table
make_otu_table.py -i $2/nosingles_otus.txt -o $2/97_otus.biom
rm -f  $2/97_otus.tsv  # delete old OTU .tsv if present
biom convert -i $2/97_otus.biom -o $2/97_otus.tsv --to-tsv
