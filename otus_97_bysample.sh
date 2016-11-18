#!/bin/bash

# Takes input arguments: 1=path to input fasta file; 2=output directory path

# Split into individual sample files
rm -rf $2
split_sequence_file_on_sample_ids.py -i $1 -o $2

# OTU clustering at 97% similarity for each sample (in parallel)
clust97() {
	s=$(basename $1 | cut -d. -f1)
	dir=$(dirname $1)
	pick_otus.py -i $1 -s 0.97 --denovo_otu_id_prefix $s'_denovo' --optimal_uclust -o $dir
}
export -f clust97
parallel clust97 ::: $2/*.fasta

# Pick representative sequence set for each sample
repset97() {
	s=$(basename $1 | cut -d_ -f1)
	dir=$(dirname $1)
	pick_rep_set.py -i $1 -m most_abundant -f $dir/$s'.fasta' -o $dir/$s'_rep_set.fasta'
}
export -f repset97
rm -f $2/all_rep_set_otus.txt  # remove all_rep_set_otus if it already exists
parallel repset97 ::: $2/*_otus.txt

# Merge all rep sets together
rm -f $2/all_rep_set.fasta  # remove all_rep_set if it already exists
rm -f $2/all_rep_set_rep_set.fasta  # remove all_rep_set_rep_set if it already exists
cat $2/*_rep_set.fasta > $2/all_rep_set.fasta

# Cluster all 97% rep sets at 100% identity
pick_otus.py -i $2/all_rep_set.fasta -s 1.0 --optimal_uclust -o $2

# Get rep set of 100% OTUs
pick_rep_set.py -i $2/all_rep_set_otus.txt -f $2/all_rep_set.fasta -o $2/all_rep_set_rep_set.fasta

# Concatenate 97% OTU maps and merge with 100% OTU map
cat $2/*_otus.txt > $2/all_97_otus.txt
merge_otu_maps.py -i $2/all_97_otus.txt,$2/all_rep_set_otus.txt -o $2/merged_otu_map.txt

# Remove singletons
awk '$3 ~ /./ {print}' $2/merged_otu_map.txt > $2/nosingles_otus.txt

# Make OTU table
make_otu_table.py -i $2/nosingles_otus.txt -o $2/97_otus_bysample.biom
rm -f  $2/97_otus_bysample.tsv  # delete old OTU .tsv if present
biom convert -i $2/97_otus_bysample.biom -o $2/97_otus_bysample.tsv --to-tsv
