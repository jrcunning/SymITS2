#!/bin/bash

# Split into individual sample files
rm -rf data/fasta/bysample
split_sequence_file_on_sample_ids.py -i data/fasta/combined_seqs_trimmed.fasta -o data/fasta/bysample

# OTU clustering at 97% similarity for each sample (in parallel)
clust97() {
	s=`echo $1 | grep -Eoh '/[0-9]{1,4}' | cut -c 2-`
	pick_otus.py -i $1 -s 0.97 --denovo_otu_id_prefix $s'_denovo' --optimal_uclust -o data/otus_97_bysample
}
export -f clust97
parallel clust97 ::: data/fasta/bysample/*.fasta

# Pick representative sequence set for each sample
repset97() {
	s=`echo $1 | grep -Eoh '/[0-9]{1,4}' | cut -c 2-`
	pick_rep_set.py -i $1 -m most_abundant -f 'data/fasta/bysample/'$s'.fasta' \
	-o 'data/otus_97_bysample/'$s'_rep_set.fasta'
}
export -f repset97
parallel repset97 ::: data/otus_97_bysample/[0-9]*_otus.txt

# Merge all rep sets together
cat data/otus_97_bysample/[0-9]*_rep_set.fasta > data/otus_97_bysample/all_rep_set.fasta

# Cluster all 97% rep sets at 100% identity
pick_otus.py -i data/otus_97_bysample/all_rep_set.fasta -s 1.0 --optimal_uclust -o data/otus_97_bysample

# Get rep set of 100% OTUs
pick_rep_set.py -i data/otus_97_bysample/all_rep_set_otus.txt -f data/otus_97_bysample/all_rep_set.fasta \
-o data/otus_97_bysample/all_rep_set_rep_set.fasta

# Concatenate 97% OTU maps and merge with 100% OTU map
cat data/otus_97_bysample/*_otus.txt > data/otus_97_bysample/all_97_otus.txt
merge_otu_maps.py -i data/otus_97_bysample/all_97_otus.txt,data/otus_97_bysample/all_rep_set_otus.txt -o data/otus_97_bysample/merged_otu_map.txt

# Remove singletons
awk '$3 ~ /./ {print}' data/otus_97_bysample/merged_otu_map.txt > data/otus_97_bysample/nosingles_otus.txt

# Make OTU table
make_otu_table.py -i data/otus_97_bysample/nosingles_otus.txt -o data/otus_97_bysample/97_otus_bysample.biom
rm -f  data/otus_97_bysample/97_otus_bysample.tsv  # delete old OTU .tsv if present
biom convert -i data/otus_97_bysample/97_otus_bysample.biom -o data/otus_97_bysample/97_otus_bysample.tsv --to-tsv
