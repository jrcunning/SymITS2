#!/bin/bash

rm -rf data/otus_97_byspecies
# Split into individual sample files
rm -rf data/fasta/bysample
split_sequence_file_on_sample_ids.py -i data/fasta/combined_seqs_trimmed.fasta -o data/fasta/bysample

rm -rf data/fasta/byspecies; mkdir data/fasta/byspecies
# Get SampleIDs for each species
cd data/fasta/bysample
awk -F"\t" '{print $1".fasta" > $6"_IDs.txt"}' ../../mapping_file.txt
rm -f Species_IDs.txt
for spp in `ls *_IDs.txt`; do
	cat $(cat $spp) > ../byspecies/${spp%%_*}.fasta
done
cd ../../..

# OTU clustering at 97% similarity for each species (in parallel)
clust97() {
	s=`basename $1`
	pick_otus.py -i $1 -s 0.97 --denovo_otu_id_prefix $s'_denovo' --optimal_uclust -o data/otus_97_byspecies
}
export -f clust97
parallel -j3 clust97 ::: data/fasta/byspecies/*.fasta

# Pick representative sequence set for each sample
repset97() {
	filename=`basename $1`
	s=${filename%_*}
	pick_rep_set.py -i $1 -m most_abundant -f 'data/fasta/byspecies/'$s'.fasta' \
	-o 'data/otus_97_byspecies/'$s'_rep_set.fasta'
}
export -f repset97
parallel repset97 ::: data/otus_97_byspecies/*_otus.txt

cat data/otus_97_byspecies/*otus.txt > data/otus_97_byspecies/all_97_otus.txt
# Merge all rep sets together
cat data/otus_97_byspecies/*_rep_set.fasta > data/otus_97_byspecies/all_rep_set.fasta

# Cluster all 97% rep sets at 100% identity
pick_otus.py -i data/otus_97_byspecies/all_rep_set.fasta -s 1.0 --optimal_uclust -o data/otus_97_byspecies

# Get rep set of 100% OTUs and assign taxonomy
pick_rep_set.py -i data/otus_97_byspecies/all_rep_set_otus.txt -f data/otus_97_byspecies/all_rep_set.fasta \
-o data/otus_97_byspecies/all_rep_set_rep_set.fasta

# Concatenate 97% OTU maps and merge with 100% OTU map
merge_otu_maps.py -i data/otus_97_byspecies/all_97_otus.txt,data/otus_97_byspecies/all_rep_set_otus.txt -o data/otus_97_byspecies/merged_otu_map.txt

# Remove singletons
awk '$3 ~ /./ {print}' data/otus_97_byspecies/merged_otu_map.txt > data/otus_97_byspecies/nosingles_otus.txt

# Make OTU table
make_otu_table.py -i data/otus_97_byspecies/nosingles_otus.txt -o data/otus_97_byspecies/97_otus_byspecies.biom
rm -f  data/otus_97_byspecies/97_otus_byspecies.tsv  # delete old OTU .tsv if present
biom convert -i data/otus_97_byspecies/97_otus_byspecies.biom -o data/otus_97_byspecies/97_otus_byspecies.tsv --to-tsv
