#!/bin/bash

# use accession numbers to download sequences from NCBI in parallel
getseq() {
	n=$1
	acn=${n##*_}
	echo $1 > data/dbseqs/${acn}.fasta
	curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acn}&rettype=fasta" >> data/dbseqs/${acn}.fasta
}
export -f getseq
rm -rf data/dbseqs; mkdir data/dbseqs
parallel -a data/accn_nos.txt getseq

# Put all sequences together and reformat to keep only subtype_accn as sequence name, and put sequence all on one line
cat data/dbseqs/*.fasta | \
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | \
awk '!/^>gi/ {print}' > data/ITS2db_raw.fasta

# Clean up
rm -rf data/dbseqs
