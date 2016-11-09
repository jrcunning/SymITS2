#!/bin/bash

# Generate config files for merging reads with illumina-utils
mkdir -p data/merge
cp data/fastq_list.txt data/fastq/fastq_list.txt
iu-gen-configs data/fastq/fastq_list.txt -o data/merge

# Merge paired reads for each sample (in parallel)
parallel iu-merge-pairs {} --min-overlap-size 150 --enforce-Q30-check --marker-gene-stringent ::: data/merge/*.ini

# Filter sequences - keep only those with 3 mismatches or less (in parallel)
parallel iu-filter-merged-reads {} --max-mismatches 3 ::: data/merge/*_MERGED


