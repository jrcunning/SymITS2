# Add QIIME labels
add_qiime_labels.py -m data/mapping_file.txt -i data/merge/ -c InputFileName -o data/fasta

# Chimera checking
identify_chimeric_seqs.py -i data/fasta/combined_seqs.fna --suppress_usearch61_ref -m usearch61 -o data/fasta/usearch61_chimeras
filter_fasta.py -f data/fasta/combined_seqs.fna -o data/fasta/combined_seqs_chimera_filtered.fasta -s data/fasta/usearch61_chimeras/chimeras.txt -n

# Trim primer sequences using cutadapt
# Trim forward primers (allow error rate of 15% (3 indels/mismatches)) and keep only trimmed sequences
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 data/fasta/combined_seqs_chimera_filtered.fasta --trimmed-only -o data/fasta/trimF.fasta
# Trim forward primers again (may be multiple internal primer sequences), but do not discard sequences that do not contain primer
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 data/fasta/trimF.fasta -o data/fasta/trimF2.fasta
# Trim reverse primers using cutadapt and keep only trimmed sequences
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 data/fasta/trimF2.fasta --trimmed-only -o data/fasta/trimF2_trimR.fasta
# Trim reverse primers again for any remaining internally and remove sequences shorter than 250 bp
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 data/fasta/trimF2_trimR.fasta -m 250 -o data/fasta/combined_seqs_trimmed.fasta

# Remove intermediate files
rm data/fasta/trimF.fasta
rm data/fasta/trimF2.fasta
rm data/fasta/trimF2_trimR.fasta