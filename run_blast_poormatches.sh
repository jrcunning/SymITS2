# $1 is path to nw_tophits, $2 is path to representative sequences

dir=$(dirname $1)

# Get names of otus whose representative sequence is <90% similar to best hit
awk '$3 < 90 {print $2}' $1 | sort -u > $dir/poor_matches.txt
# Get corresponding sequences
awk 'NR==FNR {a[">"$1]; next} $1 in a {print; getline; print}' $dir/poor_matches.txt $2 > $dir/poor_matches_seqs.fasta

# Blast sequences to NCBI nr databse
blastn -db nr -remote -query $dir/poor_matches_seqs.fasta -outfmt '6 qseqid qcovs stitle sseqid' -max_target_seqs 1 | sort -u -k1,1 > $dir/poormatch_IDs1.txt

join -a 1 $dir/poor_matches.txt $dir/poormatch_IDs1.txt > $dir/poormatch_IDs.txt

# Clean up
rm $dir/poor_matches.txt
rm $dir/poor_matches_seqs.fasta
rm $dir/poormatch_IDs1.txt
