#!/bin/bash

dir=$(dirname $1)

# BLAST each sequence to the reference ITS2 database (trimmed, dereplicated)
formatdb -p F -i $2
blastall -p blastn -i $1 -d $2 -e 0.00001 -b 10 -m 8 > $dir/blast_results.txt
# Sort by query sequence and then E-value in descending order, then get only the first hit for each query sequence
# sort -k1,1 -k11,11g $dir/blast_results.txt | sort -u -k1,1 > $dir/blast_tophits.txt

# Identify the lowest E-value obtained for each sequence
sort -k1,1 -k11,11g $dir/blast_results.txt | sort -u -k1,1 | cut -f1,11 > $dir/lowevals.txt
# Get all BLAST hits that obtained the lowest E-value for each sequence (keeps ties/multiple hits with same E-val)
rm -f $dir/blast_tophits_withties.txt
while read n; do
  otu=`echo $n | cut -d" " -f1`
  eval=`echo $n | cut -d" " -f2`
  awk '$1 ~ /'"$otu"'$/ && $11 ~ /'"$eval"'$/ {print}' $dir/blast_results.txt >> $dir/blast_tophits_withties.txt
done < $dir/lowevals.txt

# Merge hits with equivalent results...(if score/length/eval/mismatch/gaps are identical for multiple hits)
# by appending multiple hits to the end of the line
awk '
BEGIN { FS="\t" }
{
    curr = $1 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9 FS $10 FS $11 FS $12
    if (curr == prev) {
        rec = rec "\t" $2
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }
' $dir/blast_tophits_withties.txt > $dir/blast_tophits_merged.txt

# Move column two (first blast hit) to last column position so all hits are together, and collapse multiple hits into a single field separated by ";"
awk 'BEGIN{FS=OFS="\t"} {a=$2; for (i=2;i<NF; i++) $i=$(i+1); $NF=a}1' $dir/blast_tophits_merged.txt | \
awk -v ORS= '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11; for(i=12;i<=NF;i++) print (i%1?OFS:";") $i; print "\n"}' | \
sed -e 's/;/ /' > $dir/blast_tophits.tsv

# Clean up 
rm -f $2.*
rm $dir/lowevals.txt
rm $dir/blast_tophits_withties.txt
rm $dir/blast_tophits_merged.txt