# FILTER OUT OTUS THAT ARE NOT SYMBIODINIUM
library(stringr)
library(Biostrings)
library(phyloseq)

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
# If two arguments not provided, return an error
if (length(args) < 3) {
  stop("must specify 1: phyloseq object (.RData); 2: sequences; 3: output filename", call.=FALSE)
}

# Get tax data from phyloseq object
load(args[1])
tax <- data.frame(tax_table(phy), stringsAsFactors=F)

# Identify taxa with poor matches to reference database
poortax <- subset(tax, as.numeric(sim) < 70)$otu

# Get corresponding sequences
seqs <- readDNAStringSet(args[2])
names(seqs) <- gsub(" .*$", "", names(seqs))
poorseqs <- subset(seqs, names(seqs) %in% poortax)

writeXStringSet(poorseqs, filepath="data/poorseqs.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

# Blast poor matching sequences to NCBI nr databse, get top hit
system("/Users/jrcunning/ncbi-blast-2.2.31+/bin/blastn -db nr -remote -query data/poorseqs.fasta -outfmt '6 qseqid qcovs stitle sseqid' -max_target_seqs 1 | sort -u -k1,1 > data/poorseqs_blast_results.txt")

# Read BLAST results
poorseqs_blast <- readLines("data/poorseqs_blast_results.txt")
# If the top hit from NCBI does not contain the string "Symbiodinium", then this sequence is assumed to not be Symbiodinium.
poorseqs_sym <- data.frame(otu=str_extract(poorseqs_blast, "denovo[^\t]*"),
                           symbio=str_detect(poorseqs_blast, "Symbiodinium"),   # TRUE if Symbiodinium
                           stringsAsFactors = FALSE)

symbio <- merge(data.frame(otu=poortax), poorseqs_sym, all=T)
symbio[is.na(symbio)] <- FALSE  # If NA then no BLAST hit found... assume NOT Symbiodinium

# Remove OTUs that are not Symbiodinium
phy.f <- prune_taxa(!taxa_names(phy) %in% symbio[symbio$symbio==FALSE, "otu"], phy)

# Save filtered phyloseq object
save(phy.f, file=args[3])
