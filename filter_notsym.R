# FILTER OUT OTUS THAT ARE NOT SYMBIODINIUM
library(stringr)
library(Biostrings)
library(phyloseq)

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
# If two arguments not provided, return an error
if (length(args) < 3) {
  stop("must specify 
       1: phyloseq object (.RData); 
       2: sequences; 
       3: output filename
       4: local nt blast database", call.=FALSE)
}

# Get tax data from phyloseq object
load(args[1])
tax <- data.frame(tax_table(phy), stringsAsFactors=F)

# Identify taxa with poor matches to reference database
poortax <- subset(tax, as.numeric(sim) < 90)$otu

# Get corresponding sequences
seqs <- readDNAStringSet(args[2])
names(seqs) <- gsub(" .*$", "", names(seqs))
poorseqs <- subset(seqs, names(seqs) %in% poortax)

poorseqsfile <- file.path(dirname(args[2]), "poorseqs.fasta")
writeXStringSet(poorseqs, filepath=poorseqsfile, append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")





# Blast poor matching sequences to local NCBI nt databse, get top hit, write to file
outfile <- file.path(dirname(args[2]), "poorseqs_blast_results.txt")
system(paste("blastn -db", args[4],
            "-query", poorseqsfile, "-num_threads 100",
            "-outfmt '6 qseqid qcovs stitle sseqid' -max_target_seqs 1 | sort -u -k1,1 >",
            outfile))

# Split poor-matching sequences into files with 100 seqs each for blasting to nt database
# poorseqs <- readDNAStringSet("~/Documents/Academia/HIMB/USVI/STJ2015/data/otus_97_bysample/poorseqs.fasta")
# poorseqsbins <- split(poorseqs, cut(1:length(poorseqs), breaks=length(poorseqs) %/% 100 + 1, labels=seq(1,length(poorseqs) %/% 100 + 1)))
# sapply(names(poorseqsbins), function(x) {
#   writeXStringSet(poorseqsbins[[x]], filepath = file.path(dirname(args[2]), paste0("poorseqs", x,".fasta")))
# })
#
# outfile <- file.path(dirname(args[2]), "poorseqs_blast_results.txt")
# system(paste("rm", outfile))
# poorseqsfiles <- list.files(path=file.path(dirname(args[2])), pattern="^poorseqs[0-9]", full.names=T)
# 
# library(parallel) # Set up for parallel processing
# cl <- makeCluster(detectCores())  # Initiate cluster
# clusterEvalQ(cl, library(Biostrings))  # Make Biostrings library available to cluster
# clusterExport(cl, c("poorseqsfiles", "outfile"))  # Export objects to cluster
# 
# parSapply(cl, poorseqsfiles, function(x) {
#   system(paste("blastn -db", args[4],
#                "-query", x,
#                "-outfmt '6 qseqid qcovs stitle sseqid' -max_target_seqs 1 | sort -u -k1,1 >>",
#                outfile))
# })
# 
# stopCluster(cl)  # Stop cluster


# Read BLAST results
poorseqs_blast <- readLines(outfile)
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
