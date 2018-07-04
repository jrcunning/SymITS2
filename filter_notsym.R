# FILTER OUT OTUS THAT ARE NOT SYMBIODINIUM
library(stringr)
library(Biostrings)
library(phyloseq)
library(parallel) # Set up for parallel processing
library(tidyverse)

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
# If four arguments not provided, return an error
if (length(args) < 4) {
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

# Split poor-matching sequences into X equal-sized files where X is number of cores available
ncores <- detectCores()
nseqs.core <- ceiling(length(poorseqs) / ncores)
poorseqs.split <- split(poorseqs, f = factor(ceiling(seq_along(poorseqs)/nseqs.core)))
for (i in 1:length(poorseqs.split)) {
  writeXStringSet(poorseqs.split[[i]], 
                  filepath = file.path(dirname(args[2]), paste0("poorseqs", i, ".fasta")),
                  append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
}


# TODO: ADD MINIMUM EVAL OR OTHER FILTER ON BLAST SEARCH?
# Blast poor matching sequences to local NCBI nt databse, get top hit, write to file
# outfile <- file.path(dirname(args[2]), "poorseqs_blast_results.txt")
# system(paste("blastn -db", args[4],
#             "-query", poorseqsfile, "-num_threads 100",
#             "-outfmt '6 qseqid qcovs stitle sseqid' -max_target_seqs 1 | sort -u -k1,1 >",
#             outfile))

outfile <- file.path(dirname(args[2]), "poorseqs_blast_results.txt")
system(paste("rm", outfile))
poorseqsfiles <- list.files(path=file.path(dirname(args[2])), pattern="^poorseqs[0-9]", full.names=T)


cl <- makeCluster(detectCores())  # Initiate cluster
clusterEvalQ(cl, library(Biostrings))  # Make Biostrings library available to cluster
clusterExport(cl, c("poorseqsfiles", "outfile"))  # Export objects to cluster

# create commands function
cmdCreate <- function(infile, outfile){
  paste("blastn -db ", args[4], 
        " -query ", infile, 
        " -out ", outfile, 
        " -outfmt '6 qseqid qcovs stitle sseqid'",
        " -num_threads 2",
        " -max_target_seqs 1",
        sep = "")
}

# create actual commands
cmds <- mapply(FUN = cmdCreate, infile = poorseqsfiles, outfile = paste0(poorseqsfiles, "_blastout"))

parSapply(cl, cmds, system)

stopCluster(cl)  # Stop cluster

# Additionally, closest match subject sequences were screened for two sequences in particular (JN406302 and JN406301), which are mis-annotated as Symbiodinium (highly divergent from any other Symbiodinium sequences) in the ‘nt’ database. Thus, query sequences matching these sequences were annotated as ‘Unclutured eukaryote.’ 

# # Read BLAST results
blastoutfiles <- list.files(path = dirname(args[2]), pattern = "*_blastout", full.names = TRUE)
poorseqs_blast <- blastoutfiles %>%
  map(readLines) %>%
  reduce(rbind)

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
