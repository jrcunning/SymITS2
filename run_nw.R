# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
# If two arguments not provided, return an error
if (length(args) < 2) {
  stop("must specify query sequences (.fasta) and reference sequences (.fasta)", call.=FALSE)
}

# Import query sequences
library(Biostrings)
otus <- readDNAStringSet(args[1])
names(otus) <- gsub(" .*$", "", names(otus))

# Import reference sequences
ref <- readDNAStringSet(args[2])

# import substitution matrix
EDNAFULL <- as.matrix(read.table("data/EDNAFULL", skip=8, header=T, row.names=1))

# Build function to find which reference sequences produces best alignment score using
# Needleman-Wunsch algorithm with default EMBOSS parameters for DNA sequences.
# -type="global-local" forces whole reference (pattern) sequence to be used, but subject/query may be 
# a concescutive subsequence. This is because queries may include slightly longer region of gene 
# than references because references were trimmed by o-smart-trim. This means that query sequence 
# can have more bases than the reference at beginning or end of the sequence (=end gaps in 
# alignment) without penalties.
get.best <- function(x) {
  scores <- pairwiseAlignment(subject=x, pattern=ref, type="global-local", scoreOnly=TRUE,
                           substitutionMatrix=EDNAFULL,
                           gapOpening=10.0, gapExtension=0.5)
  best <- cbind(names(x), names(ref)[which(scores==max(scores))])
  return(best)
}

# Set up for parallel computation
library(parallel)
cl <- makeCluster(detectCores() - 1)  # Initiate cluster
clusterEvalQ(cl, library(Biostrings))  # Make Biostrings library available to cluster
clusterExport(cl, c("EDNAFULL", "ref", "otus", "get.best"))  # Export objects to cluster

# Run function in parallel to get best hit for each sequence
best <- parSapply(cl, otus, FUN=get.best)
stopCluster(cl)  # Stop cluster

# Collect results in a data frame
results <- data.frame(otu=rep(names(best), sapply(best, length)),
                      hit=unlist(best), stringsAsFactors = FALSE)

# Align each sequence with its bets hits from database to get number of mismatches, indels, and % similarity
# Align each otu sequence to its best hits
psa <- pairwiseAlignment(subject=otus[results$otu], pattern=ref[results$hit], type="global-local",
                         substitutionMatrix=EDNAFULL,
                         gapOpening=10.0, gapExtension=0.5)

# Calculate percent similarity to reference with indels counting as single difference regardless of length
results <- within(results, {
  score <- score(psa)
  len <- width(pattern(psa))
  mis <- nmismatch(psa)
  ins <- unlist(lapply(as.list(insertion(psa)), length))
  del <- unlist(lapply(as.list(deletion(psa)), length))
  sim <- round((len-mis-ins-del)/len, 4) * 100
})

# Merge identical hits
results2 <- with(results, aggregate(hit ~ otu + sim + del + ins + mis + len + score, FUN=paste0, collapse=";"))

# Write results to .tsv file
write.table(results2, file=paste(dirname(args[1]), "/", "nw_tophits.tsv", sep=""), 
            sep="\t", quote=FALSE)

