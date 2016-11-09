# Script to build phyloseq object

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
# If two arguments not provided, return an error
if (length(args) < 4) {
  stop("must specify 1: tax data (nw_tophits.tsv); 2: sample data (mapping_file.txt); 3: OTU table (.tsv); 4: output name", call.=FALSE)
}

# Import taxonomic assignment data from nw
read.nw <- function(file) {
  nw <- read.table(file, stringsAsFactors=FALSE)
  nw <- nw[order(nw$otu),]
  return(nw)
}

tax <- read.nw(args[1])

# Deal with identical taxonomic assignments (because some reference sequences are not unique...)
# Create names for groups of identical sequences based on members of group...
ident <- readLines("data/ITS2db_trimmed_notuniques_otus.txt")
ident <- gsub("denovo[0-9]*\t", "", ident)
ident <- strsplit(ident, split="\t")
ident2 <- lapply(ident, str_match_all, pattern="[A-I]{1}[0-9]{1,3}.*[_/]")
ident2 <- lapply(ident2, function(g) gsub(pattern="\\..*_$", "", x=unlist(g)))
ident2 <- lapply(ident2, function(g) gsub(pattern="_$", "", g))
ident2 <- lapply(ident2, function(g) unlist(strsplit(g, split="/")))
subtypes <- lapply(ident2, function(x) levels(factor(unlist(x))))
subtypes <- lapply(subtypes, function(s) paste(paste(s, collapse="/"), "_multiple", sep=""))
names(ident) <- subtypes
ident <- melt(do.call(rbind, ident))
ident <- unique(ident[order(ident[,1]), c(3,1)])

# Replace any sequence name in taxonomy assignment that is a member of a group of identical sequences with the name of the group
collapse.idents <- function(df) {
  within(df, {
    for (i in 1:nrow(ident)) {
      hit <- gsub(ident[i,1], ident[i,2], hit)
    }
  })
}

tax <- collapse.idents(tax)
tax <- with(tax, tax[order(otu, -sim), ])
tax <- tax[!duplicated(tax$otu), ]
rownames(tax) <- tax$otu

# Import sample data
sam <- sample_data(read.delim(args[2], sep="\t", header=T, row.names=1))

# Import OTU tables
otu <- otu_table(read.table(args[3], header=T, check.names=F, row.names=1,
                              skip=1, sep="\t", comment.char=""), taxa_are_rows=T)



# Build phyloseq objects
phy <- phyloseq(otu, tax_table(as.matrix(tax)), sam)

# Save phyloseq object
save(phy, file=file.path("data", paste0(args[4], ".RData")))

