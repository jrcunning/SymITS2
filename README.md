# SymITS2

This is a collection of scripts for bioinformatic analysis of *Symbiodinium* ITS2 data. To utilize these scripts, you will need to have working installations of [qiime](http://qiime.org), [GNU parallel](https://www.gnu.org/software/parallel/), [blast (>2.2.31)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [R](https://www.r-project.org), [gawk](https://www.gnu.org/software/gawk/manual/html_node/This-Manual.html) (can be downloaded on Linux machines via 'sudo apt install gawk'), and the R packages [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)).

## Scripts

- **fetch_seqs.sh**  
   This script downloads sequence data from NCBI given a list of named sequences with accession numbers (e.g., >C1_AF333515)
- **prep_ITS2db.sh**  
   This script takes downloaded *Symbiodinium* sequence data and prepares an ITS2 reference database by 1) trimming primers (using [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html)), 2) separating by clade, aligning (using [muscle](http://drive5.com/muscle/), and trimming sequences to equal length (using [oligotyping](https://github.com/merenlab/oligotyping)), 3) dereplicating (using [usearch](http://www.drive5.com/usearch/)).
- **merge_reads.sh**  
   This script parallelizes [illumina-utils](https://github.com/merenlab/illumina-utils) software to merge paired Illumina reads (.fastq)
- **qc_trim_reads.sh**  
   This script adds QIIME labels, removes chimeric sequences, and trims primer sequences using [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html).
- **otus_97.sh**  
   This script 1) clusters the entire sequence dataset at 97% similarity implemented in QIIME using uclust, 2) removes singletons, 3) picks the most abundant sequence variants as the representative set, and 4) generates .biom and .tsv OTU tables.
- **otus_100.sh**  
   This script 1) clusters the entire sequence dataset at 100% similarity implemented in QIIME using uclust, 2) removes singletons, 3) picks the most abundant sequence variants as the representative set, and 4) generates .biom and .tsv OTU tables.
- **otus_97_bysample.sh**  
   This script 1) clusters the sequences from each sample independently at 97% similarity in parallel, 2) picks the most abundant sequence variants as the representative set, 3) merges OTUs that are 100% identical across samples, 4) removes singletons, and 5) generates .biom and .tsv OTU tables.  
- **run_nw.R**  
   This script takes two command line arguments (query sequences and reference sequences in .fasta format) and performs global aligment (in parallel) of each query against the reference database to identify the best match. The alignment score, percent similarity, and number of matches, mismatches, and indels are calculated for each query's best match, and the results are written to a .tsv file.
- **build_phyloseq.R**  
   This script takes five command line arguments (1: taxonomic assignment output from run_nw.R, 2: sample metadata, 3: OTU table, 4: duplicate reference taxa names, 5: output filename) and builds an R [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html) object, saving as .RData file.
- **filter_notsym.R**  
   This script imports a phyloseq object, identifies OTUs whose best taxonomic assignment from the reference database was \<90% similar, and BLASTs the corresponding sequences to the NCBI nr database. If the top hit returned does not contain the string "Symbiodinium", then the OTU is assumed to not be *Symbiodinium*. A filtered phyloseq object without these OTUs is returned and saved as .RData file.
   takes three command line arguments (1: phyloseq object .RData, 2: sequence data .fasta, 3: output filename).
