# SymITS2

This is a collection of scripts for bioinformatic analysis of *Symbiodinium* ITS2 data. To utilize these scripts, you will need to have working installations of [qiime](http://qiime.org), [GNU parallel](https://www.gnu.org/software/parallel/), [blast (>2.2.31)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [R](https://www.r-project.org), and the R packages [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)). Many of the scripts are designed to make use of parallel processing on capable systems.

## Scripts

- **fetch_seqs.sh**  
   This script downloads sequence data from NCBI given a list of named sequences with accession numbers (e.g., >C1_AF333515)
- **prep_ITS2db.sh**  
   This script takes downloaded *Symbiodinium* sequence data and prepares an ITS2 reference database by 1) trimming primers (using [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html)), 2) separating by clade, aligning (using [muscle](http://drive5.com/muscle/), and trimming sequences to equal length (using [oligotyping](https://github.com/merenlab/oligotyping)), 3) dereplicating (using [usearch](http://www.drive5.com/usearch/)).
- **merge_reads.sh**  
   This script parallelizes [illumina-utils](https://github.com/merenlab/illumina-utils) software to merge paired Illumina reads (.fastq)
- **qc_trim_reads.sh**  
   
- **otus_97_bysample.sh**  
- **otus_100.sh**  
- **run_nw.R**  
- **build_phyloseq.R**  
- **filter_notsym.R**  
