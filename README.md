# BinQuasi

This package provides code to call peaks in ChIP-seq data with biological replicates using the BinQuasi algorithm.

## Installation

The BinQuasi package for R can be installed from Github using devtools following the code below.

```R
devtools::install_github("emilygoren/BinQuasi", args = "--preclean", build_vignettes = TRUE)
library(BinQuasi)
```

## Data Preprocessing

BinQuasi accepts sorted and indexed BAM files (note that it does not perform genome alignment of raw reads). If your BAM files are not indexed and sorted, we recommend using [samtools](http://www.htslib.org/).


## Peak Calling

Once installed, BinQuasi calls peaks with the function "BQ()." Below is code to run BinQuasi with all default settings, where the sorted and indexed BAM files are stored in the directory specified by "fpath" under the file names "C1.bam", "
C2.bam" and "I1.bam", "I2.bam" for ChIP and input files, respectively. 

```R
fpath <- paste0(system.file(package = 'BinQuasi'), '/extdata/')
results <- BQ(fpath, 
              ChIP.files = c('C1.bam', 'C2.bam'), 
              control.files = c('I1.bam', 'I2.bam'))
head(results$peaks)
```

See the package documentation for information on changing the default settings.
```R
?BQ
```