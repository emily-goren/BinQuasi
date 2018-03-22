# BinQuasi

This package provides code to call peaks in ChIP-seq data with biological replicates using the algorithm of Goren, Liu, Wang and Wang (2018).

## Installation

The BinQuasi package for R can be installed from Github using devtools following the code below.

```R
install.packages("devtools")
library(devtools)
install_github("emilygoren/BinQuasi", args = "--preclean")
library(BinQuasi)
```

## Data Preprocessing

BinQuasi accepts sorted and indexed BAM files (note that it does not perform genome alignment of raw reads).


## Peak Calling

Once installed, BinQuasi calls peaks with the function "BQ()." Below is code to run BinQuasi with all default settings, where the sorted and indexed BAM files are stored in the directory "/Users/username/mybamfolder/" under the file names "chip1.bam", "chip2.bam" and "input1.bam", "input2.bam" for ChIP and input files, respectively. See the package documentation for information on changing the default settings.

```R
dir <- '/Users/username/mybamfolder/'
results <- BQ(dir, ChIP.files = c('chip1.bam', 'chip2.bam'), control.files = c('inp1.bam', 'inp2.bam'))
head(results$peaks)
```
