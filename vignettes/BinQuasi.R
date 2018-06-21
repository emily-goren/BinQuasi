## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(BinQuasi)
fpath <- paste0(system.file(package = 'BinQuasi'), '/extdata/')
results <- BQ(fpath, 
              ChIP.files = c('C1.bam', 'C2.bam'), 
              control.files = c('I1.bam', 'I2.bam'))
head(results$peaks)

## ------------------------------------------------------------------------
?BQ

## ------------------------------------------------------------------------
# Sort peaks by p-value
opeaks <- results$peaks[order(results$peaks$P.val),]
# Name the peaks by rank
opeaks$name <- paste0('BQ_Peak_', 1:nrow(opeaks))
# Save as .bed file, setting the scores to be -log10(p-value)
bedout <- data.frame(chrom = opeaks$chr,
                     chromStart = opeaks$start,
                     chromEnd = opeaks$end,
                     name = opeaks$name,
                     score = -log10(opeaks$P.val),
                     strand = c(rep(".",  nrow(opeaks))))
head(bedout)
write.table(bedout, file="BinQuasiPeaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

