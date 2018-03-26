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

