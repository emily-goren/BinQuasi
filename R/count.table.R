
# Cross correlation of for a given chromosome. See MaSC paper (Bioinformatics, 2013).
xcorr <- function(gdata) { # gdata is GRanges object for a single chromosome.
  chr <- gdata@seqnames@values
  B <- gdata@seqinfo@seqlengths[gdata@seqinfo@seqnames == chr]
  mus <- table(gdata@strand) / B
  mu.f <-mus['+']
  mu.g <- mus['-']
  V.f <- mu.f * (1 - mu.f)
  V.g <- mu.g * (1 - mu.g)
  pos.strand <- gdata[gdata@strand == '+']
  neg.strand <- gdata[gdata@strand == '-']
  f <- pos.strand@ranges@start
  g <- neg.strand@ranges@start + neg.strand@ranges@width - 1
  ds <- seq(50, 600, by = 5)
  rd <- sapply(ds, function(d) {
    f.d <- f[f <= (B-d)]
    g.d <- g[g >= (1+d)]
    g.d <- g.d - d
    fgsum <- length(intersect(f.d, g.d))
    num <- fgsum / (B-d) - mu.f * mu.g
    den <- sqrt(V.f * V.g)
    return(num / den)
  })
  d <- ds[which.max(rd)]
  return(d)
}



# Function to estimate the fragment length using cross correlation.
fragment.length <- function(dir, ChIP.files, control.files) {
  # Get bam file info.
  bamFiles <- paste0(dir, c(ChIP.files, control.files))
  n <- length(bamFiles)
  bfList <- BamFileList(bamFiles)
  samHeader <- lapply(bfList, function(x) scanBamHeader(x$path))
  
  # Extract chromosomes.
  chromosomes <- lapply(samHeader, function(x) x[[1]]$targets)
  if (length(unique(chromosomes)) != 1)
    warning("The supplied bam files were not all aligned to the same chromosome set. Using chromosomes from the first bam file to estimate fragment length.")
  chromosomes <- chromosomes[[1]]
  chr.ranges <- GRanges(seqnames = names(chromosomes), ranges = IRanges(1, chromosomes))
  param <- lapply(seq_along(chr.ranges), function(i) ScanBamParam(which=chr.ranges[i]))
  
  # Estimate fragment length for each chromosome within each replicate.
  fragLen <- sapply(1:n, function(i) { # Apply over replicates.
    ests <- sapply(param, function(x) { # Apply over chromosomes.
      rds <- readGAlignments(bfList[[i]]$path, param = x) # Convert bam to GAlignments.
      rds <- unique(GRanges(rds)) # Dedupe and convert to GRanges.
      est <- xcorr(rds)
      return(est)
    })
    frag <- mean(ests)
    return(frag)
  })
  
  fl <- round(fragLen, 0)
  
  return(fl)
}


# Functions to estimate bin width. See Shimazaki and Shinomoto 2007.
cost.func <- function(gdata, d) {
  chr <- gdata@seqnames@values
  B <- gdata@seqinfo@seqlengths[gdata@seqinfo@seqnames == chr]
  start <- seq(1, B, by = d) # Window start locations.
  end <- start + d - 1 # Window end locations.
  bins <- IRanges(start = start, end = end)
  k <- countOverlaps(GRanges(seqnames = chr, ranges = bins), gdata) # Bin counts.
  kbar <- mean(k)
  N <- sum(k)
  v <- (sum((k - kbar)^2)) / length(k)
  num <-  2*kbar - v
  den <- (N*d)^2
  cost <- -(log(abs(num)) - log(den))
  return(cost)
}

bin.width <- function(dir, ChIP.files, frag.length) {
  # Get bam file info.
  bamFiles <- paste0(dir, ChIP.files)
  n <- length(bamFiles)
  bfList <- BamFileList(bamFiles)
  samHeader <- lapply(bfList, function(x) scanBamHeader(x$path))
  
  # Extract largest chromosome from each sample.
  # chrs <- lapply(samHeader, function(x) {
  #   t <- x[[1]]$targets
  #   return(t[which.max(t)])
  # })
  
  # Extract chromosomes from each sample.
  chrs <- lapply(samHeader, function(x) x[[1]]$targets)
  chr.ranges <- lapply(chrs, function(x)  GRanges(seqnames = names(x), ranges = IRanges(1, x)))
  param <- lapply(chr.ranges, function(x) ScanBamParam(which=x))
  rds <- lapply(1:n, function(i) { # Apply over samples.
    rds.ga <- readGAlignments(bfList[[i]]$path, param = param[[i]]) # Convert bam to GAlignments.
    rds.gr <- unique(GRanges(rds.ga)) # Dedupe and convert to GRanges.
    rds.out <- resize(rds.gr, frag.length[i]) # Extend reads to fragment length.
    return(rds.out)
  })
  min.bin <- min(frag.length)/2
  if (min.bin < 50)
    min.bin <- 50
  candidate.bins <- seq(min.bin, 1000, by = 10)
  chr.all <- unique(sapply(chrs, names))
  est.all <- lapply(chr.all, function(thischr) {
    lapply(rds, function(x) {
      thisxchr <- seqnames(x)
      if (thischr %in% thisxchr) {
          w <- lapply(candidate.bins, function(d) cost.func(x[thisxchr == thischr], d))
          best <- candidate.bins[which.min(unlist(w))]
      } else {
        best <- Inf
      }
      return(best) })
  })
  est.all <- sapply(est.all, function(e) min(unlist(e)))
  names(est.all) <- chr.all
  est <- min(est.all)
  return(est)
}

#' 
#' Create a matrix of ChIP-seq count data
#' 
#' @description Create a matrix of ChIP-seq count data from sorted bam files
#'   using a non-overlapping genomic partition. Used within the main peak calling
#'   function, \code{\link{BQ}}.
#'   
#' @param dir Directory where the sorted bam files (and their corresponding 
#' bam indices) are saved.
#' @param ChIP.files File names (with file extensions) of the ChIP sample files
#'   in sorted bam format.
#' @param control.files File names (with file extensions) of the input/control
#'   sample files in sorted bam format. 
#' @param bin.size Window size, constant across
#'   all samples, used to generate a non-overlapping partition for counts. If
#'   \code{NULL}, an estimate will be used (see details).
#' @param frag.length Average length of the ChIP fragments in each sample
#'   provided. Reads are extended to this length from their 3' ends. If
#'   \code{NULL}, cross correlation will be used to estimate the fragment length
#'   of each sample (see details).
#' @param minimum.count The count threshold used for filtering out windows with
#'   sparse counts. Any genomic window with counts less than this value across
#'   all samples will be removed.
#'   
#' @return A list containing: 
#' \item{counts}{Data frame with rows corresponding
#'   to genomic windows and columns for the chromosomes, start and end
#'   locations, as well as a column for the counts of each sample.} 
#'   \item{bin.size}{The bin size used to create the genomic partition.} 
#'   \item{fragment.length}{Vector of the fragment lengths used to extend the
#'   reads in each sample.} 
#'   \item{filter}{Count threshold used to create the
#'   counts data frame. Windows with counts summed across all samples that fall
#'   below this value were removed.}
#'   
#' @details This function creates a count table of ChIP sequencing data
#'   (supplied as sorted bam files) using a non-overlapping partition across 
#'   the genome.
#'   
#'   The fragment length (if not provided) is estimated using the 
#'   cross-correlation method of Ramachandran et al (2013). A fragment length
#'   is estimated for each sample, after removing duplicate reads, by taking the
#'   average over all chromosomes in the sample. Estimation is performed at 5 bp
#'   resolution and restricted to a minimum fragment length of 50 bp and maximum
#'   of 600 bp.
#'   
#'   The bin size (if not provided) is selected using a procedure by Shimazaki
#'   and Shinomoto (2007) based on minimizing the mean-integrated squared error
#'   for a time-dependent Poisson point process. This procedure is applied to
#'   each ChIP sample (at 5 bp resolution, restricted to a minimum of 50 bp and
#'   maximum of 1000 bp), and the minimum across all ChIP samples is returned as
#'   the bin size.
#'   
#'   For a given sample and window, the count is determined as the number of
#'   fragments overlapping the window.
#'   
#' @references Shimazaki and Shinomoto (2007)  "A method for selecting the bin
#' size of a time histogram" \emph{Neural computation}, \bold{19}(6), 1503-27.
#' 
#' Ramachandran, Palidwor, Porter,  and Perkins (2013) "MaSC:
#' mappability-sensitive cross-correlation for estimating mean fragment length
#' of single-end short-read sequencing data" \emph{Bioinformatics} \bold{29}(4),
#' 444-50.
#' 
#' @author Emily Goren (\email{emily.goren@gmail.com}).
#'   
#' @examples
#' \dontrun{
#' fpath <- paste0(system.file(package = 'BinQuasi'), '/extdata/')
#' d <- count.table(dir = fpath,
#'                  ChIP.files = c('C1.bam', 'C2.bam'),
#'                  control.files = c('I1.bam', 'I2.bam'),
#'                  bin.size = 60, frag.length = c(101, 300, 150, 10),
#'                  minimum.count = 20)
#'                  head(d$counts)
#' }
#' 
#' 
#' @export
#' 
count.table <- function(dir, ChIP.files, control.files, bin.size = NULL, frag.length = NULL, minimum.count = 20) {
  
  # If shift size is not specified or wrong length, estimate using cross correlation.
  N <- length(c(ChIP.files, control.files))
  if (is.null(frag.length)) {
    message('Fragment length not provided. Estimating fragment length using cross correlation... please wait...')
    frag.length <- fragment.length(dir, ChIP.files, control.files)
  }
  if (length(frag.length) == 1) {
    frag.length <- rep(frag.length, N)
  }
  if (length(frag.length) != N) {
    message('Fragment length vector differs in length from number of samples. Estimating fragment length using cross correlation.')
    frag.length <- fragment.length(dir, ChIP.files, control.files)
  }
  if (any(round(frag.length, 0) != frag.length)) {
    message('Fragment length(s) are not integer valued. Rounding to nearest integer.')
    frag.length <- round(frag.length, 0)
  }
  if (any(frag.length < 0)) {
    message('Negative fragment length(s) provided. Estimating fragment length using cross correlation.')
    frag.length <- fragment.length(dir, ChIP.files, control.files)
  }
  
  #########
  # If binsize is not specified or not a reasonable value, estimate using cost function.
  if (is.null(bin.size)) {
    message('Bin size not provided. Estimating bin size... please wait...')
    bin.size <- bin.width(dir, ChIP.files, frag.length)
  }
  if (bin.size < 50) {
    message('The bin size is too small (less than 50bp). Estimating bin size.')
    bin.size <- bin.width(dir, ChIP.files, frag.length)
  }
  if (bin.size > 1000) {
    message('The bin size is too large (greater than 1000bp). Estimating bin size.')
    bin.size <- bin.width(dir, ChIP.files, frag.length)
  }
  
  if (bin.size != floor(bin.size))  {
    message('The bin size specified is not integer valued. Rounding to nearest integer.')
    bin.size <- round(bin.size, 0)
  }
  
  # Read in all bam files and create count table.
  bam.files <- paste0(dir, c(ChIP.files, control.files))
  cts.all <- windowCounts(bam.files, spacing = bin.size,
                          ext = list(frag.length, NA), shift = 0,
                          filter = minimum.count, bin = FALSE)
  cts <- assay(cts.all)
  colnames(cts) <- c(ChIP.files, control.files)
  
  cts.rr <- rowRanges(cts.all)
  start <- cts.rr@ranges@start
  end <- start + cts.rr@ranges@width - 1
  chr <- as.character(cts.rr@seqnames)
  
  # Counts file.
  cts.out <- data.frame(chr, start, end, cts)
  out <- list(counts = cts.out, 
              bin.size = bin.size, 
              fragment.length = frag.length, 
              filter = minimum.count)
  return(out)
}
