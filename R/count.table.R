# Cross correlation of for a given chromosome. See MaSC paper (Bioinformatics, 2013).
xcorr <- function(gdata) { # gdata is GRanges object for a single chromosome.
    chr <- gdata@seqnames@values
    B <- gdata@seqinfo@seqlengths[gdata@seqinfo@seqnames == chr]
    mus <- table(strand(gdata)) / B
    mu.f <-mus[1]
    mu.g <- mus[2]
    V.f <- mu.f * (1 - mu.f)
    V.g <- mu.g * (1 - mu.g)
    pos.strand <- gdata[strand(gdata) == '+']
    neg.strand <- gdata[strand(gdata) == '-']
    f <- start(pos.strand)
    g <- end(neg.strand)
    ds <- seq(50, 600, by = 5)
    rd <- lapply(ds, function(d) {
                     f.d <- f[f <= (B-d)]
                     g.d <- g[g >= (1+d)]
                     g.d <- g.d - d
                     fgsum <- length(intersect(f.d, g.d))
                     num <- fgsum / (B-d) - mu.f * mu.g
                     den <- sqrt(V.f * V.g)
                     return(num / den)
                 })
    d <- ds[which.max(unlist(rd))]
    return(d)
}



# Function to estimate the fragment length using cross correlation.
fragment.length <- function(dir, ChIP.files, control.files) {
    # Get bam file info.
    bamFiles <- paste0(dir, c(ChIP.files, control.files))
    n <- length(bamFiles)
    bfList <- BamFileList(bamFiles)
    samHeader <- lapply(bfList, function(x) scanBamHeader(path(x)))

    # Extract chromosomes.
    chromosomes <- lapply(samHeader, function(x) x[[1]]$targets)
    if (length(unique(chromosomes)) != 1)
        warning("The supplied bam files were not all aligned to the same chromosome set. Using chromosomes from the first bam file to estimate fragment length.")
    chromosomes <- chromosomes[[1]]
    chr.ranges <- GRanges(seqnames = names(chromosomes), ranges = IRanges(1, chromosomes))
    param <- lapply(chr.ranges, function(x) ScanBamParam(which=x))

    # Estimate fragment length for each chromosome within each replicate.
    fragLen <- lapply(1:n, function(i) { # Apply over replicates.
                           out <- lapply(param, function(x) { # Apply over chromosomes.
                                         rds <- readGAlignments(path(bfList[[i]]), param = x) # Convert bam to GAlignments.
                                         rds <- unique(GRanges(rds)) # Dedupe and convert to GRanges.
                                         est <- xcorr(rds)
                                         return(est)
                                     })
                           ests <- unlist(out)
                           frag <- mean(ests)
                           return(frag)
                       })

    fl <- round(unlist(fragLen), 0)

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
    samHeader <- lapply(bfList, function(x) scanBamHeader(path(x)))

    # Extract largest chromosome from each sample.
    chrs <- lapply(samHeader, function(x) {
                              t <- x[[1]]$targets
                              return(t[which.max(t)])
                          })
    chr.ranges <- lapply(chrs, function(x)  GRanges(seqnames = names(x), ranges = IRanges(1, x)))
    param <- lapply(chr.ranges, function(x) ScanBamParam(which=x))
    rds <- lapply(1:n, function(i) { # Apply over samples.
                      rds.ga <- readGAlignments(path(bfList[[i]]), param = param[[i]]) # Convert bam to GAlignments.
                      rds.gr <- unique(GRanges(rds.ga)) # Dedupe and convert to GRanges.
                      rds.out <- resize(rds.gr, frag.length[i]) # Extend reads to fragment length.
                      return(rds.out)
                  })
    min.bin <- min(frag.length)/2
    if (min.bin < 50)
        min.bin <- 50
    candidate.bins <- seq(min.bin, 1000, by = 10)
    est <- lapply(rds, function(x) {
                      w <- lapply(candidate.bins, function(d) cost.func(x, d))
                      best <- candidate.bins[which.min(unlist(w))]
                      return(best) })
    est <- min(unlist(est))
    return(est)
}




# Function to construct table of counts.
count.table <- function(dir, ChIP.files, control.files, bin.size = NULL, frag.length = NULL, minimum.count = 20) {

    # If shift size is not specified or wrong length, estimate using cross correlation.
    N <- length(c(ChIP.files, control.files))
    if (is.null(frag.length)) {
        message('Fragment length not provided. Estimating fragment length using cross correlation.')
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
        message('Bin size not provided. Estimating bin size.')
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
                            ext = frag.length, shift = 0,
                            filter = minimum.count, bin = FALSE)
    cts <- assay(cts.all)
    colnames(cts) <- c(ChIP.files, control.files)

    start <- rowRanges(cts.all)@ranges@start
    end <- start + rowRanges(cts.all)@ranges@width -1
    chr <- as.character(seqnames(rowRanges(cts.all)))

    # Counts file.
    cts.out <- data.frame(chr, start, end, cts)
    out <- list(counts = cts.out, bin.size = bin.size, bin.size, fragment.length = frag.length, filter = minimum.count)
    return(out)
}
