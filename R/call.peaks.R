
get.simes <- function(bins, regions, pvals) { ## Supply bins, regions as GRanges objects
    simes <- function(p) min(length(p) * p / rank(p))
    hits <- findOverlaps(bins, regions)
    nbins <- length(regions)
    out <- lapply(1:nbins, function(j) {
                        idx <- hits@from[hits@to == j]
                        p.simes <- simes(pvals[idx])
                        return(p.simes)
                    })
    return(unlist(out))
}

call.peaks <- function(window.pvals, start, end, chromosomes, alpha=0.05) {
    if (alpha <= 0 | alpha >= 0.5)
        stop("Please specify a significance level between 0 and 0.5")
    K <- length(window.pvals)
    if (!(K == length(start) & K == length(end) & K == length(chromosomes)))
        stop("The length of the p-value vector, start locations, end locations, and chromosomes must match.")
    q <- p.adjust(window.pvals, method = 'fdr')
    sig <- q < alpha # Which windows are significant?
    if (sum(sig) == 0)
        stop("No windows are significant at the specified alpha. No peaks called.")
    bins <- GRanges(chromosomes, IRanges(start = start, end = end))
    regions <- reduce(GRanges(chromosomes[sig], IRanges(start = start[sig], end = end[sig]))) # candidate regions
    out <- data.frame(start = start(regions),
                      end = end(regions),
                      width = width(regions),
                      chr = as.character(seqnames(regions)),
                      P.val = get.simes(bins, regions, window.pvals))
    out$Q.val <- p.adjust(out$P.val, method = 'fdr')
    called <- out$Q.val < alpha
    return(out[called, ])
}
