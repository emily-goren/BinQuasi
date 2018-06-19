simes <- function(p) min(length(p) * p / rank(p))

get.simes <- function(bins, regs, pvals) { 
  ## Supply bins, regions as GRanges objects
  hits <- findOverlaps(bins, regs)
  out <- sapply(seq_along(regs), function(j) {
    idx <- hits@from[hits@to == j]
    p.simes <- simes(pvals[idx])
    return(p.simes)
  })
  return(out)
}

#' 
#' Call peaks from a list of window-level p-values
#' 
#' @description Call peaks from a list of p-values corresponding to window-level
#'   tests on a genomic partition of ChIP-seq counts. Used within the main peak
#'   calling function, \code{\link{BQ}}.
#'   
#' @details After correcting for multiple testing using the adjustment specified
#'   by \code{method}, windows that are significant according to the threshold 
#'   \code{alpha} are merged if adjacent and retained as candidate regions. 
#'   Simes' procedure is used to combine the window-level p-values in each 
#'   candidate region into a region-level p-value. The Benjamini-Hochberg 
#'   procedure is applied to the resulting candidate regions and those that 
#'   exceed the significance threshold \code{alpha} are returned as peaks.
#'   
#' @param window.pvals Vector of p-values, with each element corresponding to a 
#'   window of a genomic partition. Typically obtained from the 
#'   \code{\link{QL.fit}} and \code{\link{QL.results}} functions.
#' @param method Correction method applied to \code{window.pvals}. Must be one 
#'   of \code{"BH"}, \code{"BY"}, or \code{"none"} to specify Benjamini-Hochberg, 
#'   Benjamini-Yekutieli, or no adjustment, respectively.
#' @param start Vector of the genomic start locations corresponding to the 
#'   supplied p-values.
#' @param end Vector of the genomic end locations corresponding to the supplied 
#'   p-values.
#' @param chromosomes Vector of the chromosome names corresponding to the 
#'   supplied p-values.
#' @param alpha The desired significance threshold in (0, 0.5).
#'   
#' @return The called peaks as a dataframe with variables: 
#'   \item{start}{Genomic start locations of the called peaks.} 
#'   \item{end}{Genomic end locations of the called peaks.} 
#'   \item{width}{Width of the called peaks.} 
#'   \item{chr}{Chromosomes of the called peaks.} 
#'   \item{P.val}{p-values of the  called peaks (aggregated from the windows 
#'   comprising the peak using Simes' procedure).} 
#'   \item{Q.val}{q-values of the called peaks (computing using the
#'   Benjamini-Hochberg procedure).}
#'   
#' @author Emily Goren (\email{emily.goren@gmail.com}).
#'   
#' @references Benjamini and Hochberg (1995) "Controlling the false discovery 
#'   rate: a practical and powerful approach to multiple testing" \emph{Journal
#'   of the Royal Statistical Society Series B}, \bold{57}: 289-300.
#'   
#'   Benjamini and Yekutieli (2001) "The control of the false discovery rate in 
#'   multiple testing under dependency" \emph{Annals of Statistics}. \bold{29}: 
#'   1165-1188.
#'   
#'   Simes (1986) "An improved Bonferroni procedure for multiple tests of 
#'   significance" \emph{Biometrika}, \bold{73}(3): 751-754.
#'   
#' @examples
#' # Example for a single chromosome.
#' start <- seq(1, 1e6, by = 200)
#' end <- start + 200 - 1
#' chromosomes <- rep('chr1', length(start))
#' p <- c(runif(length(start) - 10), rep(1e-12, 10))
#' called <- call.peaks(p, "BH", start, end, chromosomes)
#' called
#' 
#' @export
#' 

call.peaks <- function(window.pvals, method=c("BY", "BH", "none"), start, end, chromosomes, alpha=0.05) {
  if (alpha <= 0 | alpha >= 0.5)
    stop("Please specify a significance level between 0 and 0.5")
  K <- length(window.pvals)
  if (!(K == length(start) & K == length(end) & K == length(chromosomes)))
    stop("The length of the p-value vector, start locations, end locations, and chromosomes must match.")
  q <- p.adjust(window.pvals, method = method)
  sig <- q < alpha # Which windows are significant?
  if (sum(sig) == 0) {
    message("No windows are significant at the specified alpha using this binwidth and count filtering. No peaks called.")
    out <- data.frame(start = NA, end =  NA, width = NA, chr = NA, P.val = NA, Q.val = NA)
    called <- FALSE
  } else {
    bins <- GRanges(chromosomes, IRanges(start = start, end = end))
    regs <- reduce(GRanges(chromosomes[sig], IRanges(start = start[sig], end = end[sig]))) # candidate regions
    out <- data.frame(start = regs@ranges@start,
                      end = regs@ranges@start + regs@ranges@width - 1,
                      width = regs@ranges@width,
                      chr = as.character(regs@seqnames),
                      P.val = get.simes(bins, regs, window.pvals))
    out$Q.val <- p.adjust(out$P.val, method = 'BH')
    called <- (out$Q.val < alpha)
  }
  ans <- out[called,]
  return(ans)
}