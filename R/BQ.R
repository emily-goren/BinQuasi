#' 
#' Call peaks in replicated ChIP-seq data using BinQuasi
#' 
#' @description Use the BinQuasi algorithm to call peaks using ChIP-seq data with biological replicates.
#' 
#' @param dir Directory where the sorted bam files (and their corresponding bam 
#' indices) are saved.
#' @param ChIP.files File names (with file extensions) of the ChIP sample files 
#'   in sorted bam format.
#' @param control.files File names (with file extensions) of the control/input 
#'   sample files in sorted bam format.
#' @param alpha The desired significance threshold used to call peaks. Must be 
#'   in (0, 0.5).
#' @param bin.size Window size (constant across all samples) used to generate a 
#'   partition for counts. If \code{NULL}, it will be estimated based on 
#'   Shimazaki and Shinomoto (2007).
#' @param frag.length Average length of the ChIP fragments in each sample 
#'   provided. Reads are extended to this length in the 5'-to-3' direction. If 
#'   \code{NULL}, cross correlation will be used to estimate the fragment 
#' @param minimum.count The count threshold used for filtering out windows with 
#'   sparse counts. Any genomic window with a total count, across all samples, 
#'   less than this value will be removed.
#' @param Model Must be one of \code{"Poisson"} or \code{"NegBin"}, specifying use of a 
#'   quasi-Poisson or quasi-negative binomial model, respectively.
#' @param print.progress logical. If \code{TRUE}, updates are provided regarding
#'   which window (row number) is being analyzed. Updates occur frequently to 
#'   start then eventually occur every 5000 windows.
#' @param method Must be one of \code{"QL"}, \code{"QLShrink"}, or \code{"QLSpline"}, 
#' specifying which method of Lund, Nettleton, McCarthy and Smyth (2012) should be used to 
#'   compute p-values.
#' @param p.window.adjust FDR control method applied to the windows. Must be 
#' either \code{"BH"} or \code{"BY"} to specify the procedure of Benjamini-Hochberg 
#' or Benjamini-Yekutieli, respectively.
#' @param Dispersion Must be one of \code{"Deviance"} or \code{"Pearson"}, specifying which 
#'   type of estimator should be used for estimating the quasi-likelihood 
#'   dispersion parameters.
#' @param log.offset A vector of log-scale, additive factors used to adjust 
#'   estimated log-scale means for differences in library sizes across samples. 
#'   Commonly used offsets include \code{log.offset=log(colSums(counts))} and 
#'   \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}. 
#'   If \code{NULL}, the later offset is used.
#' @param NBdisp Used only when \code{Model="NegBin"}. Must be one of \code{"trend"}, 
#'   \code{"common"}, or a vector of non-negative real numbers with length equal to 
#'   \code{nrow(counts)}. Specifying \code{NBdisp="trend"} or 
#'   \code{NBdisp="common"} will use {\link{estimateGLMTrendedDisp}} or 
#'   \code{\link{estimateGLMCommonDisp}}, respectively, from the package 
#'   \code{\link{edgeR}} to estimate negative binomial dispersion parameters for each 
#'   window.  Estimates obtained from other sources can be used by entering 
#'   \code{NBdisp} as a vector containing the negative binomial dispersion value
#'   to use for each window when fitting the quasi-likelihood model.
#' @param bias.fold.tolerance A numerical value no smaller than 1. If the bias 
#'   reduction of maximum likelihood estimates of (log) fold change is likely to
#'   result in a ratio of fold changes greater than this value, then bias 
#'   reduction will be performed on such windows. Setting 
#'   \code{bias.fold.tolerance=Inf} will completely disable bias reduction; 
#'   setting \code{bias.fold.tolerance=1} will always perform bias reduction. 
#'   See \code{\link{NBDev}} or \code{\link{PoisDev}} for details.
#'   
#' @details This function calls peaks in replicated ChIP-seq data.
#' 
#' @references Shimazaki and Shinomoto (2007)  "A method for selecting the bin 
#' size of a time histogram" \emph{Neural computation}, \bold{19}(6), 1503-27.
#' 
#' Ramachandran, Palidwor, Porter, and Perkins (2013) "MaSC: 
#' mappability-sensitive cross-correlation for estimating mean fragment length 
#' of single-end short-read sequencing data" \emph{Bioinformatics} \bold{29}(4),
#' 444-50.
#' 
#' Benjamini and Hochberg (1995) "Controlling the false discovery rate: a 
#' practical and powerful approach to multiple testing" \emph{Journal of the 
#' Royal Statistical Society Series B}, \bold{57}: 289-300.
#' 
#' Benjamini and Yekutieli (2001) "The control of the false discovery rate in
#' multiple testing under dependency" \emph{Annals of Statistics}. \bold{29}:
#' 1165-1188.
#' 
#' Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression
#' in RNA-sequence data using quasi-likelihood with shrunken dispersion 
#' estimates" \emph{SAGMB}, \bold{11}(5).
#' 
#' @author Emily Goren (\email{emily.goren@gmail.com})
#' 
#' @examples 
#' \dontrun{
#' # Fit a quasi-negative binomial model using all default settings.
#' fpath <- paste0(system.file(package = 'BinQuasi'), '/extdata/')
#' results <- BQ(fpath, ChIP.files = c('C1.bam', 'C2.bam'), control.files = c('I1.bam', 'I2.bam'))
#' head(results$peaks)
#' }
#' 
#' @return A list containing: 
#'   \item{peaks}{Dataframe of the 
#'   called peaks with columns for the start and end location, width, 
#'   chromosome, p-value, and q-value computed using the Benjamini and Hochberg 
#'   method.} 
#'   \item{bin.size}{The window width used to create the counts dataframe.}
#'   \item{fragment.length}{Vector of the fragment lengths used to extend the 
#'   reads in each sample.} 
#'   \item{filter}{The count threshold used to create the
#'   counts dataframe. Windows with counts below this value were removed.}
#'   
#' @export
#' 

BQ <- function(dir, ChIP.files, control.files, alpha=0.05,
               bin.size=NULL, frag.length=NULL,
               minimum.count=20, Model="NegBin",
               print.progress=TRUE, method="QLShrink",
               p.window.adjust="BY",
               Dispersion="Deviance", log.offset=NULL,
               NBdisp="trend", bias.fold.tolerance=1.10) {
  
  if (alpha <= 0 | alpha >= 0.5)
    stop("Please specify a significance level between 0 and 0.5")
  
  if(!(Dispersion %in% c("Deviance","Pearson")))
    stop("Unidentified Dispersion: Dispersion must be either 'Deviance' or 'Pearson'.")
  
  if(!(p.window.adjust %in% c("BY","BH")))
    stop("Unidentified method of window-level FDR control: p.window.adjust must be either 'BY' or 'BH'.")
  
  if(!(method %in% c("QL", "QLShrink", "QLSpline")))
    stop("Unidentified method for computing p-values: method must be either 'QL', 'QLShrink', or 'QLSpline'.")
  
  
  if (!(Model %in% c("NegBin", "Poisson")))
    stop("Unidentified Model: Model must be either 'NegBin' or 'Poisson'.")
  
  n.C <- length(ChIP.files); n.I <- length(control.files)
  cts <- count.table(dir, ChIP.files, control.files, bin.size, frag.length, minimum.count)
  message('Using bin size of ', cts$bin.size, ' bp \n')
  Len <- cts$fragment.length
  Lenname <- c(ChIP.files, control.files)
  for (i in 1:length(Len))
    message('Using estimated fragment length for ', Lenname[i], ' equal to ', Len[i], ' bp \n')
  counts <- cts$counts
  Y <- counts[ , 4:(3+n.C+n.I)]
  
  
  if (Model == "NegBin") {
    if (!NBdisp[1] %in% c("trend", "common") & length(NBdisp) != nrow(Y))
      stop("Unidentified NegBin Dispersion: NBdisp must be set as 'trend' or 'common' to estimate negative binomial dispersion from data using GLM edgeR (McCarthy et al., 2012),\n or it must be a vector providing negative binomial dispersion parameter value to use for each window.")
    
    if (length(NBdisp) == nrow(Y) & !is.numeric(NBdisp))
      stop("NBdisp contains non-numeric values.\n\tAll negative binomial dispersion parameters must be non-negative real numbers.")
    
    if (length(NBdisp) == nrow(Y) & any(NBdisp < 0))
      stop("NBdisp contains negative values.\nAll negative binomial dispersion parameters must be non-negative real numbers.")
  }
  
  
  design.matrix <- cbind(rep(1, n.C + n.I),
                         c(rep(1, n.C), rep(0, n.I)))
  fit <- QL.fit(counts = Y,
                design.matrix = design.matrix,
                chip.col.indicator = c(0,1),
                log.offset = log.offset,
                NBdisp = NBdisp,
                print.progress = print.progress,
                bias.fold.tolerance = bias.fold.tolerance,
                Model = Model)
  res <- QL.results(fit = fit, Dispersion = Dispersion)
  if (method == 'QL') {
    ps <- res$P.values$QL
  } else if (method == 'QLShrink') {
    ps <- res$P.values$QLShrink
  } else if (method == 'QLSpline') {
    ps <- res$P.values$QLSpline
  }
  called <- call.peaks(window.pvals = ps,
                       method = p.window.adjust,
                       start = counts$start,
                       end = counts$end,
                       chromosomes = counts$chr,
                       alpha = alpha)
  out <- list(peaks = called,
              bin.size = cts$bin.size,
              fragment.length = cts$fragment.length,
              filter = cts$filter)
  return(out)
}
