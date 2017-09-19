% File src/library/QuasiSeq/man/BQ.Rd
\name{BQ}
\alias{BQ}
\title{Call Peaks in Replicated ChIP-seq Data via BinQuasi}
\description{Use the BinQuasi algorithm to call peaks using ChIP-seq data with
biological repliates.}
\usage{BQ(dir, ChIP.files, control.files, alpha=0.05,
          bin.size=NULL, frag.length=NULL, minimum.count=20,
          Model="NegBin", print.progress=TRUE, method="QLShrink",
          Dispersion="Deviance", log.offset=NULL, NBdisp="trend", bias.fold.tolerance=1.10)}
\arguments{
  \item{dir}{Directory where the sorted bam files are saved.}
  \item{ChIP.files}{File names (with file extensions) of the ChIP sample
  files in sorted bam format.}
  \item{control.files}{File names (with file extensions) of the
    input/control sample files in sorted bam format.}
   \item{alpha}{The desired significance threshold used to call
     peaks. Must be in (0, 0.5).}
  \item{bin.size}{Window size, constant across all samples, used to generate a partition for
    counts. If \code{NULL}, it will be estimated based on Shimazaki
    and Shinomoto (2007).}
  \item{frag.length}{Average length of the ChIP fragments in each sample
    provided. Reads are extended to this length in the 5'-to-3' direction. If \code{NULL}, cross correlation will be used to
    estimate the fragment lengths. See Ramachandran, Palidwor, Porter, and Perkins (2013).}
  \item{minimum.count}{The count threshold used for filtering out
    windows with sparse counts. Any genomic window with counts less than
    this value will be removed.}
	\item{Model}{Must be one of "Poisson" or "NegBin", specifying
	  use of a quasi-Poisson or quasi-negative binomial model,
	  respectively.}
	\item{print.progress}{logical. If \code{TRUE}, updates are
	  provided regard what window (row number) is being analyzed.
	  Updates occur frequently to start then eventually occur every
	  5000 windows.}
	\item{method}{Must be one of "QL", "QLShrink", or "QLSpline",
	  specifing which method of Lund, Nettleton, McCarthy and Smyth
	  (2012) should be used to compute p-values.}
	\item{Dispersion}{Must be one of "Deviance" or "Pearson",
	  specifying which type of estimator should be used for
	  estimating the quasi-likelihood dispersion parameters.}
	\item{log.offset}{A vector of log-scale, additive factors used
     to adjust estimated log-scale means for differences in library
     sizes across samples.  Commonly used offsets include
     \code{log.offset=log(colSums(counts))} and
     \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}.  If \code{NULL}, the later offset is used.}
     \item{NBdisp}{Used only when \code{Model="NegBin"}. Must be one of "trend", "common" or a vector of non-negative real numbers with length equal to \code{nrow(counts)}. Specifying \code{NBdisp="trend"} or \code{NBdisp="common"} will use {\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}}, respectively, from the package \code{edgeR} to estimate negative binomial dispersion parameters for each window.  Estimates obtained from other sources can be used by entering \code{NBdisp} as a vector containing the negative binomial dispersion value to use for each window when fitting the quasi-likelihood model.}
	\item{bias.fold.tolerance}{A numerical value no smaller than
	  1. If the bias reduction of maximum likelihood estimates of
	  (log) fold change is likely to result in a ratio of fold
	  changes greater than this value, then bias reduction will be
	  performed on such genes. Setting
	  \code{bias.fold.tolerance=Inf} will completely disable bias
	  reduction;  setting \code{bias.fold.tolerance=1} will always
	  perform bias reduction.  See \code{\link{NBDev}} or \code{\link{PoisDev}} for
	  details.}
}
\value{list containing:
  \item{peaks}{dataframe of the called peaks with
       columns for the start and end locations, width,
       chromosome, p-value, and q-value computed using the Benjamini and Hochberg method.}
     \item{bin.size}{the bin.size used to create the counts dataframe.}
     \item{fragment.length}{vector of the fragment lengths used to
       extend the reads in each sample.}
     \item{filter}{the count threshold used to create the counts
       dataframe. Windows with counts below this value were removed.}
     }
\details{
This function calls peaks in replicated CHiP-seq data using the methods of Goren, Wang, and Liu (to
be submitted).}

\references{
Shimazaki and Shinomoto (2007)  "A method for selecting the bin size of a
time histogram" \emph{Neural computation}, \bold{19}(6), 1503-27.

Ramachandran, Palidwor, Porter, and Perkins (2013)
"MaSC: mappability-sensitive cross-correlation for estimating mean fragment length of single-end short-read sequencing data" \emph{Bioinformatics} \bold{29}(4), 444-50.

Benjamini and Hochberg (1995)
"Controlling the false discovery rate: a practical and powerful approach to multiple testing" \emph{Journal of the Royal Statistical Society Series B}, \bold{57}: 289-300.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" \emph{SAGMB}, \bold{11}(5).
}

\author{Emily Goren \email{emily.goren@gmail.com}}
\examples{
\dontrun{
# Fit a quasi-Poisson model to 4 samples.
peaks <- BQ(dir = '/Users/username/mybamfolder/',
            ChIP.files = c('chip1.bam', 'chip2.bam'),
            control.files = c('inp1.bam', 'inp2.bam'),
            bin.size = 60, frag.length = c(101, 300, 150, 10),
            minimum.count = 20, Model = "Poisson"))
}
}