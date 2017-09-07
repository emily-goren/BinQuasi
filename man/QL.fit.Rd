% File src/library/QuasiSeq/man/QL.fit.Rd
     \name{QL.fit}
     \alias{QL.fit}
     \title{Fit quasi-likelihood models to ChIP-seq data partitioned into a count matrix}
     \description{
       Fit a quasi-likelihood model, where the ChIP regression
     coefficient is constrained to be nonnegative, to ChIP-seq data
     partitioned into a count matrix using the methods detailed in
     Goren, Wang, and Liu (to be submitted).}
     \usage{
    QL.fit(counts,design.matrix, chip.col.indicator, log.offset=NULL,
           Model="NegBin", print.progress=TRUE, NBdisp="trend",
           bias.fold.tolerance=1.10, ...)
     }
     \arguments{
      \item{counts}{A matrix of integer counts obtained by counting
     reads accross a genomic partition. Each row contains observations
     from a single window. Each column contains observations from a single sample (either ChIP or control/input).}
    \item{design.matrix}{A design matrix for the full model, including a
     column that indicates whether the observation is a ChIP sample. The number
     of rows must be \code{ncol(counts)}.  The number of columns must be
     at least two, usually an intercept and an indicator whether the
     sample is ChIP (1) or input/control (0). Means are modeled with a
     log link function.}
     \item{chip.col.indicator}{A binary vector of length
     \code{ncol(design.matrix)} that indicates which column of the
     design.matrix corresponds to the ChIP indicator.}
	\item{log.offset}{A vector of log-scale, additive factors used
     to adjust estimated log-scale means for differences in library
     sizes across samples.  Commonly used offsets include
     \code{log.offset=log(colSums(counts))} and
     \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}.  If not provided, the later offset is used.}
	\item{Model}{Must be one of "Poisson" or "NegBin", specifying use of a quasi-Poisson or quasi-negative binomial model, respectively.}
	\item{print.progress}{logical. If \code{TRUE}, updates are provided regard what window (row number) is being analyzed.  Updates occur frequently to start then eventually occur every 5000 windows.}
      \item{NBdisp}{Used only when \code{Model="NegBin"}. Must be one of "trend", "common" or a vector of non-negative real numbers with length equal to \code{nrow(counts)}. Specifying \code{NBdisp="trend"} or \code{NBdisp="common"} will use {\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}}, respectively, from the package \code{edgeR} to estimate negative binomial dispersion parameters for each window.  Estimates obtained from other sources can be used by entering \code{NBdisp} as a vector containing the negative binomial dispersion value to use for each window when fitting quasi-likelihood model.}
	\item{bias.fold.tolerance}{A numerical value no smaller than 1. If the bias reduction of maximum likelihood estimates of (log) fold change is likely to result in a ratio of fold changes greater than this value, then bias reduction will be performed on such genes. Setting \code{bias.fold.tolerance=Inf} will completely disable bias reduction;  setting \code{bias.fold.tolerance=1} will always perform bias reduction.  See \code{\link{NBDev}} for details. }
	\item{...}{arguments to be passed to the function \code{\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}} from the package \code{edgeR}.}
	}

\value{list containing:
	\item{"LRT"}{matrix providing unadjusted two-sided likelihood
	ratio test statistics.  Each column contains statistics from a
	single hypothesis test that the ChIP coefficient is equal to
	zero versus not equal to zero, applied separately to each window.}
	\item{"phi.hat.dev"}{vector providing unshrunken, deviance-based
	estimates of quasi-dispersion (phi) for each window.}
	\item{"phi.hat.pearson"}{vector providing unshrunken, Pearson estimates of quasi-dispersion (phi) for each window.}
	\item{"mn.cnt"}{vector providing average count for each window.}
	\item{"den.df"}{denominator degrees of freedom. Equal to the number of samples minus the number of fitted parameters in the full model.}
	\item{"num.df"}{vector of numerator degrees of freedom for each test, computed as the difference in the number of fitted parameters between the full and reduced models for each test.}
	\item{"Model"}{Either "Poisson" or "NegBin", specifying which model (quasi-Poisson or quasi-negative binomial, respectively) was used.}
	\item{"nb.disp"}{only appears when \code{Model="NegBin"}.  Vector providing negative binomial dispersion parameter estimate used during fitting of quasi-negative binomial model for each gene.}
	\item{fitted.values}{matrix of fitted mean values without constraints.}
	\item{coefficients}{matrix of estimated coefficients for each
	window.  Note that these are given on the log
	scale. (i.e. intercept coefficient reports log(average count)
	and non-intercept coefficients report estimated (and bias
	reduced, when appropriate) log fold-changes.)}
	                                %  Windows with at
					%  least one zero count and
					%  initial absolute coefficient
					%  estimates greater than three
					%  undergo the bias correction
					%  of Kosmidis & Firth (2009) to
					%  moderate extreme coefficient
					%  estimates.
	\item{LRT.constrained}{same as LRT, but uses the constrained MLE
	  in the full model. Each column contains statistics from a
	single hypothesis test that the ChIP coefficient is equal to
	zero versus greater than zero, applied separately to each window.}
	\item{phi.hat.dev.constrained}{same as phi.hat.dev, but subject to
	constraint that the ChIP coefficient is nonnegative.}
	\item{phi.hat.pearson.constrained}{same as phi.hat.pearson, but subject to
	constraint that the ChIP coefficient is nonnegative.}
	\item{fitted.values.constrained}{same as fitted.values, but subject to
	constraint that the ChIP coefficient is nonnegative.}
	\item{coefficients.constrained}{same as coefficients, but subject to
	constraint that the ChIP coefficient is nonnegative.}
}

\references{
Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear models" \emph{Biometrika}, \bold{96}, 793--804.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" \emph{SAGMB}, \bold{11}(5).

McCarthy, Chen and Smyth (2012) "Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation" \emph{Nucleic Acids Res.} \bold{40}(10), 4288--97.
}
\examples{
\dontrun{
set.seed(5)
####################################################
# Simulate data three replicates with one chromosome
####################################################
reps <- 3
chr.length <- 1e5
window.width <- 200
regions <- chr.length / window.width
start <- seq(1, chr.length, by = window.width)
end <- start + window.width - 1
n.peaks <- 100 # No. of true peaks
peak.idx <- sample.int(regions, n.peaks)
samples <- c(paste0('C', 1:reps), paste0('I', 1:reps))
# Set parameters
beta0 <- runif(regions, log(10), log(100))
beta1 <- rep(0, regions); beta1[peak.idx] <- VGAM::rpareto(n.peaks, shape = 5, scale = log(5))
# Set means
mu.ChIP <- exp(beta0 + beta1)
mu.input <- exp(beta0)
# Negative binomial dispersion parameter
phi <- geoR::rinvchisq(regions, df = 5)
# Draw read counts using a negative binomial distribution
C <- lapply(1:reps, function(r) MASS::rnegbin(regions,
                                              mu = mu.ChIP,
                                              theta = 1/phi))
I <- lapply(1:reps, function(r) MASS::rnegbin(regions,
                                              mu = mu.input,
                                              theta = 1/phi))
counts <- do.call('cbind', append(C, I))
colnames(counts) <- samples
rownames(counts) <- start


####################################################
# Fit quasi-negative binomial model to each window.
####################################################
design.matrix  <- cbind(rep(1, reps*2), # Intercept
                        rep(c(1,0), each = reps)) # Indicates ChIP sample
chip.col.indicator <- c(0,1) # Second column of design matrix indicates ChIP sample
fit <- QL.fit(counts, design.matrix, chip.col.indicator, Model = 'NegBin') # Fit model at window level.
}}
\author{Steve Lund (\email{lundsp@gmail.com}), Emily Goren (\email{emily.goren@gmail.com})}

     \keyword{ChIP-seq, quasi-likelihood}

