% File src/library/QuasiSeq/man/QL.results.Rd
     \name{QL.results}
     \alias{QL.results}
     \title{Obtain p-values and q-values using results from \code{QL.fit}}
     \description{
       Obtain significance results, at the candidate region level, for
     quasi-likelihood model fits to ChIP-seq data partitioned into
     counts using the methods detailed in Goren, Wang, and
     Liu (to be submitted).}
   \usage{
QL.results(fit, Dispersion="Deviance", one.sided=TRUE,
           spline.df=NULL, Plot=FALSE, padj=TRUE)
	     }
     \arguments{
      \item{fit}{The list returned by the function QL.fit}
	\item{Dispersion}{Must be one of "Deviance" or "Pearson", specifying which type of estimator should be used for estimating quasi-likelihood dispersion parameter.}
      \item{one.sided}{logical. If TRUE, a one-sided test for the ChIP
     coefficient is reported. Otherwise, if FALSE, a two-sided test is reported.}
	\item{spline.df}{Optional.  User may specify the degrees of freedom to use when fitting a cubic spline to log-scale(estimated dispersion) versus the log(average count).  Default uses cross-validation in \code{sreg} function to pick appropriate degrees of freedom.}
	\item{Plot}{logical.  If TRUE, the estimated dispersion versus
     the average count are plotted on a log-scale with the corresponding
     cubic spline fit overlaid. }
   \item{padj}{logical. If TRUE, Benjamini & Hochberg's adjustment for multiple comparisons is applied.}
	}

\value{list containing:
	\item{"P.values"}{list of matrices providing p-values for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{pvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.}
	\item{"Q.values"}{list of matrices providing q-values for the QL, QLShrink and QLSpline methods, respectively.  The i^th column of each element of \code{qvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.
 Q-values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
	\item{"F.stat"}{list of matrices providing F-statistics for the
	QL, QLShrink and QLSpline methods, respectively.  The i^th
	column of each element of \code{F.stat} corresponds to the
	hypothesis test assigned in the i^th row of \code{test.mat}.}
	\item{"d0"}{ vector containing estimated additional denominator degrees of freedom  gained from shrinking dispersion estimates in the QLShrink and QLSpline procedures, respectively.}
}

\author{Steve Lund \email{lundsp@gmail.com}, Emily Goren \email{emily.goren@gmail.com}}

\references{
Benjamini and Hochberg (1995)
"Controlling the false discovery rate: a practical and powerful approach to multiple testing" \emph{Journal of the Royal Statistical Society Series B}, \bold{57}: 289-300.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" \emph{SAGMB}, \bold{11}(5).
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
window.results <- QL.results(fit)

# Number of significant windows.
sum(window.results$Q.values$QLShrink < 0.05)
}
}

     \keyword{ChIP-seq, quasi-likelihood}



