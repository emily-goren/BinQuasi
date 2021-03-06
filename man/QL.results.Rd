% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QL.results.R
\name{QL.results}
\alias{QL.results}
\title{Obtain p- and q-values using results from \code{QL.fit}}
\usage{
QL.results(fit, Dispersion = "Deviance", one.sided = TRUE,
  spline.df = NULL, Plot = FALSE, padj = TRUE)
}
\arguments{
\item{fit}{The list returned by the function \code{QL.fit}.}

\item{Dispersion}{Must be one of \code{"Deviance"} or \code{"Pearson"}, specifying which
type of estimator should be used for estimating the quasi-likelihood
dispersion parameter.}

\item{one.sided}{logical. If TRUE, a one-sided test for the ChIP coefficient
is reported. Otherwise, if FALSE, a two-sided test is reported.}

\item{spline.df}{Optional. User may specify the degrees of freedom to use when
fitting a cubic spline to log-scale(estimated dispersion) versus the
log(average count). Default uses cross-validation in \code{sreg} function to
pick appropriate degrees of freedom.}

\item{Plot}{logical. If TRUE, the estimated dispersion versus the average
count are plotted on a log-scale with the corresponding cubic spline fit
overlaid.}

\item{padj}{logical. If TRUE, Benjamini & Hochberg's adjustment for multiple
comparisons is applied.}
}
\value{
A list containing: 
 \item{P.values}{List of matrices providing
 p-values for the QL, QLShrink and QLSpline methods, respectively. The i-th
 column of each element of \code{pvals} corresponds to the hypothesis test
 assigned in the i-th window.} 
 \item{Q.values}{List of
 matrices providing q-values for the QL, QLShrink and QLSpline methods,
 respectively. The i-th column of each element of \code{qvals} corresponds to
 the hypothesis test assigned in the i-th window.} 
 \item{F.stat}{List of matrices providing F-statistics for the QL, QLShrink
 and QLSpline methods, respectively. The i-th column of each element of
 \code{F.stat} corresponds to the hypothesis test assigned in the i-th window.} 
 \item{d0}{Vector containing estimated additional
 denominator degrees of freedom  gained from shrinking dispersion estimates
 in the QLShrink and QLSpline procedures, respectively.}
}
\description{
Obtain significance results for quasi-likelihood models fit to
 ChIP-seq data partitioned into counts.
}
\details{
Obtain significance results from an object fitted using \code{\link{QL.fit}}.
Used within the main peak calling function, \code{\link{BQ}}.
}
\examples{
 set.seed(5)
 ####################################################
 # Simulate data three replicates with one chromosome
 ####################################################
 reps <- 3
 chr.length <- 1e5
 window.width <- 200
 K <- chr.length / window.width
 start <- seq(1, chr.length, by = window.width)
 end <- start + window.width - 1
 n.peaks <- 100 # No. of true peaks
 peak.idx <- sample.int(K, n.peaks)
 samples <- c(paste0('C', 1:reps), paste0('I', 1:reps))
 # Set parameters
 beta0 <- runif(K, log(10), log(100))
 beta1 <- rep(0, K); beta1[peak.idx] <- log(5) / runif(n.peaks)^(1/5)
 # Set means
 mu.ChIP <- exp(beta0 + beta1)
 mu.input <- exp(beta0)
 # Negative binomial dispersion parameter
 phi <- 1/rchisq(K, df = 5)
 # Draw read counts using a negative binomial distribution
 C <- lapply(1:reps, function(r) rpois(K, (mu.ChIP * rgamma(K, 1/phi))/(1/phi)))
 I <- lapply(1:reps, function(r) rpois(K, (mu.input * rgamma(K, 1/phi))/(1/phi)))
 counts <- do.call('cbind', append(C, I))
 colnames(counts) <- samples
 rownames(counts) <- start
 
 ####################################################
 # Fit quasi-negative binomial model to each window.
 ####################################################
 design.matrix  <- cbind(rep(1, reps*2), # Intercept
                         rep(c(1,0), each = reps)) # Indicates ChIP sample
 chip.col.indicator <- c(0,1) # Second column of design matrix indicates ChIP sample
 fit <- QL.fit(counts, design.matrix, chip.col.indicator, 
               log.offset = rep(1, ncol(counts)), Model = 'NegBin') 
 window.results <- QL.results(fit)

# Number of significant windows.
sum(window.results$Q.values$QLShrink < 0.05)

# Compare significant windows to truth.
res <- as.numeric(window.results$Q.values$QLShrink < 0.05)
# Number of true positives
TP <- sum(res[peak.idx] == 1)
TP
# Number of false negatives
FN <- n.peaks - TP
FN
# Number of false positives
FP <- sum(res) - TP
FP
# Number of true negatives
TN <- (K - sum(res)) - FN
TN

}
\references{
Goren, Liu, Wang and Wang (2018) "BinQuasi: a peak detection method for 
ChIP-sequencing data with biological replicates" \emph{Bioinformatics}.

Benjamini and Hochberg (1995) "Controlling the false discovery
rate: a practical and powerful approach to multiple testing" \emph{Journal of
the Royal Statistical Society Series B}, \bold{57}: 289-300.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression
in RNA-sequence data using quasi-likelihood with shrunken dispersion
estimates" \emph{SAGMB}, \bold{11}(5).
}
\author{
Emily Goren (\email{emily.goren@gmail.com}) based on modifications of
 code by Steve Lund.
}
