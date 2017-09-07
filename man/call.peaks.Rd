% File src/library/QuasiSeq/man/call.peaks.Rd
     \name{call.peaks}
     \alias{call.peaks}
     \title{Call peaks from a list of p-values}
     \description{
       Call peaks from p-values corresponding to a test on windows of a
       genomic partition of ChIP-seq counts. After correcting for
     multiple testing (using the Benjamini and Hochberg method), windows
     that are significant according to the threshold provided are
     retained as candidate regions after merging adjacent significatn
     windows. Simes' adjustment is used to produce region-level p- and
     q-values. Regions with q-values smaller than the significance
     threshold are returned as peaks.}
     \usage{
         call.peaks(window.pvals, start, end, chromosomes, alpha=0.05)
}
     \arguments{
      \item{window.pvals}{A vector of p-values, with each element
     corresponding to a window of a genomic partition. Typically
     obtained from the \code{QL.fit} and \code{QL.results} functions.}
	\item{start}{Vector of the genomic start locations corresponding
	  to the supplied p-values.}
	\item{end}{Vector of the genomic end locations corresponding
	  to the supplied p-values.}
	\item{chromosomes}{Vector of the chromosome names corresponding
     to the supplied p-values.}
   \item{alpha}{The desired significance threshold in (0, 0.5).}
	}

\value{the called peaks as a data.frame with variables:
	\item{"start"}{Genomic start locations of the called peaks.}
	\item{"end"}{Genomic end locations of the called peaks.}
	\item{"width"}{Width of the called peaks.}
	\item{"chr"}{Chromosomes of the called peaks.}
	\item{"P.val"}{p-values of the called peaks (adjusted using
	Simes' procedure).}
	\item{"Q.val"}{q-values of the called peaks (computing using the
	Benjamini and Hochberg procedure).}
}

\author{Emily Goren \email{emily.goren@gmail.com}}

\references{
Benjamini and Hochberg (1995)
"Controlling the false discovery rate: a practical and powerful approach to multiple testing" \emph{Journal of the Royal Statistical Society Series B}, \bold{57}: 289-300.

Simes (1986)
"An improved Bonferroni procedure for multiple tests of significance"
\emph{Biometrika}, \bold{73}(3): 751-754.

}
\examples{
# Example for a single chromosome.
start <- seq(1, 1e6, by = 200)
end <- start + 200 - 1
chromosomes <- rep('chr1', length(start))
p <- c(runif(length(start) - 10), rep(0.0001, 10))
called <- call.peaks(p, start, end, chromosomes, alpha = 0.1)
called
}

     \keyword{ChIP-seq, peak caller}



