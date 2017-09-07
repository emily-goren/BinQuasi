% File src/library/QuasiSeq/man/PoisDev.Rd
     \name{PoisDev}
     \alias{PoisDev}
     \title{Compute Poisson deviances (up to a constant) for given design matrix}
     \description{A function called within \code{QL.fit} to compute Poisson deviances of each window for given design matrix}
     \usage{
PoisDev(counts,design,log.offset,print.progress=TRUE,
        bias.fold.tolerance=1.10,chip.col.indicator)
     }
    \arguments{
       \item{counts}{A matrix of integer counts obtained by counting
     reads accross a genomic partition. Each row contains observations
     from a single window. Each column contains observations from a single sample (either ChIP or control/input).}
    \item{design}{A design matrix for the full model, including a
     column that indicates whether the observation is a ChIP sample (1)
     or control/input (0). The number
     of rows must be \code{ncol(counts)}. Means are modeled with a
     log link function.}
     \item{chip.col.indicator}{A binary vector of length
     \code{ncol(design.matrix)} that indicates which column of the full
     design matrix corresponds to the ChIP indicator.}
	\item{log.offset}{A vector of log-scale, additive factors used
     to adjust estimated log-scale means for differences in library
     sizes across samples.  Commonly used offsets include
     \code{log.offset=log(colSums(counts))} and
     \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}.  If not provided, the later offset is used.}
	\item{print.progress}{logical. If TRUE, the function will
     provide an update on what window (row number) is being analyzed.
     Updates occur frequently to start then eventually occur every 5000
     genes.}
	\item{bias.fold.tolerance}{A numerical value no smaller than
     1. If the bias reduction of maximum likelihood estimates of (log)
     fold change is likely to result in a ratio of fold changes greater
     than this value, then bias reduction will be performed on such
     genes. Setting \code{bias.fold.tolerance=Inf} will completely
     disable bias reduction;  setting \code{bias.fold.tolerance=1} will
     always perform bias reduction (see details). }
}
\value{list containing:
\item{dev}{vector containing the deviance for each gene under a Poisson model fit to design matrix (or vector, for one-factor experiments) specified by \code{design}. This vector is passed along within the \code{QL.fit} function.}
  \item{means}{matrix of fitted mean values for each gene}
  \item{parms}{matrix of estimated coefficients for each gene}
 	\item{dev.constrained}{same as dev, subject to the
	  constraint that the ChIP coefficient is nonnegative. If
	  fitting a reduced model, this will be NA.}
	\item{means.constrained}{same as means, subject to the
	  constraint that the ChIP coefficient is nonnegative. If
	  fitting a reduced model, this will be NA.}
	\item{parms.constrained}{same as parms, subject to the
	constraint that the ChIP coefficient is nonnegative. If
	  fitting a reduced model, this will be NA.}
}


\details{
This functions fits, for each row of \code{counts}, a poisson log-linear model. Asymptotic biases of regression coefficients (i.e., log fold changes) are then estimated by a plug-in estimate [eqn. (15.4) of McCullagh and Nelder, 1989] from the last iteration of iteratively re-weighted least squares (IWLS) procedure. The fitted response values are then compared with or without such a bias term. If the ratio of fitted response values are larger than \code{bias.fold.tolerance} for any observation, then the bias-reduction (not bias-correction) procedure according to Firth (1993) and Kosmidis & Firth (2009) is applied to such rows of \code{counts}, by adjusting the score equation with a term based on the observed information. Such bias-reduced estimates are more advantageous than directly subtracting the estimated bias from the maximum likelihood estimates as the latter may not exist (e.g., when all counts in one treatment group are zeros).
}

\references{
Firth (1993) "Bias reduction of maximum likelihood estimates" \emph{Biometrika}, \bold{80}, 27--38.

Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear models" \emph{Biometrika}, \bold{96}, 793--804.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" emph{SAGMB}, \bold{11}(5).

McCullagh and Nelder (1989) "Generalized Linear Models", 2nd edition. London: Chapman and Hall.
}

\author{Steve Lund (\email{lundsp@gmail.com}), Emily Goren (\email{emily.goren@gmail.com})}

     \keyword{ChIP-seq, quasi-likelihood, Poisson}





