#' 
#' Fit quasi-likelihood models to replicated ChIP-seq data partitioned into a
#' count matrix
#' 
#' @description Fit constrained quasi-likelihood models to ChIP-seq data 
#'   partitioned into a count matrix.
#'   
#' @details A wrapper for \code{\link{PoisDev}} or \code{\link{NBDev}}, 
#'   depending on whether quasi-Poisson or quasi-negative binomial models are 
#'   requested. See the respective functions for details. Used within the main
#'   \code{\link{BQ}} peak calling function.
#'   
#' @param counts A matrix of integers obtained by counting reads across a 
#'   genomic partition. Each row contains observations from a single window of 
#'   the genomic partition. Each column contains observations from a single
#'   sample (either ChIP or control/input).
#' @param design.matrix A design matrix for the full model, including a column 
#'   that indicates whether the observation is a ChIP sample. The number of rows
#'   must be \code{ncol(counts)}.  The number of columns must be at least two, 
#'   usually an intercept and an indicator whether the sample is ChIP (1) or 
#'   input/control (0). Means are modeled with a log link function.
#' @param chip.col.indicator A binary vector of length 
#'   \code{ncol(design.matrix)} that indicates which column of
#'   \code{design.matrix} corresponds to the ChIP indicator.
#' @param log.offset A vector of log-scale, additive factors used to adjust
#'   estimated log-scale means for differences in library sizes across samples. 
#'   Commonly used offsets include \code{log.offset=log(colSums(counts))} and 
#'   \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}. 
#'   If \code{NULL}, the later offset is used.
#' @param Model Must be one of \code{"Poisson"} or \code{"NegBin"}, specifying
#'   use of a quasi-Poisson or quasi-negative binomial model, respectively.
#' @param print.progress logical. If \code{TRUE}, updates are provided regarding
#'   which window (row number) is being analyzed.  Updates occur frequently to
#'   start then eventually occur every 5000 windows.
#' @param NBdisp Used only when \code{Model="NegBin"}. Must be one of \code{"trend"},
#'   \code{"common"} or a vector of non-negative real numbers with length equal to
#'   \code{nrow(counts)}. Specifying \code{NBdisp="trend"} or 
#'   \code{NBdisp="common"} will use {\link{estimateGLMTrendedDisp}} or
#'   \code{\link{estimateGLMCommonDisp}}, respectively, from the package
#'   \code{\link{edgeR}} to estimate negative binomial dispersion parameters for each
#'   window. Estimates obtained from other sources can be used by entering
#'   \code{NBdisp} as a vector containing the negative binomial dispersion value
#'   to use for each window when fitting the quasi-likelihood model.
#' @param bias.fold.tolerance A numerical value no smaller than 1. If the bias
#'   reduction of maximum likelihood estimates of (log) fold change is likely to
#'   result in a ratio of fold changes greater than this value, then bias
#'   reduction will be performed on such windows. Setting
#'   \code{bias.fold.tolerance=Inf} will completely disable bias reduction; 
#'   setting \code{bias.fold.tolerance=1} will always perform bias reduction. 
#'   Estimates that are projected into the constrained space are not bias-reduced.
#' @param ... Arguments to be passed to the function
#'   \code{\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}}
#'   from the package \code{\link{edgeR}}.
#'   
#' @return  A list containing: 
#'   \item{LRT}{Matrix providing unadjusted two-sided
#'   likelihood ratio test statistics.  Each column contains statistics from a 
#'   single hypothesis test that the ChIP coefficient is equal to zero versus
#'   not equal to zero, applied separately to each window.} 
#'   \item{phi.hat.dev}{Vector providing unshrunken, deviance-based estimates of
#'   the quasi-dispersion parameter for each window.} 
#'   \item{phi.hat.pearson}{Vector providing unshrunken, Pearson estimates of 
#'   the quasi-dispersion parameter for each window.} 
#'   \item{mn.cnt}{Vector of the average count for each window.} 
#'   \item{den.df}{Denominator degrees of
#'   freedom. Equal to the number of samples minus the number of fitted
#'   parameters in the full model.} 
#'   \item{num.df}{Vector of numerator degrees of
#'   freedom for each test, computed as the difference in the number of fitted
#'   parameters between the full and reduced models.} 
#'   \item{Model}{Either
#'   "Poisson" or "NegBin", specifying which model (quasi-Poisson or
#'   quasi-negative binomial, respectively) was used.} 
#'   \item{nb.disp}{Only
#'   appears when \code{Model="NegBin"}. Vector providing negative binomial
#'   dispersion parameter estimates used during fitting of quasi-negative
#'   binomial model for each window.} 
#'   \item{fitted.values}{Matrix of fitted mean
#'   values without constraints.} 
#'   \item{coefficients}{Matrix of estimated
#'   coefficients for each window. Note that these are given on the log scale.
#'   (i.e., intercept coefficients report log(average count) and non-intercept
#'   coefficients report estimated (and bias reduced, when appropriate) log
#'   fold-changes.)} 
#'   \item{LRT.constrained}{Same as \code{LRT}, but uses the
#'   constrained MLE in the full model. Each column contains statistics from a
#'   single hypothesis test that the ChIP coefficient is equal to zero versus
#'   greater than zero, applied separately to each window.} 
#'   \item{phi.hat.dev.constrained}{Same as \code{phi.hat.dev}, but subject to the
#'   constraint that the ChIP coefficient is non-negative.} 
#'   \item{phi.hat.pearson.constrained}{Same as \code{phi.hat.pearson}, but subject to
#'   the constraint that the ChIP coefficient is non-negative.} 
#'   \item{fitted.values.constrained}{Same as \code{fitted.values}, but subject to the 
#'   constraint that the ChIP coefficient is non-negative.} 
#'   \item{coefficients.constrained}{Same as \code{coefficients}, but subject to the
#'   constraint that the ChIP coefficient is non-negative.}
#'   
#' @references
#' 
#' Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear
#' models" \emph{Biometrika}, \bold{96}, 793--804.
#' 
#' Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression
#' in RNA-sequence data using quasi-likelihood with shrunken dispersion
#' estimates" \emph{SAGMB}, \bold{11}(5).
#' 
#' McCarthy, Chen and Smyth (2012) "Differential expression analysis of
#' multifactor RNA-Seq experiments with respect to biological variation"
#' \emph{Nucleic Acids Res.} \bold{40}(10), 4288--97.
#' 
#' 
#' @examples
#'  set.seed(5)
#'  ####################################################
#'  # Simulate data three replicates with one chromosome
#'  ####################################################
#'  reps <- 3
#'  chr.length <- 1e5
#'  window.width <- 200
#'  K <- chr.length / window.width
#'  start <- seq(1, chr.length, by = window.width)
#'  end <- start + window.width - 1
#'  n.peaks <- 100 # No. of true peaks
#'  peak.idx <- sample.int(K, n.peaks)
#'  samples <- c(paste0('C', 1:reps), paste0('I', 1:reps))
#'  # Set parameters
#'  beta0 <- runif(K, log(10), log(100))
#'  beta1 <- rep(0, K); beta1[peak.idx] <- log(5) / runif(n.peaks)^(1/5)
#'  # Set means
#'  mu.ChIP <- exp(beta0 + beta1)
#'  mu.input <- exp(beta0)
#'  # Negative binomial dispersion parameter
#'  phi <- 1/rchisq(K, df = 5)
#'  # Draw read counts using a negative binomial distribution
#'  C <- lapply(1:reps, function(r) rpois(K, (mu.ChIP * rgamma(K, 1/phi))/(1/phi)))
#'  I <- lapply(1:reps, function(r) rpois(K, (mu.input * rgamma(K, 1/phi))/(1/phi)))
#'  counts <- do.call('cbind', append(C, I))
#'  colnames(counts) <- samples
#'  rownames(counts) <- start
#'  head(counts)
#'  
#'  ####################################################
#'  # Fit quasi-negative binomial model to each window.
#'  ####################################################
#'  design.matrix  <- cbind(rep(1, reps*2), # Intercept
#'                          rep(c(1,0), each = reps)) # Indicates ChIP sample
#'  chip.col.indicator <- c(0,1) # Second column of design matrix indicates ChIP sample
#'  fit <- QL.fit(counts, design.matrix, chip.col.indicator, 
#'                log.offset = rep(1, ncol(counts)), Model = 'NegBin') 
#'  # Look at fitted values
#'  counts.fitted <- fit$fitted.values.constrained
#'  head(round(counts.fitted, 2))
#'  
#' @author Emily Goren (\email{emily.goren@gmail.com}) based on modifications of
#'   code by Steve Lund.
#'   
#' @export
#' 

QL.fit <- function(counts, design.matrix, chip.col.indicator, log.offset = NULL, Model = "NegBin", print.progress = TRUE,
                   NBdisp = "trend", bias.fold.tolerance=1.10, ...) {
  if(is.data.frame(counts)) counts<-as.matrix(counts)
  ### Note: First element of design.list should pertain to overall full model.  This is the design used to obtain
  ### dispersion estimates for quasi-likelihood models.
  if(any(rowSums(counts)==0)) stop(paste(sum(rowSums(counts)==0)," genes have 0 counts across all samples.  Please remove genes with zero total counts before analyzing.",sep=""))
  ### Check for errors
  if (any(round(counts) != counts))
    stop("Count data contains non-integers.")
  if (any(counts < 0))
    stop("Count data contains negative counts.")
  
  if (!Model %in% c("NegBin", "Poisson"))
    stop("Unidentified Model: Model must be either 'NegBin' or 'Poisson'.")
  if (Model == "NegBin") {
    
    if (!NBdisp[1] %in% c("trend", "common") & length(NBdisp) != nrow(counts))
      stop("Unidentified NegBin Dispersion: NBdisp must be set as 'trend' or 'common' to estimate negative binomial dispersion from data using GLM edgeR (McCarthy et al., 2012),\n or it must be a vector providing negative binomial dispersion parameter value to use for each window.")
    
    if (length(NBdisp) == nrow(counts) & !is.numeric(NBdisp))
      stop("NBdisp contains non-numeric values.\n\tAll negative binomial dispersion parameters must be non-negative real numbers.")
    
    if (length(NBdisp) == nrow(counts) & any(NBdisp < 0))
      stop("NBdisp contains negative values.\nAll negative binomial dispersion parameters must be non-negative real numbers.")
  }
  ### Default log offset to 0.75 quantile.
  if (is.null(log.offset)) {
    size <- apply(counts, 2, quantile, 0.75)
    log.offset <- log(size)
  }
  
  ### Only one inequality constraint, for the ChIP coefficient, is supported.
  if (is.null(chip.col.indicator))
    stop("Vector indicating which column of the full design matrix corresponds to the ChIP coefficient was not provided.")
  if (length(chip.col.indicator) != ncol(design.matrix))
    stop("Vector indicating which column of the full design matrix corresponds to the ChIP coefficient differs in length from number of design matrix rows.")
  if (sum(!(chip.col.indicator %in% c(0,1))) != 0)
    stop("Vector indicating which column of the full design matrix corresponds to the ChIP coefficient is not binary.")
  if (sum(chip.col.indicator) != 1)
    stop("Vector indicating which column of the full design matrix corresponds to the ChIP coefficient indicates more than one row.")
  
  ## Currently only supports a full model versus a reduced model, which lacks the chip coefficient.
  if (!is.matrix(design.matrix))
    stop("Design matrix provided is not a matrix. Please provide a matrix.")
  if (nrow(design.matrix) != ncol(counts))
    stop('The number of rows of the design matrix does not match that of the count table.')
  design.list <- rep(list(NA),2)
  design.list[[1]] <- design.matrix
  design.list[[2]] <- design.matrix[ , !chip.col.indicator]
  
  
  ### Fit model and evaluate deviance under each design provided in design.list
  deviance.list <- deviance.list.constrained <- vector("list", length(design.list))
  p <- NULL
  n <- ncol(counts)  # p is used to store the d.f. for each model (it will be a vector) and n is the total number of samples (also the number of observations for each gene)
  
  for (jj in 1:length(design.list)) {
    design <- design.list[[jj]]
    
    if (is.vector(design)) {
      p <- c(p, length(unique(design)))  ## Record the d.f. for current model
      
      ### Check for errors if the current model design is specified as a vector
      if (p[jj] > p[1])
        stop(paste("Full model design must be first element in 'design.list'.\n'p' for element", jj, "is larger than 'p' for first element,\nindicating first element does not provide full model design."))
      if (length(design) != n)
        stop(paste("Element", jj, "in 'design.list' has length", length(design), ".\nDesign vectors must have length",
                   n, "(to match number of columns in data)."))
    }
    
    if (is.matrix(design)) {
      p <- c(p, ncol(design))
      
      ### Check for errors if the current model design is specified as a matrix
      ### if(prod(design[,1]==1)!=1)
      ### stop(paste('The first column of matrix in element',jj,'of 'design.list' is not a column of 1s for the
      ### intercept. Please include intercept.'))
      if (nrow(design) != n)
        stop(paste("Element", jj, "in 'design.list' has", nrow(design), "rows.\nDesign matrices must have",
                   n, "rows (to match number of columns in data)."))
      if (p[jj] > p[1])
        stop(paste("Full model design must be first element in 'design.list'.\n'p' for element", jj, "is larger than 'p' for first element,\nindicating first element does not provide full model design."))
    }
    
    ### Analyze using quasi-negative binomial model, if chosen
    if (Model == "NegBin") {
      if (is.vector(design)) {
        if (length(unique(design)) > 1)
          design <- model.matrix(~as.factor(design))
        if (length(unique(design)) == 1)
          design <- matrix(1, ncol(counts), 1)
      }
      if (jj == 1) {
        ### Obtain negative binomial dispersion estimates from edgeR, if requested
        if(NBdisp %in% c("trend", "common")){
          
          
          ### If provided, use prespecified library size normalization
          if (!is.null(log.offset))
            d <- DGEList(counts = counts, group = design[, 2], lib.size = exp(log.offset))
          
          ### If requested, use gene-specific trended dispersion estimates from GLM edgeR (McCarthy et al., 2012).
          if (NBdisp == "trend")
            nb.disp <- estimateGLMTrendedDisp(d, design)$trended.dispersion
          
          ### If requested, use common dispersion estimate from GLM edgeR (McCarthy et al., 2012).
          if (NBdisp == "common")
            nb.disp <- rep(estimateGLMCommonDisp(d, design, ...)$common.dispersion, nrow(counts))
        }
        ### If provided, use prespecified dispersion estimates.
        if (length(NBdisp) == nrow(counts)) {
          if (is.numeric(NBdisp) & !any(NBdisp < 0))
            nb.disp <- NBdisp
        }
      }
      
      ### Analyze genes with positive dispersion parameters using quasi-negative binomial model
      if (any(nb.disp > 0))
        res <- NBDev(counts[nb.disp > 0, ], design, log.offset, nb.disp[nb.disp > 0], print.progress, bias.fold.tolerance=bias.fold.tolerance, chip.col.indicator=chip.col.indicator)
      
      ### If present, analyze genes for which nb.disp==0 using quasi-Poisson model
      if (any(nb.disp == 0)) {
        res2 <- PoisDev(counts[nb.disp == 0, ], design, log.offset, print.progress,  bias.fold.tolerance=bias.fold.tolerance, chip.col.indicator=chip.col.indicator)
        means <- dev <- means.constrained <- dev.constrained <- rep(NA, nrow(counts))
        parms <- parms.constrained <- matrix(NA, nrow(counts), ncol(design))
        means[nb.disp == 0] <- res2$means
        dev[nb.disp == 0] <- res2$dev
        parms[nb.disp == 0, ] <- res2$parms
        means.constrained[nb.disp == 0] <- res2$means.constrained
        dev.constrained[nb.disp == 0] <- res2$dev.constrained
        parms.constrained[nb.disp == 0, ] <- res2$parms.constrained
        if (any(nb.disp > 0)) {
          means[nb.disp > 0] <- res$means
          dev[nb.disp > 0] <- res$dev
          parms[nb.disp > 0, ] <- res$parms
          means.constrained[nb.disp > 0] <- res$means.constrained
          dev.constrained[nb.disp > 0] <- res$dev.constrained
          parms.constrained[nb.disp > 0, ] <- res$parms.constrained
        }
        res <- list(dev = dev, means = means, parms = parms,
                    dev.constrained = dev.constrained, means.constrained = means.constrained, parms.constrained = parms.constrained)      }
    }
    
    ### Analyze using quasi-Poisson model, if chosen
    if (Model == "Poisson") res <- PoisDev(counts, design, log.offset, print.progress,  bias.fold.tolerance=bias.fold.tolerance, chip.col.indicator=chip.col.indicator)
    ### Record means and parameter estimate from full model
    if (jj == 1) {
      means <- res$means
      parms <- res$parms
      means.constrained <- res$means.constrained
      parms.constrained <- res$parms.constrained
    }
    deviance.list[[jj]] <- res$dev
    deviance.list.constrained[[jj]] <- res$dev.constrained
    
    if (jj == 2) {
      means.null <- res$means
    }
    deviance.list[[jj]] <- res$dev
    deviance.list.constrained[[jj]] <- res$dev.constrained
  }
  
  LRT <- LRT.constrained <- num.df <- NULL
  
  if (length(design.list) == 2) {
    ### Compute likelihood ratio test statistics. Compare reduced model to the full model
    
    test.mat <- cbind(1, 2)
    rownames(test.mat) <- "Full vs Reduced"
    
    for (i in 1:nrow(test.mat)) {
      i1 <- test.mat[i, 1]
      i2 <- test.mat[i, 2]
      num.df <- c(num.df, abs(p[i2] - p[i1]))
      LRT <- cbind(LRT, -(deviance.list[[i2]] - deviance.list[[i1]])/(p[i2] - p[i1]))
      LRT.constrained <- cbind(LRT.constrained, -(deviance.list[[i2]] - deviance.list.constrained[[i1]])/(p[i2] - p[i1]))
    }
    colnames(LRT) <- rownames(test.mat)
  }
  den.df <- (n - p[1])
  
  ### Compute deviance dispersion estimate
  phi.hat.dev <- deviance.list[[1]]/den.df
  phi.hat.dev.constrained <- deviance.list.constrained[[1]]/den.df
  phi.hat.dev.null <- deviance.list[[2]]/den.df
  
  ### Compute Pearson dispersion estimate
  if (Model == "NegBin") {
    phi.hat.pearson <- (means - counts)^2/(means + means^2 * nb.disp)
    phi.hat.pearson.constrained <- (means.constrained - counts)^2/(means.constrained + means.constrained^2 * nb.disp)
    phi.hat.pearson.null <- (means.null - counts)^2/(means.null + means.null^2 * nb.disp)
  }
  if (Model == "Poisson") {
    phi.hat.pearson <- (means - counts)^2/means
    phi.hat.pearson.constrained <- (means.constrained - counts)^2/means.constrained
    phi.hat.pearson.null <- (means.null - counts)^2/means.null
  }
  phi.hat.pearson[means == 0] <- 0
  phi.hat.pearson <- rowSums(phi.hat.pearson)/den.df
  phi.hat.pearson.constrained[means.constrained == 0] <- 0
  phi.hat.pearson.null[means.null == 0] <- 0
  phi.hat.pearson.constrained <- rowSums(phi.hat.pearson.constrained)/den.df
  phi.hat.pearson.null <- rowSums(phi.hat.pearson.null)/den.df
  
  if (Model == "Poisson")
    return(list(LRT = LRT, LRT.constrained = LRT.constrained, phi.hat.dev = phi.hat.dev, phi.hat.pearson = phi.hat.pearson, mn.cnt = rowMeans(counts),
                den.df = den.df, num.df = num.df, Model = Model, fitted.values = means, coefficients = parms,
                phi.hat.dev.constrained = phi.hat.dev.constrained, phi.hat.pearson.constrained = phi.hat.pearson.constrained,
                fitted.values.constrained = means.constrained, coefficients.constrained = parms.constrained))
  if (Model == "NegBin")
    return(list(LRT = LRT, LRT.constrained = LRT.constrained, phi.hat.dev = phi.hat.dev, phi.hat.pearson = phi.hat.pearson, mn.cnt = rowMeans(counts),
                den.df = den.df, num.df = num.df, Model = Model, NB.disp = nb.disp, fitted.values = means, coefficients = parms,
                phi.hat.dev.constrained = phi.hat.dev.constrained, phi.hat.pearson.constrained = phi.hat.pearson.constrained,
                fitted.values.constrained = means.constrained, coefficients.constrained = parms.constrained))
  
}
