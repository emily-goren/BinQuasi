
#'
#'Fit a negative binomial GLM for a given design matrix
#'
#'@description A function called within \code{\link{QL.fit}} to fit a negative
#'  binomial GLM to each window for a given design matrix.
#'  
#'@param counts A matrix of integers obtained by counting reads across a 
#'  genomic partition. Each row contains observations from a single window of
#'  the genomic partition. Each column contains observations from a single
#'  sample (either ChIP or control/input).
#'@param design A design matrix for the full model, including a column that 
#'  indicates whether the observation is a ChIP sample (1) or control/input (0).
#'  The number of rows must be \code{ncol(counts)}. Means are modeled with a log
#'  link function.
#'@param chip.col.indicator A binary vector of length \code{ncol(design.matrix)}
#'  that indicates which column of the full design matrix corresponds to the 
#'  ChIP indicator.
#'@param log.offset A vector of log-scale, additive factors used to adjust 
#'  estimated log-scale means for differences in library sizes across samples. 
#'  Commonly used offsets include \code{log.offset=log(colSums(counts))} and 
#'  \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}. If
#'  \code{NULL}, the later offset is used.
#'@param nb.disp Estimated negative binomial dispersion parameters obtained from
#'  either \code{\link{estimateGLMTrendedDisp}} or \code{\link{estimateGLMCommonDisp}} in 
#'  package \code{\link{edgeR}}.  These estimates are treated as known and are
#'  used to compute deviances.
#'@param print.progress logical. If TRUE, the function will provide an update on
#'  what window (row number) is being analyzed. Updates occur frequently to 
#'  start then eventually occur every 5000 windows.
#'@param bias.fold.tolerance A numerical value no smaller than 1. If the bias 
#'  reduction of maximum likelihood estimates of (log) fold change is likely to 
#'  result in a ratio of fold changes greater than this value, then bias 
#'  reduction will be performed on such windows. Setting 
#'  \code{bias.fold.tolerance=Inf} will completely disable bias reduction; 
#'  setting \code{bias.fold.tolerance=1} will always perform bias reduction (see
#'  details). If the constrained estimate differs from the unconstrained 
#'  estimate, bias reduction is not performed.
#'  
#'@return A list containing: 
#'  \item{dev}{Vector containing the deviance for each
#'  window under a negative binomial model fit to design matrix specified by
#'  \code{design}. This vector is passed along within the \code{\link{QL.fit}}
#'  function.} 
#'  \item{means}{Matrix of fitted mean values for each window.} 
#'  \item{parms}{Matrix of estimated coefficients for each window. Note that
#'  these are given on the log scale. (i.e., intercept coefficient reports
#'  log(average count) and non-intercept coefficients report estimated (and bias
#'  reduced, when appropriate) log fold-changes.)} 
#'  \item{dev.constrained}{Same as \code{dev}, subject to the
#'  constraint that the ChIP coefficient is non-negative. If fitting a reduced
#'  model with no ChIP coefficient, this will be \code{NA}.} 
#'  \item{means.constrained}{Same as \code{means}, subject to
#'  the constraint that the ChIP coefficient is non-negative. If fitting a
#'  reduced model with no ChIP coefficient, this will be \code{NA}.} 
#'  \item{parms.constrained}{Same as \code{parms},
#'  subject to the constraint that the ChIP coefficient is non-negative. If
#'  fitting a reduced model with no ChIP coefficient, this will be \code{NA}.}
#'  
#'@details This functions fits, for each row of \code{counts}, a negative 
#'  binomial log-linear model through the GLM framework with the over-dispersion
#'  parameter fixed.
#'  
#'  Asymptotic biases of regression coefficients (i.e., log fold changes) are 
#'  then estimated by a plug-in estimate [eqn. (15.4) of McCullagh and Nelder, 
#'  1989] from the last iteration of iteratively reweighted least squares (IWLS)
#'  procedure. The fitted response values are then compared with or without such
#'  a bias term. If the ratio of fitted response values are larger than 
#'  \code{bias.fold.tolerance} for any observation and the unconstrained
#'  estimate equals the consrained estimate, then the bias-reduction (not 
#'  bias-correction) procedure according to Firth (1993) and Kosmidis & Firth 
#'  (2009) is applied to such rows of \code{counts}, by adjusting the score 
#'  equation with a term based on the observed information. Such bias-reduced
#'  estimates are more advantageous than directly subtracting the estimated bias
#'  from the maximum likelihood estimates as the latter may not exist (e.g.,
#'  when all counts in the control/input group are zero).
#'
#'@references Firth (1993) "Bias reduction of maximum likelihood estimates"
#'\emph{Biometrika}, \bold{80}, 27--38.
#'
#'Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear
#'models" \emph{Biometrika}, \bold{96}, 793--804.
#'
#'Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression
#'in RNA-sequence data using quasi-likelihood with shrunken dispersion
#'estimates" emph{SAGMB}, \bold{11}(5).
#'
#'McCullagh and Nelder (1989) "Generalized Linear Models", 2nd edition. London:
#'Chapman and Hall.
#'
#'@author Emily Goren (\email{emily.goren@gmail.com}) based on modifications of
#'  code by Steve Lund and Long Qu.
#'  
#'@export
#'

NBDev <- function(counts, design, log.offset, nb.disp, print.progress=TRUE, bias.fold.tolerance=1.10, chip.col.indicator) {
  
  n <- ncol(counts)
  
  if (is.null(log.offset))
    log.offset <- rep(0, ncol(counts))
  est.offset <- exp(log.offset)
  
  deviance <- deviance.constrained <- rep(NA, nrow(counts))
  means <- means.constrained <- matrix(NA, nrow(counts), ncol(counts))
  parms <- parms.constrained <- matrix(NA, nrow(counts), ncol(design))
  
  design.df=as.data.frame(design)
  glm.ctrl=glm.control(epsilon = 1e-08, maxit = 1500L, trace = FALSE)
  fbrNBglm.ctrl=fbrNBglm.control(coefOnly=TRUE, infoParms=list(j=1,k=1,m=1), maxit=1500L, tol=1e-8, standardizeX=TRUE)
  nbLogFamily=negbin("log", 1)
  logFcCorrection=abs(log(bias.fold.tolerance))
  log10=log(10)
  
  ### For each gene and given design matrix, fit glm using provided negative binomial dispersion estimate
  for (gn in 1:nrow(counts)) {
    ### If wanted, provide running progress update (eventually once every 5000 genes)
    if (gn %in% c(2, 10, 100, 500, 1000, 2500, 5000 * (1:200)) & print.progress)
      print(paste("Analyzing Window #", gn))
    
    #### For 2000 Fly genes, glm takes roughly 9 seconds (optim took roughly 21 seconds)
    res <- glmsolve(formula = counts[gn, ] ~ . - 1 + offset(log.offset), family = update.fbrNBfamily(nbLogFamily, overdisp=nb.disp[gn]), data = design.df, control = glm.ctrl, x=TRUE)
    
    if (ncol(design) == length(chip.col.indicator)) {
      if (res$coefficients[which(chip.col.indicator == 1)] < 0) {
        Amat = t(chip.col.indicator)
        res.constrained  <- orglm.fit(x = design.df, y = counts[gn, ], offset = log.offset,
                                      family = update.fbrNBfamily(nbLogFamily, overdisp=nb.disp[gn]),
                                      nec = 0, constr = Amat, rhs = 0, control = glm.ctrl)
      } else {
        res.constrained <- res
      }
      parms.constrained[gn, ] <- res.constrained$coefficients
      deviance.constrained[gn] <- res.constrained$deviance
      means.constrained[gn, ] <- as.vector(exp(design[,,drop=FALSE] %*% parms.constrained[gn, ]) * est.offset)
    }
    
    parms[gn, ] <- res$coefficients
    
    ### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
    deviance[gn] <- res$deviance
    
    
    ### When a group has only zero counts, the MLE doesn't exist.  For these genes, we use the method Kosmidis & Firth(2009)
    ### to moderate the parameter estimates.
    ### For 2000 Fly genes, fbr takes roughly 41 seconds
    #if (any(counts[gn, ] == 0) && any(abs(coefficients(glm.fitted)) > 3)) {
    glm.fitted = res
    tmp.bias = coef(glm.fitted, type='bias')
    tmp.nonNA=which(!is.na(tmp.bias) & !is.na(glm.fitted$coefficients))
    this.x=model.matrix(glm.fitted)
    if(any(abs(this.x[,tmp.nonNA,drop=FALSE]%*%tmp.bias[tmp.nonNA]) > logFcCorrection)) {
      fbrNBglm.ctrl$start = glm.fitted$coefficients - pmax.int(-log10, pmin.int(log10, tmp.bias))
      wt.idx=which(glm.fitted$weights>0)
      fbrglm.fit <- fbrNBglm.fit(x=this.x[wt.idx,,drop=FALSE], y=counts[gn, wt.idx],offset=log.offset[wt.idx], family=nbLogFamily, control = fbrNBglm.ctrl)
      ## note: As fbrNBglm.fit is only called right after a update.fbrfamily call, the overdisp parameter is already up-to-date
      parms[gn, ] <- fbrglm.fit
      if (ncol(design) == length(chip.col.indicator)) {
        if (!(parms[gn, which(chip.col.indicator == 1)] < 0)) {
          parms.constrained[gn, ] <- parms[gn, ]
        }
      }
      tmp.nonNA = which(!is.na(fbrglm.fit))
    }
    ### Save optimized means (used in Pearson's dispersion estimator)
    ### NOTE: This line was originally before the Firth bias correction block.
    means[gn, ] <- as.vector(exp(design[,,drop=FALSE] %*% parms[gn, ]) * est.offset)
    
  }
  return(list(dev=deviance, means=means, parms=parms,
              dev.constrained=deviance.constrained, means.constrained=means.constrained, parms.constrained=parms.constrained))
}



