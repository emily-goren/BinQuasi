#' 
#' Compute Poisson deviances for a given design matrix
#' 
#' @description A function called within \code{\link{QL.fit}} to compute Poisson
#' deviances of each window for a given design matrix.
#' 
#' @param counts A matrix of integer counts obtained by counting reads across a 
#'   genomic partition. Each row contains observations from a single window of
#'   the genomic partition. Each column contains observations from a single
#'   sample (either ChIP or control/input).
#' @param design A design matrix for the full model, including a column that
#'   indicates whether the observation is a ChIP sample (1) or control/input
#'   (0). The number of rows must be \code{ncol(counts)}. Means are modeled with
#'   a log link function.
#' @param chip.col.indicator A binary vector of length
#'   \code{ncol(design.matrix)} that indicates which column of the full design
#'   matrix corresponds to the ChIP indicator.
#' @param log.offset A vector of log-scale, additive factors used to adjust
#'   estimated log-scale means for differences in library sizes across samples. 
#'   Commonly used offsets include \code{log.offset=log(colSums(counts))} and 
#'   \code{log.offset=log(apply(counts[rowSums(counts)!=0,],2,quantile,.75))}. 
#'   If \code{NULL}, the later offset is used.
#' @param print.progress logical. If TRUE, the function will provide an update
#'   on which window (row number) is being analyzed. Updates occur frequently to
#'   start then eventually occur every 5000 windows. 
#' @param bias.fold.tolerance A numerical value no smaller than 1. If the bias
#'   reduction of maximum likelihood estimates of (log) fold change is likely to
#'   result in a ratio of fold changes greater than this value, then bias
#'   reduction will be performed on such windows. Setting
#'   \code{bias.fold.tolerance=Inf} will completely disable bias reduction; 
#'   setting \code{bias.fold.tolerance=1} will always perform bias reduction
#'   (see details). If the constrained estimate differs from the unconstrained
#'   estimate, bias reduction is not performed.
#'   
#' @return A list containing: 
#'   \item{dev}{Vector containing the deviance for each
#'   window under a Poisson model fit to design matrix specified by \code{design}. This vector is passed
#'   along within the \code{QL.fit} function.} 
#'   \item{means}{Matrix of fitted mean values for each window.} 
#'   \item{parms}{Matrix of estimated coefficients for each window.} 
#'   \item{dev.constrained}{Same as \code{dev}, subject to the
#'   constraint that the ChIP coefficient is non-negative. If fitting a reduced
#'   model with no ChIP coefficient, this will be \code{NA}.} 
#'   \item{means.constrained}{Same as \code{means}, subject to
#'   the constraint that the ChIP coefficient is non-negative. If fitting a
#'   reduced model with no ChIP coefficient, this will be \code{NA}.} 
#'   \item{parms.constrained}{Same as \code{parms},
#'   subject to the constraint that the ChIP coefficient is non-negative. If
#'   fitting a reduced model with no ChIP coefficient, this will be \code{NA}.}
#'   
#' @details This functions fits, for each row of \code{counts}, a Poisson
#' log-linear model. 
#' 
#' Asymptotic biases of regression coefficients (i.e., log fold changes) are
#' then estimated by a plug-in estimate [eqn. (15.4) of McCullagh and Nelder,
#' 1989] from the last iteration of iteratively reweighted least squares (IWLS)
#' procedure. The fitted response values are then compared with or without such
#' a bias term. If the ratio of fitted response values are larger than
#' \code{bias.fold.tolerance} for any observation and the unconstrained
#'  estimate equals the constrained estimate, then the bias-reduction (not
#' bias-correction) procedure according to Firth (1993) and Kosmidis & Firth
#' (2009) is applied to such rows of \code{counts}, by adjusting the score
#' equation with a term based on the observed information. Such bias-reduced
#' estimates are more advantageous than directly subtracting the estimated bias
#' from the maximum likelihood estimates as the latter may not exist (e.g., when
#' all counts in one treatment group are zeros).
#' 
#' When the ChIP coefficient is constrained to be non-negative, quadratic 
#' programming is applied during IWLS using \code{\link{solve.QP}}. Note that if
#' the constrained estimate of the regression coefficient differs from the 
#' unconstrained estimate for a given window, bias reduction is not performed 
#' for that window.
#' 
#' @references Firth (1993) "Bias reduction of maximum likelihood estimates"
#' \emph{Biometrika}, \bold{80}, 27--38.
#' 
#' Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear
#' models" \emph{Biometrika}, \bold{96}, 793--804.
#' 
#' Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression
#' in RNA-sequence data using quasi-likelihood with shrunken dispersion
#' estimates" emph{SAGMB}, \bold{11}(5).
#' 
#' McCullagh and Nelder (1989) "Generalized Linear Models", 2nd edition. London:
#' Chapman and Hall.
#' 
#' @author Emily Goren (\email{emily.goren@gmail.com}) based on modifications of
#'   code by Steve Lund.
#'   
#' @export
#' 

PoisDev<-function(counts,design,log.offset,print.progress=TRUE,bias.fold.tolerance=1.10, chip.col.indicator){
  n<-ncol(counts)
  logFcCorrection=abs(log(bias.fold.tolerance))
  log10=log(10)
  ### For one factor designs, MLE's for Poisson GLM have simple closed form, so
  ### we can avoid one-window-at-a-time GLM fitting.
  if(is.matrix(design)) {if(ncol(design) == 1) design <- as.vector(design)}
  if(is.vector(design)) {
    
    ##offset should NOT be given on log scale here
    offset<-rep(1,length(design)); if(!is.null(log.offset)) offset<-exp(log.offset)
    
    counts<-as.matrix(counts)
    means<-counts
    
    parms<-NULL
    for(i in 1:length(unique(design))){
      if(sum( design==unique(design)[i])>1) means[,design==unique(design)[i]]<-rowSums(counts[,design==unique(design)[i]])/sum(offset[design==unique(design)[i]])
      if(sum( design==unique(design)[i])==1) means[,design==unique(design)[i]]<-counts[,design==unique(design)[i]]/offset[design==unique(design)[i]]
      parms<-cbind(parms,means[,design==unique(design)[i]][,1])
    }
    parms<-cbind(log(parms[,1]),log(parms[,-1]/parms[,1]))
    colnames(parms)<-c("(Intercept)",unique(design)[-1])
    means<-t(t(means)*offset)
    
    ### 0's require special attention since 0^0=1, but 0*log(0)=NaN
    deviance<-means-counts
    deviance[counts!=0]<-deviance[counts!=0]+(counts*log(counts/means))[counts!=0]
    deviance<-2*rowSums(deviance)
    
    ### If there is no coefficient corresponding to ChIP/input, return 'NA' for constrained estimates.
    deviance.constrained<-rep(NA,nrow(counts))
    means.constrained<-matrix(NA,nrow(counts),ncol(counts))
    parms.constrained<-rep(NA,nrow(counts))
  }
  
  if(is.matrix(design)) {
    ### For multi-factor designs, the first column of each element (matrix) of design.list should be a column of 1's, pertaining to the intercept.
    
    deviance<-deviance.constrained<-rep(NA,nrow(counts))
    means<-means.constrained<-matrix(NA,nrow(counts),ncol(counts))
    parms<-parms.constrained<-matrix(NA,nrow(counts),ncol(design))
    
    ### For each window and given design matrix, fit GLM to find model parameters (for mean structure) that optimize quasi-likelihood
    for(i in 1:nrow(counts)){
      ### If wanted, provide running progress update (eventually once every 5000 windows)
      if(i%in%c(2,10,100,500,1000,2500,4000,5000*(1:200))&print.progress) print(paste("Analyzing Window #",i))
      
      ### Fit GLM
      res<-withCallingHandlers(
        glm(counts[i,]~design[,-1],family="poisson",offset=log.offset,method=glm.fit3)	##offset should be given on log scale here
        , simpleWarning=ignorableWarnings
      )
      
      if (length(chip.col.indicator) == ncol(design)) {
        if ( res$coefficients[which(chip.col.indicator == 1)] < 0 ) {
          ### Fit GLM with constraints.
          Amat = t(chip.col.indicator)
          res.constrained <- orglm.fit(x = design, y = counts[i, ], offset = log.offset,
                                       family = poisson(), nec = 0, constr = Amat, rhs = 0)
        } else {
          res.constrained <- res
        }}
      
      ### Save optimized means (used in Pearson's dispersion estimator)
      ### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
      
      means[i,]<-res$fitted.values
      parms[i,]<-res$coefficients
      deviance[i]<-res$deviance
      means.constrained[i,]<-res.constrained$fitted.values
      parms.constrained[i,]<-res.constrained$coefficients
      deviance.constrained[i]<-res.constrained$deviance
      
      
      ### When a group has only zero counts, the MLE doesn't exist.  For these genes, we use the method Kosmidis & Firth(2009)
      ### to moderate the parameter estimates.
      ### For 2000 Fly genes, fbr takes roughly 41 seconds
      #if (any(counts[gn, ] == 0) && any(abs(coefficients(glm.fitted)) > 3)) {
      glm.fitted = res
      tmp.bias = coef(glm.fitted, type='bias')
      tmp.nonNA=which(!is.na(tmp.bias) & !is.na(glm.fitted$coefficients))
      this.x=model.matrix(glm.fitted)
      if(any(abs(this.x[,tmp.nonNA,drop=FALSE]%*%tmp.bias[tmp.nonNA]) > logFcCorrection)) {
        st = glm.fitted$coefficients - pmax.int(-log10, pmin.int(log10, tmp.bias))
        wt.idx=which(glm.fitted$weights>0)
        fbrglm.fit <- withCallingHandlers(
          glm(counts[i, wt.idx]~this.x[wt.idx,-1,drop=FALSE],family="poisson",offset=log.offset[wt.idx],start = st, method=fbrPoisglm)	##offset should be given on log scale here
          , simpleWarning=ignorableWarnings
        )
        parms[i, ] <- fbrglm.fit$coefficients
        means[i, ] <- fbrglm.fit$fitted.values
        if (length(chip.col.indicator) == ncol(design)) {
          if (!(parms[i, which(chip.col.indicator == 1)] < 0)) {
            parms.constrained[i, ] <- parms[i, ]
            means.constrained[i, ] <- means[i, ]
          }
        }
        
        tmp.nonNA = which(!is.na(fbrglm.fit))
      }
      
      
    }
  }
  
  return(list(dev=deviance, means=means, parms=parms,
              dev.constrained=deviance.constrained, 
              means.constrained=means.constrained, 
              parms.constrained=parms.constrained))
}
