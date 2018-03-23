
#'
#'Obtain p- and q-values using results from \code{QL.fit}
#'
#'@description Obtain significance results for quasi-likelihood models fit to
#'  ChIP-seq data partitioned into counts using the methods detailed in Goren,
#'  Liu, Wang and Wang (2018).
#'  
#'@details Obtain significance results from an object fitted using \code{\link{QL.fit}}.
#'Used within the main peak calling function, \code{\link{BQ}}.
#'  
#'@param fit The list returned by the function \code{QL.fit}.
#'@param Dispersion Must be one of \code{"Deviance"} or \code{"Pearson"}, specifying which
#'  type of estimator should be used for estimating the quasi-likelihood
#'  dispersion parameter.
#'@param one.sided logical. If TRUE, a one-sided test for the ChIP coefficient
#'  is reported. Otherwise, if FALSE, a two-sided test is reported.
#'@param spline.df Optional. User may specify the degrees of freedom to use when
#'  fitting a cubic spline to log-scale(estimated dispersion) versus the
#'  log(average count). Default uses cross-validation in \code{sreg} function to
#'  pick appropriate degrees of freedom.
#'@param Plot logical. If TRUE, the estimated dispersion versus the average
#'  count are plotted on a log-scale with the corresponding cubic spline fit
#'  overlaid.
#'@param padj logical. If TRUE, Benjamini & Hochberg's adjustment for multiple
#'  comparisons is applied.
#'  
#'@return  A list containing: 
#'  \item{P.values}{List of matrices providing
#'  p-values for the QL, QLShrink and QLSpline methods, respectively. The i-th
#'  column of each element of \code{pvals} corresponds to the hypothesis test
#'  assigned in the i-th window.} 
#'  \item{Q.values}{List of
#'  matrices providing q-values for the QL, QLShrink and QLSpline methods,
#'  respectively. The i-th column of each element of \code{qvals} corresponds to
#'  the hypothesis test assigned in the i-th window.} 
#'  \item{F.stat}{List of matrices providing F-statistics for the QL, QLShrink
#'  and QLSpline methods, respectively. The i-th column of each element of
#'  \code{F.stat} corresponds to the hypothesis test assigned in the i-th window.} 
#'  \item{d0}{Vector containing estimated additional
#'  denominator degrees of freedom  gained from shrinking dispersion estimates
#'  in the QLShrink and QLSpline procedures, respectively.}
#'  
#'@author Emily Goren (\email{emily.goren@gmail.com}) based on modifications of
#'  code by Steve Lund.
#'  
#'@references Benjamini and Hochberg (1995) "Controlling the false discovery
#'rate: a practical and powerful approach to multiple testing" \emph{Journal of
#'the Royal Statistical Society Series B}, \bold{57}: 289-300.
#'
#'Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression
#'in RNA-sequence data using quasi-likelihood with shrunken dispersion
#'estimates" \emph{SAGMB}, \bold{11}(5).
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
#'  
#'  ####################################################
#'  # Fit quasi-negative binomial model to each window.
#'  ####################################################
#'  design.matrix  <- cbind(rep(1, reps*2), # Intercept
#'                          rep(c(1,0), each = reps)) # Indicates ChIP sample
#'  chip.col.indicator <- c(0,1) # Second column of design matrix indicates ChIP sample
#'  fit <- QL.fit(counts, design.matrix, chip.col.indicator, 
#'                log.offset = rep(1, ncol(counts)), Model = 'NegBin') 
#'  window.results <- QL.results(fit)
#' 
#' # Number of significant windows.
#' sum(window.results$Q.values$QLShrink < 0.05)
#' 
#' # Compare significant windows to truth.
#' res <- as.numeric(window.results$Q.values$QLShrink < 0.05)
#' # Number of true positives
#' TP <- sum(res[peak.idx] == 1)
#' TP
#' # Number of false negatives
#' FN <- n.peaks - TP
#' FN
#' # Number of false positives
#' FP <- sum(res) - TP
#' FP
#' # Number of true negatives
#' TN <- (K - sum(res)) - FN
#' TN
#'
#'@export
#'

QL.results<-function(fit,Dispersion="Deviance", one.sided=TRUE, spline.df=NULL,Plot=FALSE, padj=TRUE){
  if(!Dispersion %in% c("Deviance","Pearson")) 
    stop("Unidentified Dispersion: Dispersion must be either 'Deviance' or 'Pearson'.")
  
  LRT<-fit$LRT; phi.hat<-fit$phi.hat.dev; mn.cnt<-fit$mn.cnt; den.df<-fit$den.df; num.df<-fit$num.df; Model=fit$Model
  if(Dispersion=="Pearson") phi.hat<-fit$phi.hat.pearson
  
  if (one.sided == TRUE) LRT <- fit$LRT.constrained
  
  
  if(length(num.df)==1) num.df<-rep(num.df,ncol(LRT))
  
  ### We use Smyth's (2004) approach from LIMMA to estimate parameters for prior distributions
  ### of window specific dispersion estimates. The function below also provides resulting
  ### shrunken point estimates of dispersion
  
  #### Code for implementing Smyth's approach begins here ####
  shrink.phi<-function(phi.hat,den.df){
    phi.hat[phi.hat<=0]<-min(phi.hat[phi.hat>0])
    z<-log(phi.hat); z[z==Inf]<-max(z[z!=Inf]); z[z==-Inf]<-min(z[z!=-Inf]);mnz<-mean(z)
    
    ## solve for d0 and phi0
    d0arg<-var(z)-trigamma((den.df)/2)
    if(d0arg>0){
      dif<-function(x,y) abs(trigamma(x)-y)
      inverse.trigamma<-function(y) optimize(dif,interval=c(0,10000),y=y)$minimum
      d0<-2*inverse.trigamma(d0arg)
      phi0<-exp(mnz-digamma((den.df)/2)+digamma(d0/2)- log(d0/(den.df)))
      
      ## compute shrunken phi's
      phi.shrink<-((den.df)*phi.hat+d0*phi0)/(den.df+d0)
    }
    else{phi.shrink<-rep(exp(mnz),length(z)); d0<-Inf; phi0<-exp(mnz) }
    return(list(phi.shrink=phi.shrink,d0=d0,phi0=phi0))
  }
  #### Code for implementing Smyth's approach ends here ####
  
  phi.hat[phi.hat<0]<-min(phi.hat[phi.hat>0])
  phi.hat2<-phi.hat
  if(Model=="Poisson") phi.hat2[phi.hat<1]<-1
  
  shrink<-shrink.phi(phi.hat,den.df)
  phi.shrink<-shrink[[1]]; est.d0<-shrink[[2]];
  if(Model=="Poisson") phi.shrink[phi.shrink<1]<-1
  
  ### Fit cubic spline to capture relationship between log(y.bar) and log(phi.hat)
  y<-log(phi.hat); y[y==-Inf]<-min(y[y!=-Inf]); y[y==Inf]<-max(y[y!=Inf])
  
  spline.fit<-if(is.null(spline.df)) smooth.spline(x=log(mn.cnt), y=y) else smooth.spline(x=log(mn.cnt), y=y, df= spline.df)
  spline.pred<-predict(spline.fit,x=log(mn.cnt))$y
  
  fit.method<-"spline"
  
  #if(spline.fit$df==0){
  #	fit.method<-"locfit curve"
  #	print("Spline fitting failed.  Using locfit instead.")
  #	fit<-locfit(y~lp(log(mn.cnt)))
  #	spline.fit<-list(fitted.values=predict(fit,newdata=log(mn.cnt)),eff.df=fit$dp["df1"])
  #}
  
  ### Obtain estimate for prior degrees of freedom and scaling factor after adjusting for cubic spline fit
  y2<-phi.hat/exp(spline.pred)
  shrink<-shrink.phi(y2,den.df)
  D0<-shrink[[2]];
  phi0<-shrink[[3]]; print(paste("Spline scaling factor:",phi0))
  
  
  ### If desired, plot resulting cubic spline fit
  if(Plot){
    dev.new(height=9); nf <- layout(matrix(1:2, 2, 1),heights=c(7,2))
    
    par(mai=c(1,1.2,1,.2))
    suppressWarnings(plot(log(mn.cnt),y,xlab=expression(log(bar(y)[phantom()%.%phantom()]*phantom()[phantom()%.%phantom()])),
                          ylab=expression(log(hat(Phi))),main="Estimated Dispersion
 versus Average Count",pch=1,cex.lab=2,cex.axis=2,cex.main=2))
    lines(sort(log(mn.cnt)),spline.pred[order(mn.cnt)],col=2,lwd=3)
    
    ###Compare quantiles from estimated theoretical and empirical dispersion distributions
    
    
    RR<-NULL; sort.mn<-log(sort(mn.cnt)); qq<-c(.05,.95);ord.F<-y[order(mn.cnt)]
    bins<-c(1+round(length(y)*(0:19)/20),length(y))
    for(ii in 1:(length(bins)-1)){
      ind<-bins[ii]:bins[ii+1]
      RR<-rbind(RR,c(quantile(ord.F[ind],qq),median(sort.mn[ind])))
    }
    
    lines(sort(log(mn.cnt)),spline.pred[order(mn.cnt)]+log(phi0*qf(.95,den.df,D0)),col=4,lwd=3)
    lines(sort(log(mn.cnt)),spline.pred[order(mn.cnt)]+log(phi0*qf(.05,den.df,D0)),col=4,lwd=3)
    
    lines(RR[,3],RR[,2],col=3,lwd=3)
    lines(RR[,3],RR[,1],col=3,lwd=3)
    par(mai=rep(0,4));plot(19,19,col="white",axes=FALSE)
    
    suppressWarnings(legend("top",legend=c(paste("Fitted", fit.method,"with",signif(spline.fit$df,2),"df"),
                                           "0.05 & 0.95 Quantiles from Empirical Distribution",
                                           "0.05 & 0.95 Quantiles from Scaled F-Distribution"),lwd=3,lty=1,col=2:4,cex=1.5))
  }
  
  
  
  phi.spline<-(D0*exp(spline.pred)*phi0+(den.df)*phi.hat)/(D0+den.df)
  if(D0==Inf){
    warning("D0 estimate is infinity for QLSpline (there's little scatter in original dispersion estimates around fitted spline).
 QLSpline dispersions set to fitted cubic spline (use 'Plot=TRUE' to view) with no uncertainty.")
    phi.spline<-exp(spline.pred)
  }
  if(Model=="Poisson") phi.spline[phi.spline<1]<-1
  
  log.pval<-F.stat<-list(QL=NULL,QLShrink=NULL,QLSpline=NULL)
  
  
  for(i in 1:ncol(LRT)){
    #### Compute p-values for each hypothesis test using each approach
    log.pval[[1]]<-cbind(log.pval[[1]],pf(LRT[,i]/phi.hat2,num.df[i],den.df,lower.tail = FALSE,log.p=TRUE))
    log.pval[[2]]<-cbind(log.pval[[2]],pf(LRT[,i]/phi.shrink,num.df[i],est.d0+den.df,lower.tail = FALSE,log.p=TRUE))
    log.pval[[3]]<-cbind( log.pval[[3]], pf(LRT[,i]/phi.spline, num.df[i], D0+den.df, lower.tail = FALSE, log.p=TRUE))
    F.stat[[1]]<-cbind(F.stat[[1]],LRT[,i]/phi.hat2)
    F.stat[[2]]<-cbind(F.stat[[2]],LRT[,i]/phi.shrink)
    F.stat[[3]]<-cbind(F.stat[[3]],LRT[,i]/phi.spline)
  }
  
  pval<-log.pval
  for(ii in 1:3){
    pval[[ii]]<-exp(log.pval[[ii]])
    if (one.sided == TRUE)
      pval[[ii]]<- 0.5 *  pval[[ii]]
    colnames(F.stat[[ii]])<-colnames(log.pval[[ii]])<-colnames(pval[[ii]])<-colnames(LRT)
  }
  d0<-c(QLShrink=est.d0, QLSpline=D0)
  
  qval<-pval
  for(ii in 1:length(qval)){
    M0<-NULL; Qval<-qval[[ii]]
    for(jj in 1:ncol(Qval)){
      qval[[ii]][!is.na(Qval[,jj]),jj]=p.adjust(Qval[!is.na(Qval[,jj]),jj], method = "fdr")
    }
  }
  
  return(list(P.values=pval,log.P.values=log.pval,Q.values=qval,F.stat=F.stat,d0=d0))
}





