PoisDev<-function(counts,design,log.offset,print.progress=TRUE,bias.fold.tolerance=1.10, chip.col.indicator){
    n<-ncol(counts)
        logFcCorrection=abs(log(bias.fold.tolerance))
    log10=log(10)
### For one factor designs, MLE's for Poisson GLM have simple closed form, so we can avoid one-gene-at-a-time GLM fitting.
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

		### For each gene and given design matrix, fit GLM to find model parameters (for mean structure) that optimize quasi-likelihood
		for(i in 1:nrow(counts)){
			### If wanted, provide running progress update (eventually once every 5000 genes)
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
                dev.constrained=deviance.constrained, means.constrained=means.constrained, parms.constrained=parms.constrained))
}
