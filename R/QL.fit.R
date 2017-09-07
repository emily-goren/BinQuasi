
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
