BQ <- function(dir, ChIP.files, control.files, alpha=0.05,
               bin.size=NULL, frag.length=NULL,
               minimum.count=20, Model="NegBin",
               print.progress=TRUE, method="QLShrink",
               Dispersion="Deviance", log.offset=NULL,
               NBdisp="trend", bias.fold.tolerance=1.10) {

    if (alpha <= 0 | alpha >= 0.5)
        stop("Please specify a significance level between 0 and 0.5")

    if(!Dispersion %in% c("Deviance","Pearson"))
        stop("Unidentified Dispersion: Dispersion must be either 'Deviance' or 'Pearson'.")

    if(!method %in% c("QL", "QLShrink", "QLSpline"))
        stop("Unidentified method for computing p-values: method must be either 'QL', 'QLShrink', or 'QLSpline'.")


    if (!Model %in% c("NegBin", "Poisson"))
        stop("Unidentified Model: Model must be either 'NegBin' or 'Poisson'.")

    n.C <- length(ChIP.files); n.I <- length(control.files)
    cts <- count.table(dir, ChIP.files, control.files, bin.size, frag.length, minimum.count)
    message('Using bin size of ', cts$bin.size, ' bp \n')
    Len <- cts$fragment.length
    Lenname <- c(ChIP.files, control.files)
    for (i in 1:length(Len))
         message('Using estimated fragment length for ', Lenname[i], ' equal to ', Len[i], ' bp \n')
    counts <- cts$counts
    Y <- counts[ , 4:(3+n.C+n.I)]


    if (Model == "NegBin") {
        if (!NBdisp[1] %in% c("trend", "common") & length(NBdisp) != nrow(Y))
            stop("Unidentified NegBin Dispersion: NBdisp must be set as 'trend' or 'common' to estimate negative binomial dispersion from data using GLM edgeR (McCarthy et al., 2012),\n or it must be a vector providing negative binomial dispersion parameter value to use for each window.")

        if (length(NBdisp) == nrow(Y) & !is.numeric(NBdisp))
            stop("NBdisp contains non-numeric values.\n\tAll negative binomial dispersion parameters must be non-negative real numbers.")

        if (length(NBdisp) == nrow(Y) & any(NBdisp < 0))
            stop("NBdisp contains negative values.\nAll negative binomial dispersion parameters must be non-negative real numbers.")
}


    design.matrix <- cbind(rep(1, n.C + n.I),
                           c(rep(1, n.C), rep(0, n.I)))
    fit <- QL.fit(counts = Y,
                  design.matrix = design.matrix,
                  chip.col.indicator = c(0,1),
                  log.offset = log.offset,
                  NBdisp = NBdisp,
                  print.progress = print.progress,
                  bias.fold.tolerance = bias.fold.tolerance,
                  Model = Model)
    res <- QL.results(fit = fit, Dispersion = Dispersion)
    if (method == 'QL') {
        ps <- res$P.values$QL
    } else if (method == 'QLShrink') {
        ps <- res$P.values$QLShrink
    } else if (method == 'QLSpline') {
        ps <- res$P.values$QLSpline
    }
    called <- call.peaks(window.pvals = ps,
                         start = counts[,2],
                         end = counts[,3],
                         chromosomes = counts[,1],
                         alpha = alpha)
    out <- list(peaks = called,
                bin.size = cts$bin.size,
                fragment.length = cts$fragment.length,
                filter = cts$filter)
    return(out)
}
