#'
#' @import edgeR
#' @import quadprog
#' @import GenomicRanges
#' @import IRanges
#' @import GenomicAlignments
#' @import SummarizedExperiment
#' @import csaw
#' @import Rsamtools
#'
#' @importFrom mgcv negbin
#' @importFrom pracma broyden
#' @importFrom grDevices dev.new
#' @importFrom graphics layout legend lines par plot
#' @importFrom stats coef coefficients gaussian glm glm.control lm.fit model.matrix optimize p.adjust pf poisson predict qf smooth.spline
#' 
#' @useDynLib BinQuasi, .registration=TRUE, .fixes="C_"
#' 
NULL

# @S3method coef, glm
# @S3method unique, matrix
# @S3method duplicated, matrix