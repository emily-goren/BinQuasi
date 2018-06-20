#'
#' @name BinQuasi
#'
#' @title Analyzing Replicated ChIP Sequencing Data Using Quasi-Likelihood
#' 
#' @description Identify peaks in ChIP-seq data with biological replicates.
#' 
#' @details Identify peaks in ChIP-seq data with biological replicates using a one-sided quasi-likelihood ratio test in quasi-Poisson or quasi-negative binomial models.
#' 
#' @import edgeR
#' @import quadprog
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import SummarizedExperiment
#' @import csaw
#' @import BiocGenerics
#' @import S4Vectors
#' @import Rsamtools
#' @import RMySQL
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