#' tidylda
#'
#' Latent Dirichlet Allocation Using 'tidyverse' Conventions
#'
#' Implements the WarpLDA algorithm for Latent Dirichlet Allocation 
#' using style conventions from the 'tidyverse' and specifically 'tidymodels'. 
#' Designed to work seamlessly with 'textmineR'.
#'
#' @name tidylda-package
#' @docType package
NULL



#' @import Rcpp
#' @importFrom textmineR CalcGamma
#' @importFrom textmineR CalcProbCoherence
#' @importFrom textmineR CalcTopicModelR2
#' @importFrom textmineR TmParallelApply
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom Matrix rbind2
#' @importFrom Matrix cbind2
#' @importFrom gtools rdirichlet
#' @importFrom methods as
#' @importFrom Rcpp sourceCpp
#' @importFrom stats median
#' @useDynLib "tidylda", .registration=TRUE
NULL


