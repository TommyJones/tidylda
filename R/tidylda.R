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


#' @import Matrix
#' @import Rcpp
#' @import textmineR
#' @importFrom gtools rdirichlet
#' @importFrom methods as
#' @importFrom Rcpp sourceCpp
#' @importFrom stats median
#' @useDynLib "tidylda", .registration=TRUE
NULL


