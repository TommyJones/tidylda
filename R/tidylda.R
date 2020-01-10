#' tidylda
#'
#' Latent Dirichlet Allocation Using 'tidyverse' Conventions
#'
#' Implements the WarpLDA algorithm for Latent Dirichlet Allocation 
#' using style conventions from the 'tidyverse' and specifically 'tidymodels'. 
#' Designed to work seamlessly with 'textmineR'.
#'
#' @name tidylda
#' @docType package
NULL


#' @import Matrix
#' @import Rcpp
#' @import textmineR
#' @importFrom gtools rdirichlet
#' @importFrom Rcpp sourceCpp
#' @importFrom methods as
#' @useDynLib "tidylda", .registration=TRUE
NULL


