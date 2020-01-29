#' Latent Dirichlet Allocation Using 'tidyverse' Conventions
#' @description
#'   Implements an algorithm for Latent Dirichlet Allocation (LDA)
#'   using style conventions from the 'tidyverse' and specifically 'tidymodels'. 
#'   Also implements several novel features for LDA such as guided models and 
#'   transfer learning.
#' @name tidylda-package
#' @docType package
#' @keywords internal
NULL



#' @import Rcpp
#' @import furrr
#' @import tibble
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


