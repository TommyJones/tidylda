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
#' @importFrom gtools rdirichlet
#' @importFrom methods as
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom Matrix rbind2
#' @importFrom Matrix cbind2
#' @importFrom Rcpp sourceCpp
#' @importFrom stats median
#' @importFrom stringr str_replace_all
#' @importFrom textmineR CalcGamma
#' @importFrom textmineR CalcProbCoherence
#' @importFrom textmineR CalcTopicModelR2
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidytext cast_dfm
#' @importFrom tidytext cast_dtm
#' @useDynLib "tidylda", .registration=TRUE
NULL

#' @importFrom generics augment
#' @export
generics::augment

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics refit
#' @export
generics::refit
