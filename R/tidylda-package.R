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
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom stringr str_replace_all
#' @importFrom textmineR CalcProbCoherence
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr %>%
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

#' @importFrom generics generate
#' @export
generics::generate

