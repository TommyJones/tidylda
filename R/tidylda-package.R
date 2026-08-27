#' Latent Dirichlet Allocation Using 'tidyverse' Conventions
#' @description
#'   Implements an algorithm for Latent Dirichlet Allocation (LDA)
#'   using style conventions from the 'tidyverse' and specifically 'tidymodels'.
#'   Also implements several novel features for LDA such as guided models and
#'   transfer learning.
#'
#'   Fitting uses warpLDA (Chen et al., 2016,
#'   \doi{10.48550/arXiv.1510.08628}), a Metropolis-Hastings sampler that
#'   alternates document-ordered and word-ordered passes over the corpus so that
#'   each pass touches only a small, cache-resident working set. It replaced the
#'   collapsed Gibbs sampler in version 0.1.0. Sampling is multithreaded via the
#'   \code{threads} argument and results do not depend on the thread count.
#'
#' @section Options:
#'   \code{tidylda.max_result_size} caps the size of the object
#'   \code{\link[tidylda]{posterior.tidylda}} and
#'   \code{\link[tidylda]{tidy.tidylda}} will build, in bytes. Both return one
#'   row per cell of a matrix that grows with topics times vocabulary, so a
#'   plausible-looking call can ask for billions of rows; above the cap they
#'   raise an error naming the size and a smaller alternative rather than
#'   exhausting the session. Defaults to \code{1024^3} (1 GB). Raise it with
#'   \code{options(tidylda.max_result_size = 4 * 1024^3)}.
#' @name tidylda-package
#' @keywords internal
"_PACKAGE"

# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("tidylda is under active development. The API and behavior may change.")
# }

#' @import Rcpp
#' @importFrom gtools rdirichlet
#' @importFrom methods as
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom stringr str_replace_all
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

#' Abstracts and metadata from NIH research grants awarded in 2014
#' @name nih
#' @aliases nih_sample
#' @aliases nih_sample_dtm
#' @docType data
#' @description
#' This dataset holds information on research grants awarded by the National
#' Institutes of Health (NIH) in 2014. The data set was downloaded in
#' approximately January of 2015. It includes both 'projects' and 'abstracts' 
#' files.
#' @usage 
#' data("nih_sample")
#' @format
#' For \code{nih_sample}, a \code{\link[tibble]{tibble}} of 100 randomly-sampled
#' grants' abstracts and metadata. For \code{nih_sample_dtm}, a
#' \code{\link[Matrix]{dgCMatrix-class}} representing the document term matrix
#' of abstracts from 100 randomly-sampled grants.
#' @source
#' National Institutes of Health ExPORTER
#' \url{https://reporter.nih.gov/exporter}
NULL


