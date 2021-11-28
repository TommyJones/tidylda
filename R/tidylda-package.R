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

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("tidylda is under active development. The API and behavior may change.")
}

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

#' @importFrom generics generate
#' @export
generics::generate

#' Abstracts and metadata from NIH research grants awarded in 2014
#' @name nih
#' @aliases nih_sample
#' @aliases nih_sample_dtm
#' @docType data
#' @description
#' This dataset holds information on research grants awarded by the National
#' Institutes of Health (NIH) in 2014. The data set was downloaded in
#' approximately January of 2015 from
#' \url{https://exporter.nih.gov/ExPORTER_Catalog.aspx}. It includes both
#' 'projects' and 'abstracts' files.
#' @usage 
#' data("nih_sample")
#' @format
#' For \code{nih_sample}, a \code{\link[tibble]{tibble}} of 100 randomly-sampled
#' grants' abstracts and metadata. For \code{nih_sample_dtm}, a
#' \code{\link[Matrix]{dgCMatrix-class}} representing the document term matrix
#' of abstracts from 100 randomly-sampled grants.
#' @source
#' National Institutes of Health ExPORTER
#' \url{https://exporter.nih.gov/ExPORTER_Catalog.aspx}
NULL


