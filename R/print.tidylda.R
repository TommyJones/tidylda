#' Print Method for tidylda
#' @description Print a summary for objects of class \code{tidylda}
#' @param x an object of class \code{tidylda}
#' @param digits minimal number of significant digits
#' @param n Number of rows to show in each displayed \code{\link[tibble]{tibble}}.
#' @param ... further arguments passed to or from other methods
#' @examples
#' \donttest{
#' dtm <- textmineR::nih_sample_dtm
#'
#' lda <- tidylda(dtm = dtm, k = 10, iterations = 100)
#'
#' print(lda)
#'
#' lda
#'
#' print(lda, digits = 2)
#' }
#' @export
print.tidylda <- function(x, digits = max(3L, getOption("digits") - 3L), n = 5, ...) {

  # make this fail elegantly for legacy topic model objects

  ### check consistency of inputs ----


  ### assemble a bunch of stuff to print back at you

  s <- x$summary

  cat(
    "A Latent Dirichlet Allocation Model of ", nrow(x$phi), "topics, ",
    nrow(x$theta), " documents, and ", ncol(x$phi), " tokens:\n"
  )

  print(x$call)

  cat("\n")

  if ("r2" %in% names(x)) {
    cat("The model's R-squared is ", round(x$r2, digits = digits), "\n")
  }

  cat("The ", n, " most prevalent topics are:\n")

  print(s[order(s$prevalence, decreasing = TRUE), ], n = n)

  cat("\n")

  cat("The ", n, " most coherent topics are:\n")

  print(s[order(s$coherence, decreasing = TRUE), ], n = n)



  invisible(x)
}
