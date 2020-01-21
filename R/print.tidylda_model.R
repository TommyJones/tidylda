#' Print Method for tidylda_model
#' @description Print a summary for objects of class \code{tidylda_model}
#' @param x an object of class \code{tidylda_model}
#' @param digits minimal numer of significant digits
#' @param ... further arguments passed to or from other methods
#' @examples 
#' \dontrun{
#' dtm <- textmineR::nih_sample_dtm
#' 
#' lda <- fit_tidylda(dtm = dtm, k = 10, iterations = 100)
#' 
#' print(lda)
#' 
#' lda
#' 
#' print(lda, digits = 2)
#' }
#' @export
print.tidylda_model <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  # make this fail elegantly for legacy topic model objects
  
  ### check consistency of inputs ----
  if (class(x) != "tidylda_model")
    stop("'x' must be of class tidylda_model")
  

  ### assemble a bunch of stuff to print back at you
  
  s <- x$summary
  
  cat("A Latent Dirichlet Allocation Model with", nrow(x$phi), "topics:\n")
  
  if ("r2" %in% names(x)) {
    cat("The model's R-squared is ", round(x$r2, digits = digits), "\n")
  }
  
  cat("The most prevalent topics are:\n")
  
  tibble:::print.tbl(s[order(s$prevalence, decreasing = TRUE), ])
  
  cat("\n")
  
  cat("The most coherent topics are:\n")
  
  tibble:::print.tbl(s[order(s$coherence, decreasing = TRUE), ])
  

  
  invisible(x)
  
}
