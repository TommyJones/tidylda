#' Glance method for \code{tidylda} objects
#' @description
#'   \code{glance} constructs a single-row summary "glance" of a \code{tidylda}
#'   topic model.
#' @param x an object of class \code{tidylda}
#' @param ... other arguments passed to methods
#' @return
#'   \code{glance} returns a one-row \code{\link[tibble]{tibble}} with the 
#'   following columns:
#'   
#'   \code{num_topics}: the number of topics in the model
#'   \code{num_documents}: the number of documents used for fitting
#'   \code{num_tokens}: the number of tokens covered by the model
#'   \code{iterations}: number of total Gibbs iterations run
#'   \code{burnin}: number of burn-in Gibbs iterations run
#' @examples
#' \donttest{
#'   dtm <- textmineR::nih_sample_dtm
#'   
#'   lda <- tidylda(dtm = dtm, k = 10, iterations = 100, burnin = 75)
#'   
#'   glance(lda)
#'   
#' }
#' @export
glance.tidylda <- function(x, ...) {
  
  # get some objects to return
  if (inherits(x$call, "call")) {
    
    call_params <- as.list(x$call)
    
    if (! "burnin" %in% names(call_params))
      call_params$burnin <- NA
    
  } else {
    
    call_params <- list(iterations = NA,
                        burnin = NA)
    
  }
  
  out <- data.frame(num_topics = nrow(x$phi),
                    num_documents = nrow(x$theta),
                    num_tokens = ncol(x$phi),
                    iterations = call_params$iterations,
                    burnin = call_params$burnin,
                    stringsAsFactors = FALSE)
  
  tibble::as_tibble(out)
  
}


