#' Glance method for \code{tidylda} objects
#' @description
#'   \code{glance} constructs a single-row summary "glance" of a \code{tidylda}
#'   topic model.
#' @param x an object of class \code{tidylda}
#' @param ... other arguments passed to methods
#' @examples
#' \donttest{
#'   dtm <- textmineR::nih_sample_dtm
#'   
#'   lda <- tidylda(dtm = dtm, k = 10, iterations = 100, burnin = 75)
#'   
#'   glance(lda)
#'   
#' }
glance.tidylda <- function(x, ...) {
  
  out <- data.frame(num_topics = nrow(x$phi),
                    num_documents = nrow(x$theta),
                    vocab_size = ncol(x$phi),
                    avg_topic_coherence = mean(x$summary$coherence),
                    stringsAsFactors = FALSE)
  
  tibble::as_tibble(out)
  
}
