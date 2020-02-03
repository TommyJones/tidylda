#' Glance method for \code{tidylda} objects
#' @description
#'   \code{glance} constructs a single-row summary "glance" of a \code{tidylda}
#'   topic model.
#' @param x an object of class \code{tidylda}
#' @param ... other arguments passed to methods,currently not used
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
#' dtm <- textmineR::nih_sample_dtm
#'
#' lda <- tidylda(dtm = dtm, k = 10, iterations = 100, burnin = 75)
#'
#' glance(lda)
#' }
#' @export
glance.tidylda <- function(x, ...) {

  # get some objects to return
  if (inherits(x$call, "call")) {
    call_params <- as.list(x$call)

    if (!"burnin" %in% names(call_params)) {
      call_params$burnin <- NA
    }
  } else {
    call_params <- list(
      iterations = NA,
      burnin = NA
    )
  }

  out <- data.frame(
    num_topics = nrow(x$phi),
    num_documents = nrow(x$theta),
    num_tokens = ncol(x$phi),
    iterations = call_params$iterations,
    burnin = call_params$burnin,
    stringsAsFactors = FALSE
  )

  tibble::as_tibble(out)
}

#' Tidy a matrix from a \code{tidylda} topic model
#' @description
#' Tidy the result of a \code{tidylda} topic model
#' @param x an object of class \code{tidylda} or an individual \code{phi}, 
#'   \code{theta}, or \code{gamma} matrix.
#' @param matrix the matrix to tidy; one of \code{'phi'}, \code{'theta'}, or
#'   \code{'gamma'}
#' @param log do you want to have the result on a log scale? Defaults to \code{FALSE}
#' @param ... other arguments passed to methods,currently not used
#' @return
#'   Returns a \code{\link[tibble]{tibble}}.
#'
#'   If \code{matrix = "phi"} then the result is a table of one row per topic
#'   and token with the following columns: \code{topic}, \code{token}, \code{phi}
#'
#'   If \code{matrix = "theta"} then the result is a table of one row per document
#'   and topic with the following columns: \code{document}, \code{topic}, \code{theta}
#'
#'   If \code{matrix = "gamma"} then the result is a table of one row per topic
#'   and token with the following columns: \code{topic}, \code{token}, \code{gamma}
#' @note
#'   If \code{log = TRUE} then "log_" will be appended to the name of the third
#'   column of the resulting table. e.g "\code{phi}" becomes "\code{log_phi}".
#' @examples
#' \donttest{
#' dtm <- textmineR::nih_sample_dtm
#'
#' lda <- tidylda(dtm = dtm, k = 10, iterations = 100, burnin = 75)
#'
#' tidy_phi <- tidy(lda, matrix = "phi")
#'
#' tidy_theta <- tidy(lda, matrix = "theta")
#'
#' tidy_gamma <- tidy(lda, matrix = "gamma")
#' }
#' @export
tidy.tidylda <- function(x, matrix, log = FALSE, ...) {
  
  if (!inherits(matrix, "character") |
      !sum(c("phi", "theta", "gamma") %in% matrix) >= 1) {
    stop("matrix should be one of c('phi', 'theta', 'gamma')")
  }
  
  if (matrix == "phi") {
    out <- tidy.matrix(x = x$phi, matrix = matrix, log = log)
  } else if (matrix == "gamma") {
    out <- tidy.matrix(x = x$gamma, matrix = matrix, log = log)
  } else {
    out <- tidy.matrix(x = x$theta, matrix = matrix, log = log)
  }
  
  out
}

#' Tidy an individual matrix. Useful for predictions and called from tidy.tidylda
#' @describeIn tidy.tidylda Tidy an individual matrix. Useful for predictions and called from tidy.tidylda
#' @export
tidy.matrix <- function(x, matrix, log = FALSE, ...) {
  
  # check inputs
  if (!inherits(matrix, "character") |
      !sum(c("phi", "theta", "gamma") %in% matrix) >= 1) {
    stop("matrix should be one of c('phi', 'theta', 'gamma')")
  }
  
  if (!is.logical(log)) {
    stop("log must be logical.")
  }
  
  out <- data.frame(
    names_col = rownames(x),
    x,
    stringsAsFactors = FALSE
  )
  
  out <- tidyr::pivot_longer(
    data = out, cols = setdiff(colnames(out), "names_col"),
    names_to = "index", values_to = "value"
  )
  
  if (matrix == "phi") {
    
    colnames(out) <- c("topic", "token", "phi")
    
    out$topic <- as.numeric(out$topic)
    
  } else if (matrix == "gamma") {
    
    colnames(out) <- c("topic", "token", "gamma")
    
    out$topic <- as.numeric(out$topic)

  } else { # meanse matrix  == theta
    
    colnames(out) <- c("document", "topic", "theta")
    
    out$topic <- as.numeric(out$topic)
    
  }
  
  if (log) {
    out[[3]] <- log(out[[3]])
    
    names(out)[3] <- paste0("log_", names(out)[3])
  }
  
  out
  
}
