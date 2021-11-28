#' Draw from the marginal posteriors of a tidylda topic model
#' @description These functions are used to sample from the marginal posteriors
#'   of a \code{tidylda} topic model. This is useful for quantifying uncertainty
#'   around the parameters of \code{beta} or \code{theta}.
#' @param x For \code{posterior}, an object of class
#'   \code{tidylda}. For \code{generate}, an object of class
#'   \code{tidylda_posterior} obtained by a call to \code{posterior}.
#' @param matrix A character of either 'theta' or 'beta', indicating from which
#'   matrix to draw posterior samples.
#' @param which Row index of \code{theta}, for document, or \code{beta}, for
#'   topic, from which to draw samples. \code{which} may also be a vector of
#'   indices.
#' @param times Integer number of samples to draw.
#' @param ... Other arguments, currently not used.
#' @return 
#' \code{posterior} returns an object of class \code{tidylda_posterior}.
#' 
#' \code{generate} returns a tibble with one row per parameter per sample.
#' @details
#' To sample from the marginal posteriors of a model, you must first make a call
#' to \code{posterior} and then a call to \code{generate}.
#' 
#' \code{posterior} takes an object of class \code{tidylda} and constructs an
#' object of class \code{tidylda_posterior} which contains two matrices. The rows
#' of these matrices are Dirichlet parameters used to sample from the marginal
#' posteriors of \code{theta} and \code{beta}.
#' 
#' \code{generate} takes an object of class \code{tidylda_posterior} and samples
#' from the marginal posterior of the parameters specified by the \code{matrix}
#' argument.
#' @references 
#' Heinrich, G. (2005) Parameter estimation for text analysis. Technical report. 
#' \href{http://www.arbylon.net/publications/text-est.pdf}{http://www.arbylon.net/publications/text-est.pdf}
#' @return Returns a data frame where each row is a single sample from the posterior. 
#' Each column is the distribution over a single parameter. The variable \code{var}
#' is a facet for subsetting by document (for theta) or topic (for beta).
#' @export
#' @examples
#' \donttest{
#' # load some data
#' data(nih_sample_dtm)
#'
#' # fit a model
#' set.seed(12345)
#'
#' m <- tidylda(
#'   data = nih_sample_dtm[1:20, ], k = 5,
#'   iterations = 200, burnin = 175
#' )
#' 
#' # construct a posterior object
#' p <- posterior(m)
#' 
#' # sample from the marginal posterior corresponding to topic 1
#' t1 <- generate(
#'   x = p,
#'   matrix = "beta",
#'   which = 1,
#'   times = 100  
#' )
#' 
#' # sample from the marginal posterior corresponding to document 5
#' d5 <- generate(
#'   x = p,
#'   matrix = "theta",
#'   which = 5,
#'   times = 100
#' )
#' }
#' @export
posterior <- function(x, ...) UseMethod("posterior")

#' Posterior method for tidylda
#' @rdname posterior
#' @export
posterior.tidylda <- function(x, ...) {
  
  # get proper alpha
  alpha <- format_alpha(
    x$alpha, 
    k = nrow(x$beta)
  )
  
  # get proper eta
  eta <- format_eta(
    x$eta, 
    k = nrow(x$beta), 
    Nv = ncol(x$beta)
  )
  
  # extract dirichlet parameters for theta
  theta_par <- t(x$counts$Cd) + alpha$alpha
  
  colnames(theta_par) <- rownames(x$theta)
  rownames(theta_par) <- colnames(x$theta)
  
  # extract dirichlet parameters for beta
  beta_par <- x$counts$Cv + eta$eta
  
  rownames(beta_par) <- rownames(x$beta)
  colnames(beta_par) <- colnames(x$beta)
  
  beta_par <- t(beta_par)
  
  # prepare and return result
  result <- list(
    theta_par = theta_par,
    beta_par = beta_par
  )
  
  class(result) <- "tidylda_posterior"
  
  result
}

#' @rdname posterior
#' @export
generate.tidylda_posterior <- function(
  x,
  matrix,
  which,
  times,
  ...
) {
  
  # check inputs
  if (! matrix[1] %in% c("theta", "beta")) {
    stop("matrix must be one of 'theta' or 'beta'")
  }
  
  if (any(is.na(which)) | any(is.infinite(which))) {
    stop("NA or Inf detected! which cannot have any NA or infinite values.")
  }
  
  if (any(which <= 0)) {
    stop("Negative or 0 detected! which must contain positive integers")
  }
  
  if (length(times) > 1 | length(times) == 0) {
    stop("times must be a positive number of length 1.")
  }
  
  if (! is.numeric(times)) {
    stop("times must be a positive number.")
  }
  
  if (times <= 0) {
    stop("times must be a positive number.")
  }
  
  # extract the right x 
  if (matrix[1] == "theta") {
    obj <- x$theta_par
    # do one more check on which
    if (any(which > ncol(obj))) {
      stop("which is contains values greater than the maximum number of documents.")
    }
  } else {
    obj <- x$beta_par
    # do one more check on which
    if (any(which > ncol(obj))) {
      stop("which is contains values greater than the maximum number of topics.")
    }
  }
  
  # reformat for efficient sampling
  result <- as.data.frame(obj[, which], stringsAsFactors = FALSE)
  
  # sample
  result <- lapply(
    result,
    function(y) {
      samp <- gtools::rdirichlet(n = times, alpha = y)
      
      samp <- as.data.frame(samp, stringsAsFactors = FALSE)
      
      colnames(samp) <- rownames(obj)

      out <- as.data.frame(t(samp))
      
      colnames(out) <- 1:ncol(out)
      
      out$idx1 <- rownames(obj)
      
      out <- tidyr::pivot_longer(
        out,
        -.data$idx1,
        names_to = "sample"
      )
      
      out
    }
  )
  
  # reformat the output so it's tidy
  for (j in seq_along(result)) {
    result[[j]]$idx2 <- colnames(obj)[which[j]]
  }
  
  result <- do.call(rbind, result)

  result <- tibble::as_tibble(result)
  
  result$sample <- as.numeric(result$sample)
  
  if (matrix[1] == "theta") {
    names(result)[names(result) == "idx1"] <- "topic" 
    names(result)[names(result) == "idx2"] <- "document"
    names(result)[names(result) == "value"] <- "theta"
    
    result$topic <- as.numeric(result$topic)
    
    result <- result[, c("document", "topic", "sample", "theta")]
    
  } else {
    names(result)[names(result) == "idx1"] <- "token" 
    names(result)[names(result) == "idx2"] <- "topic" 
    names(result)[names(result) == "value"] <- "beta"
    
    result$topic <- as.numeric(result$topic)
    
    result <- result[, c("topic", "token", "sample", "beta")]
  }
  
  result
}