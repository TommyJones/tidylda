#' Draw from the marginal posteriors of a tidylda topic model
#' @description Sample from the marginal posteriors of a \code{tidylda} topic
#'   model. This is useful for quantifying uncertainty around the parameters of
#'   \code{beta} or \code{theta}.
#' @param x An object of class \code{tidylda}. 
#' @param matrix A character of either 'theta' or 'beta', indicating from which
#'   matrix to draw posterior samples.
#' @param which Row index of \code{theta}, for document, or \code{beta}, for
#'   topic, from which to draw samples. \code{which} may also be a vector of
#'   indices to sample from multiple documents or topics simultaneously.
#' @param times Integer, number of samples to draw.
#' @param ... Other arguments, currently not used.
#' @return 
#' \code{posterior} returns a tibble with one row per parameter per sample.
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
#' # sample from the marginal posterior corresponding to topic 1
#' t1 <- posterior(
#'   x = m,
#'   matrix = "beta",
#'   which = 1,
#'   times = 100  
#' )
#' 
#' # sample from the marginal posterior corresponding to documents 5 and 6
#' d5 <- posterior(
#'   x = m,
#'   matrix = "theta",
#'   which = c(5, 6),
#'   times = 100
#' )
#' }
#' @export
posterior <- function(x, ...) UseMethod("posterior")

#' Posterior method for tidylda
#' @rdname posterior
#' @export
posterior.tidylda <- function(
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
  
  # construct Dirichlet parameter(s) based on "matrix" argument
  if (matrix[1] == "theta") {
    
    # get proper alpha
    alpha <- format_alpha(
      x$alpha, 
      k = nrow(x$beta)
    )
    
    # extract dirichlet parameters for theta
    mat <- x$counts$Cd[which, ]
    
    if (length(which) == 1) {
      mat <- matrix(mat, nrow = 1)
    }
    
    dir_par <- t(mat) + alpha$alpha
    
    colnames(dir_par) <- rownames(x$theta)[which]
    rownames(dir_par) <- colnames(x$theta)
    
    
    
  } else {
    # get proper eta
    eta <- format_eta(
      x$eta, 
      k = nrow(x$beta), 
      Nv = ncol(x$beta)
    )
    
    # extract dirichlet parameters for beta
    dir_par <- x$counts$Cv[which, ] + eta$eta[which, ]
    
    if (length(which) == 1) {
      dir_par <- matrix(dir_par, nrow = 1)
    }
    
    
    rownames(dir_par) <- rownames(x$beta)[which]
    colnames(dir_par) <- colnames(x$beta)
    
    dir_par <- t(dir_par)
  }
  
  # sample
  result <- generate_sample(
    dir_par = dir_par,
    matrix = matrix,
    times = times
  )
  
  result
}

#' Generate a sample of LDA posteriors
#' @keywords internal
#' @description
#'   Helper function called by both posterior.tidylda and predict.tidylda to
#'   generate samples from the posterior.
#' @param dir_par matrix of Dirichlet hyperparameters, one column per
#' @param matrix character of "theta" or "beta", indicating which posterior
#'   matrix \code{dir_par}'s columns are from.
#' @param times Integer, number of samples to draw.
#' @return Returns a tibble with one row per parameter per sample.
generate_sample <- function(
  dir_par,
  matrix,
  times
) {
  
  result <- lapply(
    X = as.data.frame(dir_par),
    FUN = function(y) {
      samp <- gtools::rdirichlet(n = times, alpha = y)
      
      samp <- as.data.frame(samp, stringsAsFactors = FALSE)
      
      colnames(samp) <- rownames(dir_par)
      
      out <- as.data.frame(t(samp))
      
      colnames(out) <- 1:ncol(out)
      
      out$idx1 <- rownames(dir_par)
      
      out <- tidyr::pivot_longer(
        out,
        -idx1,
        names_to = "sample"
      )
      
      out
    }
  )
  
  
  # prepare and return result so it's tidy
  for (j in seq_along(result)) {
    result[[j]]$idx2 <- colnames(dir_par)[j]
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
