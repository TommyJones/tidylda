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
#' \href{https://web.archive.org/web/2020id_/http://www.arbylon.net/publications/text-est.pdf}{Archived copy}
#' (arbylon.net no longer resolves; this is the Internet Archive's copy.)
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
    #
    # Cv is words-by-topics (D17), so a topic is a COLUMN -- and the result is
    # already words-by-topics, which is what generate_sample() wants. The
    # trailing t() this used to need is gone. counts_cv() transposes a model
    # saved by an earlier version.
    eta_mat <- eta_matrix(eta, nrow(x$beta), ncol(x$beta))

    # as.matrix() because Cv is a dgCMatrix and eta is a base matrix: Matrix
    # promotes sparse + dense to a DENSE Matrix (dgeMatrix), not to a sparse one.
    # Everything downstream --- the matrix() reshape just below, the dimnames
    # assignments, and generate_sample()'s as.data.frame() --- wants a base
    # matrix, and a dgeMatrix satisfies none of them.
    dir_par <- as.matrix(counts_cv(x)[, which] + t(eta_mat[which, , drop = FALSE]))
    
    if (length(which) == 1) {
      dir_par <- matrix(dir_par, ncol = 1)
    }
    
    colnames(dir_par) <- rownames(x$beta)[which]
    rownames(dir_par) <- colnames(x$beta)
  }
  
  # sample
  # length(which) * times draws of nrow(dir_par) parameters, one row each. At
  # k = 1000, V = 1e6 and the default times = 100 that is 1e11 rows.
  check_result_size(
    n_rows = as.numeric(nrow(dir_par)) * ncol(dir_par) * times,
    n_cols = 4,
    what = paste0("posterior(matrix = \"", matrix[1], "\")"),
    suggestion = paste0(
      "Lower `times`, or pass fewer values in `which` and reduce each result ",
      "before requesting the next -- taking slices only bounds memory if you ",
      "do not keep them all."
    )
  )

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
  
  # The long frame is built straight from the draw matrix. This used to go
  # rdirichlet -> as.data.frame -> t() -> as.data.frame -> pivot_longer, four
  # copies of every block: 12.2 MB of draws became a 582 MB peak.
  #
  # rdirichlet returns times by n, so as.vector() walks it column-major --- all
  # samples of the first parameter, then all samples of the second. That is
  # exactly the order pivot_longer produced, which is why idx1 repeats `each =
  # times` and sample cycles fastest. `sample` is built as character because the
  # caller below converts it with as.numeric(), as it did when pivot_longer
  # supplied the column names.
  #
  # RNG order is unchanged: one rdirichlet(n = times, .) per column of dir_par,
  # in column order.
  result <- lapply(
    seq_len(ncol(dir_par)),
    function(j) {
      samp <- gtools::rdirichlet(n = times, alpha = dir_par[, j])
      
      data.frame(
        idx1   = rep(rownames(dir_par), each = times),
        sample = rep(as.character(seq_len(times)), times = nrow(dir_par)),
        value  = as.vector(samp),
        idx2   = colnames(dir_par)[j],
        stringsAsFactors = FALSE
      )
    }
  )
  
  result <- dplyr::bind_rows(result)
  
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
