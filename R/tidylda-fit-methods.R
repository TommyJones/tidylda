#' Fit a Latent Dirichlet Allocation topic model
#' @description Fit a Latent Dirichlet Allocation topic model using collapsed Gibbs sampling.
#' @param data A document term matrix or term co-occurrence matrix. The preferred
#'   class is a \code{\link[Matrix]{dgCMatrix-class}}. However there is support
#'   for any \code{\link[Matrix]{Matrix-class}} object as well as several other
#'   commonly-used classes such as \code{\link[base]{matrix}},
#'   \code{\link[quanteda]{dfm}}, \code{\link[tm]{DocumentTermMatrix}}, and
#'   \code{\link[slam]{simple_triplet_matrix}}
#' @param k Integer number of topics.
#' @param iterations Integer number of iterations for the Gibbs sampler to run.
#' @param burnin Integer number of burnin iterations. If \code{burnin} is greater than -1,
#'        the resulting "phi" and "theta" matrices are an average over all iterations
#'        greater than \code{burnin}.
#' @param alpha Numeric scalar or vector of length \code{k}. This is the prior
#'        for topics over documents.
#' @param beta Numeric scalar, numeric vector of length \code{ncol(data)},
#'        or numeric matrix with \code{k} rows and \code{ncol(data)} columns.
#'        This is the prior for words over topics.
#' @param optimize_alpha Logical. Do you want to optimize alpha every iteration?
#'        Defaults to \code{FALSE}. See 'details' below for more information.
#' @param calc_likelihood Logical. Do you want to calculate the log likelihood every iteration?
#'        Useful for assessing convergence. Defaults to \code{TRUE}.
#' @param calc_r2 Logical. Do you want to calculate R-squared after the model is trained?
#'        Defaults to \code{FALSE}. This calls \code{\link[textmineR]{CalcTopicModelR2}}.
#' @param threads Number of parallel threads, defaults to 1. See Details, below.
#' @param return_data Logical. Do you want \code{data} returned as part of the model object?
#' @param verbose Logical. Do you want to print a progress bar out to the console?
#'        Defaults to \code{FALSE}.
#' @param ... Additional arguments, currently unused
#' @return Returns an S3 object of class \code{tidylda}.
#' @details This function calls a collapsed Gibbs sampler for Latent Dirichlet Allocation
#'   written using the excellent Rcpp package. Some implementation notes follow:
#'
#'   Topic-token and topic-document assignments are not initialized based on a
#'   uniform-random sampling, as is common. Instead, topic-token probabilities
#'   (i.e. \code{phi}) are initialized by sampling from a Dirichlet distribution
#'   with \code{beta} as its parameter. The same is done for topic-document
#'   probabilities (i.e. \code{theta}) using \code{alpha}. Then an internal
#'   function is called (\code{\link[tidylda]{initialize_topic_counts}}) to run
#'   a single Gibbs iteration to initialize assignments of tokens to topics and
#'   topics to documents.
#'
#'   When you use burn-in iterations (i.e. \code{burnin = TRUE}), the resulting
#'   \code{phi} and \code{theta} matrices are calculated by averaging over every
#'   iteration after the specified  number of burn-in iterations. If you do not
#'   use burn-in iterations, then the matrices are calculated from the last run
#'   only. Ideally, you'd burn in every iteration before convergence, then average
#'   over the chain after its converged (and thus every observation is independent).
#'
#'   If you set \code{optimize_alpha} to \code{TRUE}, then each element of \code{alpha}
#'   is proportional to the number of times each topic has be sampled that iteration
#'   averaged with the value of \code{alpha} from the previous iteration. This lets
#'   you start with a symmetric \code{alpha} and drift into an asymmetric one.
#'   However, (a) this probably means that convergence will take longer to happen
#'   or convergence may not happen at all. And (b) I make no guarantees that doing this
#'   will give you any benefit or that it won't hurt your model. Caveat emptor!
#'
#'   The log likelihood calculation is the same that can be found on page 9 of
#'   \url{https://arxiv.org/pdf/1510.08628.pdf}. The only difference is that the
#'   version in \code{\link[tidylda]{tidylda}} allows \code{beta} to be a
#'   vector or matrix. (Vector used in this function, matrix used for model
#'   updates in \code{\link[tidylda]{refit.tidylda}}. At present, the
#'   log likelihood function appears to be ok for assessing convergence. i.e. It
#'   has the right shape. However, it is, as of this writing, returning positive
#'   numbers, rather than the expected negative numbers. Looking into that, but
#'   in the meantime caveat emptor once again.
#'   
#'   Parallelism, is executed using threading at the C++ level using the
#'   \code{\link[RcppThread]{RcppThread}} package. For Gibbs sampling, parallelism
#'   is implemented according to
#'   \url{http://www.jmlr.org/papers/volume10/newman09a/newman09a.pdf?wptouch_preview_theme=enabled}
#'
#' @examples
#' # load some data
#' data(nih_sample_dtm, package = "textmineR")
#'
#' # fit a model
#' set.seed(12345)
#' m <- tidylda(
#'   data = nih_sample_dtm[1:20, ], k = 5,
#'   iterations = 200, burnin = 175
#' )
#'
#' str(m)
#'
#' # predict on held-out documents using gibbs sampling "fold in"
#' p1 <- predict(m, nih_sample_dtm[21:100, ],
#'   method = "gibbs",
#'   iterations = 200, burnin = 175
#' )
#'
#' # predict on held-out documents using the dot product method
#' p2 <- predict(m, nih_sample_dtm[21:100, ], method = "dot")
#'
#' # compare the methods
#' barplot(rbind(p1[1, ], p2[1, ]), beside = TRUE, col = c("red", "blue"))
#' @export
tidylda <- function(
  data, 
  k, 
  iterations = NULL, 
  burnin = -1, 
  alpha = 0.1, 
  beta = 0.05,
  optimize_alpha = FALSE, 
  calc_likelihood = TRUE,
  calc_r2 = FALSE, 
  threads = 1,
  return_data = FALSE,
  verbose = FALSE,
  ...
) {

  # not using methods for now as I think this is cleaner
  # UseMethod("tidylda")

  # first, get the call for reproducibility
  mc <- match.call()


  tidylda_bridge(
    data = data,
    k = k,
    iterations = iterations,
    burnin = burnin,
    alpha = alpha,
    beta = beta,
    optimize_alpha = optimize_alpha,
    calc_likelihood = calc_likelihood,
    calc_r2 = calc_r2,
    threads = threads,
    return_data = return_data,
    verbose = verbose,
    mc,
    ...
  )
}


#' Bridge function for fitting \code{tidylda} topic models
#' @keywords internal
#' @description
#'   Takes in arguments from various \code{tidylda} S3 methods and fits the
#'   resulting topic model. The arguments to this function are documented in
#'   \code{\link[tidylda]{tidylda}}.
tidylda_bridge <- function(
  data, 
  k, 
  iterations, 
  burnin, 
  alpha, 
  beta,
  optimize_alpha, 
  calc_likelihood, 
  calc_r2,
  threads,
  return_data, 
  verbose,
  mc,
  ...
) {

  ### check validity of inputs ----

  # iterations and burnin acceptable?
  if (burnin >= iterations) {
    stop("burnin must be less than iterations")
  }

  # Ensure dtm is of class dgCMatrix
  dtm <- convert_dtm(dtm = data)
  
  # Ensure dtm has column names
  if (is.null(colnames(dtm))) {
    stop("data must have names for tokens. Did you pass a matrix without colnames?")
  }

  # is k formatted correctly?
  if (k < 2) {
    stop("k must be 2 or greater")
  }

  if (!is.numeric(k)) {
    stop("k must be an integer")
  }

  k <- floor(k) # in case somebody is cheeky and passes a decimal

  # iterations?
  if (is.null(iterations)) {
    stop("You must specify number of iterations")
  }

  # alpha and beta?
  alpha <- format_alpha(alpha = alpha, k = k)

  beta <- format_beta(beta = beta, k = k, Nv = ncol(dtm))

  # are you being logical
  if (!is.logical(calc_r2)) {
    stop("calc_r2 must be logical")
  }

  if (!is.logical(calc_likelihood)) {
    stop("calc_likelihood must be logical")
  }

  if (!is.logical(return_data)) {
    stop("return_data must be logical")
  }

  if (!is.logical(optimize_alpha)) {
    stop("optimize_alpha must be logical")
  }

  if (!is.logical(verbose)) {
    stop("verbose must be logical")
  }
  
  # check on threads
  if (threads > 1)
    threads <- as.integer(max(floor(threads), 1)) # prevent any decimal inputs
  
  if (threads > nrow(dtm)) {
    stop("User-supplied threads argument greater than number of documents.\n",
         "  Recommend setting threads such that nrow(dtm) / threads > 100,\n",
         "  More documents on each processor is better.")
  }
  
  if ((nrow(dtm) / threads < 100) & (threads > 1)) {
    warning("  nrow(dtm) / threads < 100.\n",
            "  If each processor has fewer than 100 documents, resulting model is likely\n",
            "  to be a poor fit. More documents on each processor is better.")
  }

  ### format inputs ----


  # initialize counts
  counts <- initialize_topic_counts(
    dtm = dtm, 
    k = 10,
    alpha = alpha$alpha, 
    beta = beta$beta,
    threads = threads
  )

  # divide into batches to enable parallel execution of the Gibbs sampler
  

  ### run C++ gibbs sampler ----
  lda <- fit_lda_c(
    Docs = counts$Docs,
    Zd_in = counts$Zd,
    Cd_in = counts$Cd,
    Cv_in = counts$Cv,
    Ck_in = counts$Ck,
    alpha_in = alpha$alpha,
    beta_in = beta$beta,
    iterations = iterations,
    burnin = burnin,
    optimize_alpha = optimize_alpha,
    calc_likelihood = calc_likelihood,
    Phi_in = counts$Cv, # this is actually ignored as freeze_topics = FALSE for initial fitting
    freeze_topics = FALSE, # this stays FALSE for initial fitting
    threads = threads,
    verbose = verbose
  )

  ### format the output ----

  result <- format_raw_lda_outputs(
    lda = lda, 
    dtm = dtm, 
    burnin = burnin,
    is_prediction = FALSE,
    alpha = alpha, 
    beta = beta,
    optimize_alpha = optimize_alpha, 
    calc_r2 = calc_r2,
    calc_likelihood = calc_likelihood,
    call = mc,
    threads = threads
  )

  ### return the result ----
  result
}
