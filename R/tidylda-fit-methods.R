#' Fit a Latent Dirichlet Allocation topic model
#' @description Fit a Latent Dirichlet Allocation topic model using warpLDA, a
#'   Metropolis-Hastings sampler.
#' @param data A document term matrix or term co-occurrence matrix. The preferred
#'   class is a \code{\link[Matrix]{dgCMatrix-class}}. However there is support
#'   for any \code{\link[Matrix]{Matrix-class}} object as well as several other
#'   commonly-used classes such as \code{\link[base]{matrix}},
#'   \code{\link[quanteda]{dfm}}, \code{\link[tm]{DocumentTermMatrix}}, and
#'   \code{\link[slam]{simple_triplet_matrix}}
#' @param k Integer number of topics.
#' @param iterations Integer number of sampling iterations to run.
#' @param burnin Integer number of burnin iterations. If \code{burnin} is greater than -1,
#'        the resulting "beta" and "theta" matrices are an average over all iterations
#'        greater than \code{burnin}.
#' @param alpha Numeric scalar or vector of length \code{k}. This is the prior
#'        for topics over documents.
#' @param eta Numeric scalar, numeric vector of length \code{ncol(data)},
#'        or numeric matrix with \code{k} rows and \code{ncol(data)} columns.
#'        This is the prior for words over topics.
#' @param optimize_alpha Deprecated as of version 0.1.0 and ignored. Accepted
#'        so that existing calls keep working; passing \code{TRUE} warns once
#'        per session. See 'details' below.
#' @param calc_likelihood Logical. Do you want to calculate the log likelihood every iteration?
#'        Useful for assessing convergence. Defaults to \code{TRUE}.
#' @param calc_r2 Logical. Do you want to calculate R-squared after the model is trained?
#'        Defaults to \code{FALSE}. See \code{\link[tidylda]{calc_lda_r2}}.
#' @param threads Number of parallel threads, defaults to 1. See Details, below.
#' @param return_data Logical. Do you want \code{data} returned as part of the model object?
#' @param verbose Logical. Do you want to print a progress bar out to the console?
#'        Defaults to \code{TRUE}.
#' @param ... Additional arguments, currently unused
#' @return Returns an S3 object of class \code{tidylda}. See \code{\link[tidylda]{new_tidylda}}.
#' @details Fitting uses **warpLDA** (Chen et al.,
#'   \url{https://arxiv.org/abs/1510.08628}), a Metropolis-Hastings sampler that
#'   alternates document-ordered and word-ordered passes over the corpus so that
#'   each pass touches only a small, cache-resident working set. It replaces the
#'   collapsed Gibbs sampler used through version 0.0.7, and is written in Rcpp
#'   and parallelized with RcppThread. Some implementation notes follow:
#'
#'   Topic-token and topic-document assignments are not initialized based on a
#'   uniform-random sampling, as is common. Instead, topic-token probabilities
#'   (i.e. \code{beta}) are initialized by sampling from a Dirichlet distribution
#'   with \code{eta} as its parameter. The same is done for topic-document
#'   probabilities (i.e. \code{theta}) using \code{alpha}. Then an internal
#'   function is called (\code{\link[tidylda]{initialize_topic_counts}}) to run
#'   a single sampling iteration to initialize assignments of tokens to topics
#'   and topics to documents.
#'
#'   When you use burn-in iterations (i.e. \code{burnin = TRUE}), the resulting
#'   \code{beta} and \code{theta} matrices are calculated by averaging over every
#'   iteration after the specified  number of burn-in iterations. If you do not
#'   use burn-in iterations, then the matrices are calculated from the last run
#'   only. Ideally, you'd burn in every iteration before convergence, then average
#'   over the chain after its converged (and thus every observation is independent).
#'
#'   \code{optimize_alpha} is deprecated as of version 0.1.0 and is ignored. It
#'   rescaled \code{alpha} by topic size each iteration, standing in for
#'   fixed-point estimation that was never written. \code{alpha} is now fixed
#'   for the whole run. The argument will be removed in a future release.
#'
#'   \strong{Two log likelihood columns, and they answer different questions.}
#'   When \code{calc_likelihood = TRUE}, the \code{log_likelihood} slot has a
#'   column of each.
#'
#'   \code{log_likelihood} is \eqn{P(tokens | \theta, \beta)}: the probability
#'   of the observed tokens under the current estimates of \code{theta} and
#'   \code{beta}. It is a plug-in quantity, so it improves monotonically as
#'   topics are added and \strong{cannot be used to choose \code{k}}.
#'
#'   \code{log_joint} is \eqn{P(tokens, topics | \alpha, \eta)}, the collapsed
#'   joint, with \code{theta} and \code{beta} analytically integrated out. This
#'   is the quantity most of the LDA literature reports. Integrating out the
#'   parameters leaves a discrete distribution, so unlike a density it is always
#'   negative, and it carries an implicit penalty for model complexity. It is
#'   also the sampler's own target, which makes it the more informative of the
#'   two for judging convergence.
#'
#'   Both condition on the current assignment of tokens to topics, so both are
#'   \emph{within-model} diagnostics. Neither is a valid basis for comparing
#'   models to each other; that requires \eqn{P(tokens | \alpha, \eta)} with the
#'   topic assignments marginalized out, which is intractable and needs the
#'   held-out estimators of Wallach et al. (2009).
#'
#'   Both are evaluated every tenth iteration by default rather than every
#'   iteration. This is not thinning: the chain advances every iteration and
#'   every post-burn-in iteration still contributes to the posterior means. Only
#'   these diagnostics, which feed nothing the sampler uses, are computed less
#'   often. Pass \code{likelihood_every = 1} to recover a value per iteration.
#'   Both accept \code{eta} as a scalar, a vector, or a matrix, here and in
#'   \code{\link[tidylda]{refit.tidylda}} and
#'   \code{\link[tidylda]{predict.tidylda}}.
#'
#'   \code{threads} sets the number of worker threads. Results are identical at
#'   any thread count, so it trades wall clock for cores and nothing else; a
#'   model fitted under \code{set.seed()} is reproducible whether it was fitted
#'   on one thread or twenty. It defaults to 1, so that \code{tidylda} never
#'   takes cores it was not asked for -- which matters when fitting many models
#'   inside your own parallel loop.
#'
#' @examples
#' # load some data
#' data(nih_sample_dtm)
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
#' # predict on held-out documents using Metropolis-Hastings "fold in"
#' p1 <- predict(m, nih_sample_dtm[21:100, ],
#'   method = "mh",
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
  eta = 0.05,
  optimize_alpha = FALSE, 
  calc_likelihood = TRUE,
  calc_r2 = FALSE, 
  threads = 1,
  return_data = FALSE,
  verbose = TRUE,
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
    eta = eta,
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
#'   resulting topic model. Most arguments to this function are documented in
#'   \code{\link[tidylda]{tidylda}}; the two below are specific to the warpLDA
#'   engine and reach this function through \code{tidylda}'s \code{...}.
#' @param likelihood_every integer. Evaluate the log likelihood every n-th
#'   iteration. Defaults to 10. The log likelihood costs
#'   \eqn{O(nnz \cdot K + VK)} and would otherwise dominate an
#'   \eqn{O(VK + N)} sampler. This is not thinning: the chain advances every
#'   iteration and every post-burn-in iteration still contributes to the count
#'   sums. Only the read-only diagnostic runs less often.
#' @param mh_steps integer. Number of Metropolis-Hastings proposals made per
#'   token per pass. Defaults to 1. Larger values mix further per iteration at
#'   a cost of \code{mh_steps * 2} bytes per token.
#' @return Returns a \code{tidylda} S3 object as documented in \code{\link[tidylda]{new_tidylda}}.
tidylda_bridge <- function(
  data, 
  k, 
  iterations, 
  burnin, 
  alpha, 
  eta,
  optimize_alpha,
  calc_likelihood,
  calc_r2,
  threads,
  return_data,
  verbose,
  mc,
  likelihood_every = 10,
  mh_steps = 1,
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

  # alpha and eta?
  alpha <- format_alpha(alpha = alpha, k = k)

  eta <- format_eta(eta = eta, k = k, Nv = ncol(dtm))

  if (!is.numeric(likelihood_every) || length(likelihood_every) != 1 ||
      likelihood_every < 1) {
    stop("likelihood_every must be a single integer >= 1")
  }

  if (!is.numeric(mh_steps) || length(mh_steps) != 1 || mh_steps < 1) {
    stop("mh_steps must be a single integer >= 1")
  }

  # D7 removed optimize_alpha: it rescaled alpha proportional to Ck each
  # iteration as a placeholder for fixed-point estimation that was never
  # implemented. The argument is retained so existing calls keep working, but it
  # no longer does anything, and silently ignoring it would be worse than saying
  # so. Fixing alpha is also what lets D19's alias table be built just once.
  # Warned once per session rather than per call: this is informational, and a
  # loop of fits should not produce a wall of identical warnings. Documented as
  # deprecated as of 0.1.0; the argument itself goes in a later release.
  if (isTRUE(optimize_alpha)) {
    rlang::warn(
      paste0(
        "optimize_alpha is deprecated as of tidylda 0.1.0 and is ignored.\n",
        "  It rescaled alpha by topic size as a stand-in for fixed-point\n",
        "  estimation that was never written; alpha is now fixed for the run.\n",
        "  The argument will be removed in a future release."
      ),
      .frequency = "once",
      .frequency_id = "tidylda_optimize_alpha_removed"
    )
  }

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
  
  # The old warning here -- that fewer than 100 documents per thread gives "a
  # poor fit" -- described an abandoned batched implementation, where the
  # partition genuinely changed the model. Under D12 the warpLDA engine seeds
  # every work item from its own index, so results are identical at any thread
  # count and the only thing `threads` buys or costs is wall clock. Nothing to
  # warn about.
  #
  # The document-count bound is gone for the same reason, and because it was
  # measuring the wrong thing: the word pass parallelizes over the vocabulary,
  # not over documents.

  ### format inputs ----


  # initialize counts
  counts <- initialize_topic_counts(
    dtm = dtm,
    k = k,
    alpha = alpha$alpha, 
    eta = eta$eta,
    threads = threads
  )

  ### run the C++ sampler ----
  lda <- fit_lda_warp(
    dtm_in = dtm,
    Cd_start = counts$Cd_start,
    alpha_in = alpha$alpha,
    eta_in = as.matrix(eta$eta), # 1 x 1 when scalar; the engine detects it (D20)
    iterations = iterations,
    burnin = burnin,
    calc_likelihood = calc_likelihood,
    Beta_in = counts$beta_initial,
    freeze_topics = FALSE,
    likelihood_every = as.integer(likelihood_every),
    mh_steps = as.integer(mh_steps),
    threads = as.integer(threads),
    verbose = verbose
  )

  ### format the output ----

  result <- new_tidylda(
    lda = lda, 
    dtm = dtm, 
    burnin = burnin,
    is_prediction = FALSE,
    alpha = alpha, 
    eta = eta,
    optimize_alpha = optimize_alpha, 
    calc_r2 = calc_r2,
    calc_likelihood = calc_likelihood,
    call = mc,
    threads = threads
  )

  ### return data if desired ----
  if (return_data) {
    result$data <- dtm
  }
  
  ### return the result ----
  result
}
