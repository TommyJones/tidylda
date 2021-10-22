#' Update a Latent Dirichlet Allocation topic model
#' @description Update an LDA model using collapsed Gibbs sampling.
#' @param object a fitted object of class \code{tidylda}.
#' @param new_data A document term matrix or term co-occurrence matrix of class dgCMatrix.
#' @param iterations Integer number of iterations for the Gibbs sampler to run.
#' @param burnin Integer number of burnin iterations. If \code{burnin} is greater than -1,
#'        the resulting "beta" and "theta" matrices are an average over all iterations
#'        greater than \code{burnin}.
#' @param prior_weight Numeric, 0 or greater or \code{NA}. The weight of the 
#'        \code{beta} as a prior from the base model. See Details, below.
#' @param additional_k Integer number of topics to add, defaults to 0.
#' @param additional_eta_sum Numeric magnitude of prior for additional topics.
#'        Ignored if \code{additional_k} is 0. Defaults to 250.
#' @param optimize_alpha Logical. Do you want to optimize alpha every iteration?
#'        Defaults to \code{FALSE}. See 'details' of documentation for
#'        \code{\link[textmineR]{FitLdaModel}}for more information.
#' @param calc_likelihood Logical. Do you want to calculate the log likelihood every iteration?
#'        Useful for assessing convergence. Defaults to \code{FALSE}.
#' @param calc_r2 Logical. Do you want to calculate R-squared after the model is trained?
#'        Defaults to \code{FALSE}. This calls \code{\link[textmineR]{CalcTopicModelR2}}.
#' @param return_data Logical. Do you want \code{new_data} returned as part of the model object?
#' @param threads Number of parallel threads, defaults to 1.
#' @param verbose Logical. Do you want to print a progress bar out to the console?
#'        Defaults to \code{TRUE}.
#' @param ... Additional arguments, currently unused
#' @return Returns an S3 object of class c("tidylda").
#' @details
#'   \code{refit} allows you to (a) update the probabilities (i.e. weights) of
#'   a previously-fit model with new data or additional iterations and (b) optionally
#'   use \code{beta} of a previously-fit LDA topic model as the \code{eta} prior
#'   for the new model. This is tuned by setting \code{beta_as_prior = FALSE} or
#'   \code{beta_as_prior = TRUE} respectively.
#'   
#'   \code{prior_weight} tunes how strong the base model is represented in the prior.
#'   If \code{prior_weight = 1}, then the tokens from the base model's training data
#'   have the same relative weight as tokens in \code{new_data}. In other words,
#'   it is like just adding training data. If \code{prior_weight} is less than 1,
#'   then tokens in \code{new_data} are given more weight. If \code{prior_weight}
#'   is greater than 1, then the tokens from the base model's training data are
#'   given more weight.
#'
#'   If \code{prior_weight} is \code{NA}, then the new \code{eta} is equal to 
#'   \code{eta} from the old model, with new tokens folded in. 
#'   (For handling of new tokens, see below.) Effectively, this just controls
#'   how the sampler initializes (described below), but does not give prior
#'   weight to the base model.
#'
#'   Instead of initializing token-topic assignments in the manner for new
#'   models (see \code{\link[tidylda]{tidylda}}), the update initializes in 2
#'   steps:
#'
#'   First, topic-document probabilities (i.e. \code{theta}) are obtained by a
#'   call to \code{\link[tidylda]{predict.tidylda}} using \code{method = "dot"}
#'   for the documents in \code{new_data}. Next, both \code{beta} and \code{theta} are
#'   passed to an internal function, \code{\link[tidylda]{initialize_topic_counts}},
#'   which assigns topics to tokens in a manner approximately proportional to 
#'   the posteriors and executes a single Gibbs iteration.
#'
#'   \code{refit} handles the addition of new vocabulary by adding a flat prior
#'   over new tokens. Specifically, each entry in the new prior is equal to the
#'   10th percentile of \code{eta} from the old model. The resulting model will
#'   have the total vocabulary of the old model plus any new vocabulary tokens.
#'   In other words, after running \code{refit.tidylda} \code{ncol(beta) >= ncol(new_data)}
#'   where \code{beta} is from the new model and \code{new_data} is the additional data.
#'
#'   You can add additional topics by setting the \code{additional_k} parameter
#'   to an integer greater than zero. New entries to \code{alpha} have a flat
#'   prior equal to the median value of \code{alpha} in the old model. (Note that
#'   if \code{alpha} itself is a flat prior, i.e. scalar, then the new topics have
#'   the same value for their prior.) New entries to \code{eta} have a shape 
#'   from the average of all previous topics in \code{eta} and scaled by
#'   \code{additional_eta_sum}.
#' @note
#'  Updates are, as of this writing, are almost-surely useful but their behaviors
#'  have not been optimized or well-studied. Caveat emptor!
#' @export
#' @examples
#' \donttest{
#' # load a document term matrix
#' data(nih_sample_dtm, package = "textmineR")
#'
#' d1 <- nih_sample_dtm[1:50, ]
#'
#' d2 <- nih_sample_dtm[51:100, ]
#'
#' # fit a model
#' m <- tidylda(d1,
#'   k = 10,
#'   iterations = 200, burnin = 175
#' )
#'
#' # update an existing model by adding documents using old model as prior
#' m2 <- refit(
#'   object = m,
#'   new_data = rbind(d1, d2),
#'   iterations = 200,
#'   burnin = 175,
#'   prior_weight = 1
#' )
#'
#' # use an old model to initialize new model and not use old model as prior
#' m3 <- refit(
#'   object = m,
#'   new_data = d2, # new documents only
#'   iterations = 200,
#'   burnin = 175,
#'   prior_weight = NA
#' )
#'
#' # add topics while updating a model by adding documents
#' m4 <- refit(
#'   object = m,
#'   new_data = rbind(d1, d2),
#'   additional_k = 3,
#'   iterations = 200,
#'   burnin = 175
#' )
#' }
refit.tidylda <- function(
  object, 
  new_data, 
  iterations = NULL, 
  burnin = -1,
  prior_weight = 1,
  additional_k = 0, 
  additional_eta_sum = 250,
  optimize_alpha = FALSE, 
  calc_likelihood = FALSE,
  calc_r2 = FALSE, 
  return_data = FALSE,
  threads = 1,
  verbose = TRUE,
  ...
) {
  
  # first, get the call for reproducibility
  mc <- match.call()
  
  ### Check inputs are of correct dimensionality ----
  
  # iterations and burnin acceptable?
  if (burnin >= iterations) {
    stop("burnin must be less than iterations")
  }
  
  # Ensure dtm is of class dgCMatrix
  dtm <- convert_dtm(dtm = new_data)
  
  # Ensure dtm has column names
  if (is.null(colnames(dtm))) {
    stop("new_data must have names for tokens. Did you pass a matrix without colnames?")
  }
  
  # is prior weight formatted correctly?
  if (! (is.na(prior_weight) | is.numeric(prior_weight))) {
    stop("prior_weight must be numeric greater than 0 or NA.")
  }
  
  if (! is.na(prior_weight)) {
    if (prior_weight <= 0) {
      stop("prior_weight must be numeric greater than 0 or NA.")
    }
  }
  
  # is k formatted correctly?
  if (!is.numeric(additional_k)) {
    stop("additional_k must be an integer >= 0")
  } else if (additional_k < 0) {
    stop("additional_k must be an integer >= 0")
  }
  
  additional_k <- floor(additional_k) # in case somebody is cheeky and passes a decimal
  
  # iterations?
  if (is.null(iterations)) {
    stop("You must specify number of iterations")
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
  
  if (! is.logical(optimize_alpha)) {
    stop("optimize_alpha must be logical")
  }
  
  ### Pull out objects used for update ----
  
  # get base model eta into a matrix for downstream formatting
  eta <- format_eta(
    eta = object$eta, 
    k = nrow(object$beta), 
    Nv = ncol(object$beta)
  )  
  
  
  # if necessary, re-scale so that new eta has the weight prescribed by prior-weight
  if (! is.na(prior_weight)) {
    w_star <- rowSums(object$counts$Cv) + rowSums(eta$eta)
    
    eta$eta <- prior_weight * w_star * object$beta
    
    eta$eta_class <- "matrix" # always a matrix for refits using eta as prior
  }
  
  dimnames(eta$eta) <- dimnames(object$beta)
  
  # beta_initial and theta_initial
  beta_initial <- object$beta
  
  theta_initial <- predict.tidylda(
    object = object,
    new_data = dtm,
    method = "dot",
    no_common_tokens = "uniform",
    threads = threads
  )
  
  # pull out alpha
  alpha <- format_alpha(
    alpha = object$alpha,
    k = nrow(object$beta)
  )
  
  ### Vocabulary alignment and new topic (if any) alignment ----
  
  # align vocab in intelligent way for adding new vocab
  total_vocabulary <- union(colnames(dtm), colnames(beta_initial))
  
  add_to_dtm <- setdiff(total_vocabulary, colnames(dtm))
  
  add_to_model <- setdiff(total_vocabulary, colnames(beta_initial))
  
  m_add_to_dtm <- matrix(0, nrow = nrow(dtm), ncol = length(add_to_dtm))
  
  colnames(m_add_to_dtm) <- add_to_dtm
  
  m_add_to_model <- matrix(0, nrow = nrow(beta_initial), ncol = length(add_to_model))
  
  colnames(m_add_to_model) <- add_to_model
  
  dtm <- cbind(dtm, m_add_to_dtm)
  

  # uniform prior over new words
  eta$eta <- cbind(eta$eta, m_add_to_model + stats::quantile(eta$eta, 0.1))

  eta$eta <- eta$eta[, colnames(dtm)]

  beta_initial <- cbind(beta_initial, m_add_to_model + stats::quantile(beta_initial, 0.1))

  beta_initial <- beta_initial[, colnames(dtm)] / rowSums(beta_initial[, colnames(dtm)])
  
  
  # add topics to eta and beta_initial
  # prior for topics inherets from eta, specifically colMeans(eta)
  # basically means that new topics are an average of old topics. If you used
  # a scalar or vector for object$eta, then prior for new topics will be
  # identical to prior for old topics. If object$eta was a matrix where rows
  # were not identical (i.e. you seeded specific topics), then your new topics
  # will have a prior that is the average of all old topics.
  m_add <- matrix(0,
                  nrow = additional_k,
                  ncol = ncol(eta$eta)
  )
  
  m_add <- t(t(m_add) + colMeans(eta$eta))
  
  m_add <- m_add / rowSums(m_add) * additional_eta_sum
  
  eta$eta <- rbind(eta$eta, m_add) # add new topics to eta
  
  beta_initial <- rbind(beta_initial, m_add / rowSums(m_add)) # new topics to beta
  
  # add topics to alpha and theta_initial
  # prior for new topics is uniform, similar to eta, it's the median of alpha
  # adding new topics to theta_inital is a little more complicated. We take the
  # median of each row of theta_initial, add that to the new topics and then
  # reweight so each row still sums to 1.
  alpha$alpha <- c(alpha$alpha, rep(median(alpha$alpha), additional_k)) # uniform prior for new topics
  
  m_add <- apply(theta_initial, 1, function(x) {
    rep(median(x), additional_k)
  })
  
  # handle cases on what m_add could be
  if (is.matrix(m_add)) { # if we add more than one topic
    
    m_add <- t(m_add)
    
    colnames(m_add) <- (max(as.numeric(colnames(theta_initial))) + 1):
      (max(as.numeric(colnames(theta_initial))) + additional_k)
    
    theta_initial <- cbind(theta_initial, m_add)
    
    theta_initial <- theta_initial / rowSums(theta_initial)
  } else if (length(m_add) == 0) { # we add no topics and get nothing back
    
    # do nothing, actually
  } else { # we add only one topic and get a vector back
    
    theta_initial <- cbind(theta_initial, m_add)
    
    theta_initial <- theta_initial / rowSums(theta_initial)
  }
  
  
  ### get initial counts to feed to gibbs sampler ----
  counts <- initialize_topic_counts(
    dtm = dtm,
    k = nrow(beta_initial),
    alpha = alpha$alpha,
    eta = eta$eta,
    beta_initial = beta_initial,
    theta_initial = theta_initial,
    freeze_topics = FALSE, # false because this is an update
    threads = threads
  )
  
  ### run C++ gibbs sampler ----
  lda <- fit_lda_c(
    Docs = counts$Docs,
    Zd_in = counts$Zd,
    Cd_in = counts$Cd,
    Cv_in = counts$Cv,
    Ck_in = counts$Ck,
    alpha_in = alpha$alpha,
    eta_in = eta$eta,
    iterations = iterations,
    burnin = burnin,
    optimize_alpha = optimize_alpha,
    calc_likelihood = calc_likelihood,
    Beta_in = object$beta, # ignored for updates as freeze_topics = FALSE
    freeze_topics = FALSE,
    threads = threads,
    verbose = verbose
  )
  
  ### Format output correctly ----
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
  
  ### return the result ----
  result
}
