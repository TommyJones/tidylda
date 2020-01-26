#' Update methods for topic models
#' @description \code{update} updates a previously-trained topic model based
#' on new data or continues training a model with its original data. Useful for 
#' updates or transfer learning.
#' @param object An existing trained topic model
#' @param ... Additional arguments to the call
#' @export
update <- function(object, ...) UseMethod("update")

#' Update a Latent Dirichlet Allocation topic model
#' @description Update an LDA model using collapsed Gibbs sampling. 
#' @param object a fitted object of class \code{tidylda_model}.
#' @param dtm A document term matrix or term co-occurrence matrix of class dgCMatrix.
#' @param iterations Integer number of iterations for the Gibbs sampler to run. 
#' @param burnin Integer number of burnin iterations. If \code{burnin} is greater than -1,
#'        the resulting "phi" and "theta" matrices are an average over all iterations
#'        greater than \code{burnin}.
#' @param optimize_alpha Logical. Do you want to optimize alpha every iteration?
#'        Defaults to \code{FALSE}. See 'details' of documentation for
#'        \code{\link[textmineR]{FitLdaModel}}for more information.
#' @param calc_likelihood Logical. Do you want to calculate the log likelihood every iteration?
#'        Useful for assessing convergence. Defaults to \code{FALSE}. 
#' @param calc_r2 Logical. Do you want to calculate R-squared after the model is trained?
#'        Defaults to \code{FALSE}. This calls \code{\link[textmineR]{CalcTopicModelR2}}.
#' @param return_data Logical. Do you want \code{dtm} returned as part of the model object?
#' @param additional_k Integer number of topics to add, defaults to 0.
#' @param phi_as_prior Logical. Do you want to replace \code{beta} with \code{phi}
#'        from the previous model as the prior for words over topics?
#' @param ... Other arguments to be passed to \code{\link[furrr]{future_map}}
#' @return Returns an S3 object of class c("tidylda_model"). 
#' @details 
#' prior + counts vs. counts only. Vocab alignment + uniform prior over new words. 
#'          Adding additional topics. works best with significant vocab overlap
#' @export
#' @examples 
#' \dontrun{
#' # load a document term matrix
#' data(nih_sample_dtm, package = "textmineR")
#' 
#' d1 <- nih_sample_dtm[1:50,]
#' 
#' d2 <- nih_sample_dtm[51:100,]
#' 
#' # fit a model
#' m <- fit_tidylda(d1, k = 10, 
#'                   iterations = 200, burnin = 175)
#' 
#' # update an existing model by adding documents
#' m2 <- update(object = m,
#'              dtm = rbind(d1, d2),
#'              iterations = 200,
#'              burnin = 175)
#'              
#' # use an old model as a prior for a new model
#' m3 <- update(object = m,
#'              dtm = d2, # new documents only
#'              phi_as_prior = TRUE,
#'              iterations = 200,
#'              burnin = 175)
#'              
#' # add topics while updating a model by adding documents
#' m4 <- update(object = m,
#'              dtm = rbind(d1, d2),
#'              additional_k = 3,
#'              iterations = 200,
#'              burnin = 175)
#'              
#' 
#' }
update.tidylda_model <- function(object, dtm, iterations = NULL, burnin = -1, 
                                   optimize_alpha = FALSE, calc_likelihood = FALSE, 
                                   calc_r2 = FALSE, return_data = FALSE, 
                                   additional_k = 0, phi_as_prior = FALSE, ...) {
  
  # first, get the call for reproducibility
  mc <- match.call()
  
  ### Check inputs are of correct dimensionality ----
  
  # object of correct class?
  if (class(object) != "tidylda_model")
    stop("object must be of class tidylda_model")
  
  # iterations and burnin acceptable?
  if (burnin >= iterations) {
    stop("burnin must be less than iterations")
  }
  
  # dtm of the correct format?
  if (! "dgCMatrix" %in% class(dtm)) {
    message("'dtm' is not of class dgCMatrix, attempting to convert...")
    
    dtm <- try(methods::as(dtm, "dgCMatrix", strict = TRUE)) # requires Matrix in namespace
    
    if (! "dgCMatrix" %in% class(dtm))
      stop("conversion failed. Please pass an object of class dgCMatrix for dtm")
  }
  
  # is k formatted correctly?
  if (! is.numeric(additional_k)){
    stop("additional_k must be an integer >= 0")
  } else if (additional_k < 0) {
    stop("additional_k must be an integer >= 0")
  }
  
  additional_k <- floor(additional_k) # in case somebody is cheeky and passes a decimal
  
  # iterations?
  if (is.null(iterations))
    stop("You must specify number of iterations")
  
  if (! is.logical(calc_likelihood))
    stop("calc_likelihood must be TRUE or FALSE")
  
  if (! is.logical(phi_as_prior))
    stop("phi_as_prior must be TRUE or FALSE")
  
  
  ### Pull out objects used for update ----
  
  # format of beta
  if (phi_as_prior) {
    
    beta <- list(beta = object$phi,
                 beta_class = "matrix")
    
    # re-scale so that beta has the same magnitude of the old beta
    
    if (is.matrix(object$beta)) {
      
      beta$beta <- beta$beta * rowSums(object$beta)
      
    } else if (is.vector(object$beta)) {
      
      beta$beta <- beta$beta * sum(object$beta)
      
    } else if (length(object$beta) == 1) {
      
      beta$beta <- beta$beta * (object$beta * ncol(object$phi))
      
    } else { # this case shouldn't happen
      
      stop("object$beta must be a numeric scalar, a numeric vector of length 
         'ncol(object$phi)', or a numeric matrix with 'nrow(object$phi)' rows 
         and 'ncol(object$phi)' columns with no missing  values and at least 
         one non-zero value.")
      
    }
    
  } else {
    
    beta <- format_beta(beta = object$beta, k = nrow(object$phi), Nv = ncol(object$phi))
    
  }
  
  dimnames(beta$beta) <- dimnames(object$phi)
  
  # phi_initial and theta_initial
  phi_initial <- object$phi
  
  theta_initial <- predict.tidylda_model(object = object,
                                           newdata = dtm,
                                           method = "dot",
                                           ...)
  
  # pull out alpha
  alpha <- format_alpha(alpha = object$alpha,
                        k = nrow(object$phi))
  
  ### Vocabulary alignment and new topic (if any) alignment ----
  
  # align vocab in intelligent way for adding new vocab
  v_diff <- setdiff(colnames(dtm), colnames(phi_initial))
  
  m_add <- matrix(0, nrow = nrow(phi_initial), ncol = length(v_diff))
  
  colnames(m_add) <- v_diff
  
  beta$beta <- cbind(beta$beta, m_add+ median(beta$beta)) # uniform prior over new words
  
  beta$beta <- beta$beta[, colnames(dtm)]
  
  phi_initial <- cbind(phi_initial, m_add + median(phi_initial)) 
  
  phi_initial <- phi_initial[, colnames(dtm)] / rowSums(phi_initial[, colnames(dtm)])
  
  # add topics to beta and phi_initial
  # prior for topics inherets from beta, specifically colMeans(beta)
  # basically means that new topics are an average of old topics. If you used
  # a scalar or vector for object$beta, then prior for new topics will be 
  # identical to prior for old topics. If object$beta was a matrix where rows
  # were not identical (i.e. you seeded specific topics), then your new topics
  # will have a prior that is the average of all old topics.
  m_add <- matrix(0, 
                  nrow = additional_k, 
                  ncol = ncol(beta$beta))
  
  m_add <- t(t(m_add) + colMeans(beta$beta)) 
  
  beta$beta <- rbind(beta$beta, m_add) # add new topics to beta
  
  phi_initial <- rbind(phi_initial, m_add / rowSums(m_add)) # new topics to phi
  
  # add topics to alpha and theta_initial
  # prior for new topics is uniform, similar to beta, it's the median of alpha
  # adding new topics to theta_inital is a little more complicated. We take the
  # median of each row of theta_initial, add that to the new topics and then
  # reweight so each row still sums to 1.
  alpha$alpha <- c(alpha$alpha, rep(median(alpha$alpha), additional_k)) # uniform prior for new topics
  
  m_add <- apply(theta_initial, 1, function(x){
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
  counts <- initialize_topic_counts(dtm = dtm, 
                                    k = nrow(phi_initial),
                                    alpha = alpha$alpha,
                                    beta = beta$beta,
                                    phi_initial = phi_initial,
                                    theta_initial = theta_initial,
                                    freeze_topics = FALSE, # false because this is an update
                                    ...) 
  
  ### run C++ gibbs sampler ----
  lda <- fit_lda_c(docs = counts$docs,
                   Nk = nrow(phi_initial),
                   alpha = alpha$alpha,
                   beta = beta$beta,
                   Cd = counts$Cd,
                   Cv = counts$Cv,
                   Ck = counts$Ck,
                   Zd = counts$Zd,
                   Phi = counts$Cv, # this is actually ignored as freeze_topics = FALSE 
                   iterations = iterations,
                   burnin = burnin,
                   freeze_topics = FALSE, # this stays FALSE for updates 
                   calc_likelihood = calc_likelihood, 
                   optimize_alpha = optimize_alpha) 
  

  
  ### Format output correctly ----
  result <- format_raw_lda(lda = lda, dtm = dtm, burnin = burnin, 
                           is_prediction = FALSE, 
                           alpha = alpha, beta = beta, 
                           optimize_alpha = optimize_alpha, calc_r2 = calc_r2, 
                           calc_likelihood = calc_likelihood, 
                           call = mc, ...)
  
  ### return the result ----
  result
  
}

