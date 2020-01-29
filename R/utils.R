################################################################################
# Functions in this file are internal to tidylda
################################################################################

#' Format \code{beta} For Input into \code{fit_lda_c}
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{beta} but the C++ Gibbs
#'   sampler in \code{\link[tidylda]{fit_lda_c}} only takes it one way. This function does the 
#'   approprate formatting. It also returns errors if the user input a malformatted
#'   \code{beta}.
#' @param beta the prior for words over topics. Can be a numeric scalar, numeric 
#'   vector, or numeric matrix.
#' @param k the number of topics.
#' @param Nv the total size of the vocabulary as inherited from \code{ncol(dtm)}
#'   in \code{\link[tidylda]{tidylda}}.
#' @return 
#'   Returns a list with two elements: \code{beta} and \code{beta_class}. 
#'   \code{beta} is the post-formatted versionof \code{beta} in the form of a 
#'   \code{k} by \code{Nv} numeric matrix. \code{beta_class} is a character 
#'   denoting whether or not the user-supplied \code{beta} was a "scalar", 
#'   "vector", or "matrix".
format_beta <- function(beta, k, Nv) {
  
  if (! is.numeric(beta) | sum(is.na(beta)) > 0 | sum(beta == 0) == length(beta))
    stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
  
  if (length(beta) == 1) { # if beta is a scalar
    
    beta <- matrix(beta, nrow = k, ncol = Nv)
    
    beta_class <- "scalar"
    
  } else if (is.vector(beta)){ # if beta is a vector
    
    if (length(beta) != Nv) # if you didn't specify this vector right
      stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    
    # otherwise let's carry on...
    # make beta a matrix to format for C++ funciton
    beta <- t(beta + matrix(0, nrow = length(beta), ncol = k))
    
    beta_class <- "vector"
    
  } else if (is.matrix(beta)) { # if beta is a matrix
    
    beta_class <- "matrix"
    
  } else { # if beta is of an un supported data type
    
    stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    
  }
  
  
  list(beta = beta,
       beta_class = beta_class)
}

#' Format \code{alpha} For Input into \code{fit_lda_c}
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{alpha} but the C++ Gibbs
#'   sampler in \code{\link[tidylda]{fit_lda_c}} only takes it one way. This function does the 
#'   approprate formatting. It also returns errors if the user input a malformatted
#'   \code{alpha}.
#' @param alpha the prior for topics over documents. Can be a numeric scalar or 
#'   numeric vector.
#' @param k the number of topics.
#' @return 
#'   Returns a list with two elements: \code{alpha} and \code{alpha_class}. 
#'   \code{alpha} is the post-formatted version of \code{alpha} in the form of a 
#'   \code{k}-length numeric vector. \code{alpha_class} is a character 
#'   denoting whether or not the user-supplied \code{alpha} was a "scalar" or
#'   "vector".
format_alpha <- function(alpha, k) {
  
  if (! is.numeric(alpha) | sum(is.na(alpha)) > 0 | sum(alpha == 0) == length(alpha))
    stop("alpha must be a numeric scalar or vector of length 'k' with no missing 
          values and at least one non-zero value")
  
  if (length(alpha) == 1 & is.numeric(alpha)) {
    alpha <- numeric(k) + alpha
    
    alpha_class <- "scalar"
    
  } else if (length(alpha) != k | ! is.vector(alpha)){
    
    stop("alpha must be a numeric scalar or numeric vector of length 'k'")
    
  } else {
    
    alpha_class <- "vector"
    
  }
  
  list(alpha = alpha,
       alpha_class = alpha_class)
  
}

#' Initialize topic counts for gibbs sampling
#' @keywords internal
#' @description
#'   Implementing seeded (or guided) LDA models and transfer learning means that
#'   we can't initialize topics with a uniform-random start. This function prepares
#'   data and then calls a C++ function, \code{\link[tidylda]{create_lexicon}}, that runs a single
#'   Gibbs iteration to populate topic counts (and other objects) used during the
#'   main Gibbs sampling run of \code{\link[tidylda]{fit_lda_c}}. In the event that 
#'   you aren't using fancy seeding or transfer learning, this makes a random
#'   initialization by sampling from Dirichlet distributions paramaterized by 
#'   priors \code{alpha} and \code{beta}.
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @param k the number of topics 
#' @param alpha the numeric vector prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_alpha}}
#' @param beta the numeric matrix prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_beta}}
#' @param phi_initial if specified, a numeric matrix for the probability of tokens
#'   in topics. Must be specified for predictions or updates as called by
#'   \code{\link[tidylda]{predict.tidylda}} or \code{\link[tidylda]{update.tidylda}}
#'   respectively.
#' @param theta_initial if specified, a numeric matrix for the probability of
#'   topics in documents. Must be specified for updates as called by 
#'   \code{\link[tidylda]{update.tidylda}}
#' @param freeze_topics if \code{TRUE} does not update counts of tokens in topics.
#'   This is \code{TRUE} for predictions.
#' @param ... other items to be passed to \code{\link[furrr]{future_map}}
#' @return 
#'   Returns a list with 5 elements: \code{docs}, \code{Zd}, \code{Cd}, \code{Cv}, 
#'   and \code{Ck}. All of these are used by \code{\link[tidylda]{fit_lda_c}}.
#'   
#'   \code{docs} is a list with one element per document. Each element is a vector
#'   of integers of length \code{sum(dtm[j,])} for the j-th document. The integer
#'   entries correspond to the zero-index column of the \code{dtm}.
#'   
#'   \code{Zd} is a list of similar format as \code{docs}. The difference is that
#'   the integer values correspond to the zero-index for topics.
#'   
#'   \code{Cd} is a matrix of integers denoting how many times each topic has
#'   been sampled in each document.
#'   
#'   \code{Cv} is similar to \code{Cd} but it counts how many times each topic
#'   has been sampled for each token.
#'   
#'   \code{Ck} is an integer vector denoting how many times each topic has been
#'   sampled overall. 
#' @note
#'   All of \code{Cd}, \code{Cv}, and \code{Ck} should be derivable by summing
#'   over Zd in various ways.
initialize_topic_counts <- function(dtm, k, alpha, beta, phi_initial = NULL, 
                                    theta_initial = NULL, freeze_topics = FALSE, 
                                    ...) {
  
  # check inputs
  
  # initialize phi if not already specified
  # this phi is used to sample topics for inital counts in the C++ function
  if (is.null(phi_initial)) {
    # phi_initial <- gtools::rdirichlet(n = k, alpha = beta)
    
    phi_initial <- apply(beta, 1, function(x){
      gtools::rdirichlet(n = 1, alpha = x)
    })
    
    phi_initial <- t(phi_initial)
  }
  
  # initialize theta if not already specified
  # if not specified (e.g. if this is a new model) make a matrix by sampling
  # from alpha. 
  if (is.null(theta_initial)) {
    
    theta_initial <- gtools::rdirichlet(n = nrow(dtm), alpha = alpha)
    
  }
  
  # initalize Cd by sampling from theta_initial. 
  # for asymmetric alpha, encodes more/less probable topics
  # we don't need to initialize Cv because we can use the probabilities in phi, 
  # along with our sampled Cd to do a single Gibbs iteration to populate all three
  # of Cd, Ck, and Cv
  cd_sampler <- function(size, prob){
    stats::rmultinom(n = 1, size = size, prob = prob)
  }
  
  # below makes use of furrr::future_map to allow for parallel processing
  # in cases where we have more than 3000 documents
  batches <- seq(1, nrow(dtm), by = 3000)
  
  # future::plan(multiprocess) 
  
  iterator <- furrr::future_map(.x = batches, .f = function(b){
    rows <- b:min(b + 2999, nrow(dtm))
    
    if (length(rows) > 1) {
      
      size <- Matrix::rowSums(dtm[rows, ])
      
      prob <- as.list(data.frame(t(theta_initial[rows, ])))
      
    } else {
      size <- sum(dtm[rows, ])
      
      prob <- theta_initial[rows, ]
    }
    
    list(size = size, prob = prob)
  }, ...)
  
  Cd_start <- furrr::future_map(.x = iterator,
                                .f = function(x) {
                                  out <- mapply(FUN = cd_sampler,
                                                size = x$size,
                                                prob = x$prob)
                                  
                                  t(out)
                                }, ...)
  
  Cd_start <- do.call(rbind, Cd_start)

  
  # Initialize objects with that single Gibbs iteration mentioned above
  # if we have more than 3000 documents, do it in parallel with furrr::future_map

  batches <- seq(1, nrow(dtm), by = 3000)
  
  lexicon <- furrr::future_map(.x = batches, 
                               .f = function(b){
    
    rows <- b:min(b + 2999, nrow(dtm))
    
    # if statement to handle single observations
    if (length(rows) == 1) {
      cd_tmp <- Cd_start
      
      dtm_tmp <- dtm
    } else {
      cd_tmp <- Cd_start[rows, ]
      
      dtm_tmp <- dtm[rows, ]
    }
    
    l <- create_lexicon(Cd = cd_tmp,
                        Phi = phi_initial,
                        dtm = dtm_tmp,
                        alpha = alpha,
                        freeze_topics = freeze_topics)
    
  }, ...)
  
  # combine 
  Zd <- Reduce("c", lapply(lexicon, function(l) l$Zd))
  
  docs <- Reduce("c", lapply(lexicon, function(l) l$docs))
  
  Cv <- Reduce("+", lapply(lexicon, function(l) l$Cv))
  
  Ck <- Reduce("+", lapply(lexicon, function(l) l$Ck))
  
  Cd <- do.call(rbind, lapply(lexicon, function(l) l$Cd))
  
  out <- list(docs = docs,
              Zd = Zd,
              Cd = Cd,
              Cv = Cv,
              Ck = Ck)
  
  out
  
}

#' Summarize a topic model consistently across methods/functions
#' @keywords internal
#' @description
#'   Summarizes topics in a model. Called by \code{\link[tidylda]{tidylda}}
#'   and \code{\link[tidylda]{update.tidylda}} and used to augment
#'   \code{\link[tidylda]{print.tidylda}}.
#' @param theta numeric matrix whose rows represent P(topic|document)
#' @param phi numeric matrix whose rows represent P(token|topic)
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @return 
#'   Returns a \code{\link[tibble]{tibble}} with the following columns: 
#'   \code{topic} is the integer row number of \code{phi}. 
#'   \code{prevalence} is the frequency of each topic throughout the corpus it 
#'     was trained on normalized so that it sums to 100. 
#'   \code{coherence} makes a call to \code{\link[textmineR]{CalcProbCoherence}}
#'     using the default 5 most-probable terms in each topic.
#'   \code{top_terms} displays the top 3 most-probable terms in each topic.   
#' @note
#'   \code{prevalence} should be proportional to P(topic). It is calculated by 
#'   weighting on document length. So, topics prevalent in longer documents get 
#'   more weight than topics prevalent in shorter documents. It is calculated
#'   by 
#'   
#'   \code{prevalence <- rowSums(dtm) * theta \%>\% colSums()}
#'     
#'   \code{prevalence <- (prevalence * 100) \%>\% round(3)}
#'   
#'   An alternative calculation (not implemented here) might have been
#'    
#'   \code{prevalence <- colSums(dtm) * t(phi) \%>\% colSums()}
#'   
#'   \code{prevalence <- (prevalence * 100) \%>\% round(3)}
summarize_topics <- function(theta, phi, dtm){
  
  # probabilistic coherence with default M = 5
  coherence <- textmineR::CalcProbCoherence(phi = phi, dtm = dtm)
  
  # prevalence of each topic, weighted by terms
  prevalence <- Matrix::rowSums(dtm) * theta
  
  prevalence <- colSums(prevalence) / sum(prevalence)
  
  prevalence <- round(prevalence * 100, 2)
  
  # top 3 terms
  top_terms <- t(textmineR::GetTopTerms(phi, 3))
  
  top_terms <- apply(top_terms, 1, function(x){
    paste(c(x, "..."), collapse = ", ")
  })
  
  # combine into a summary
  summary <- data.frame(topic = as.numeric(rownames(phi)),
                        prevalence = prevalence,
                        coherence = coherence,
                        top_terms = top_terms,
                        stringsAsFactors = FALSE)
  
  summary <- tibble::as_tibble(summary)
  
  summary
  
}

#' Format the outputs of \code{\link[tidylda]{fit_lda_c}} consistently
#' @keywords internal
#' @description
#'   Since all three of \code{\link[tidylda]{tidylda}}, 
#'   \code{\link[tidylda]{update.tidylda}}, and 
#'   \code{\link[tidylda]{predict.tidylda}} call \code{\link[tidylda]{fit_lda_c}},
#'   we need a way to format the resulting posteriors and other user-facing
#'   objects consistently. This function does that.
#' @param lda list output of \code{\link[tidylda]{fit_lda_c}}
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}
#' @param burnin integer number of burnin iterations. 
#' @param is_prediction is this for a prediction (as opposed to initial fitting,
#'   or update)? Defaults to \code{FALSE}
#' @param alpha output of \code{\link[tidylda]{format_alpha}}
#' @param beta output of \code{\link[tidylda]{format_beta}}
#' @param optimize_alpha did you optimize \code{alpha} when making a call to 
#'   \code{\link[tidylda]{fit_lda_c}}?  If \code{is_prediction = TRUE}, this 
#'   argument is ignored.
#' @param calc_r2 did the user want to calculate R-squared when calculating the
#'   the model? If \code{is_prediction = TRUE}, this argument is ignored.
#' @param calc_likelihood did you calculate the log likelihood when making a call
#'   to \code{\link[tidylda]{fit_lda_c}}?  If \code{is_prediction = TRUE}, this 
#'   rgument is ignored.
#' @param call the result of calling \code{\link[base]{match.call}} at the top of 
#'   \code{\link[tidylda]{tidylda}}.
#' @return
#'   Returns an S3 object of class \code{tidylda} with the following slots:
#'   
#'   \code{phi} is a numeric matrix whose rows are the posterior estimates 
#'     of P(token|topic)
#'   
#'   \code{theta} is a numeric matrix  whose rows are the posterior estimates of
#'     P(topic|document)
#'   
#'   \code{gamma} is a numeric matrix whose rows are the posterior estimates of
#'     P(topic|token). This is calculated by making a call to 
#'     \code{link[textmineR]{CalcGamma}} which uses Bayes's rule to calculate
#'     \code{gamma} from \code{phi}, \code{theta}, and P(document) (which is 
#'     proportional to \code{Matrix::rowSums(dtm)}).
#'   
#'   \code{alpha} is the prior for topics over documents. If \code{optimize_alpha}
#'     is \code{FALSE}, \code{alpha} is what the user passed when calling 
#'     \code{\link[tidylda]{tidylda}}. If \code{optimize_alpha} is \code{TRUE}, 
#'     \code{alpha} is a numeric vector returned in the \code{alpha} slot from a 
#'     call to \code{\link[tidylda]{fit_lda_c}}.
#'   
#'   \code{beta} is the prior for tokens over topics. This is what the user passed 
#'     when calling \code{\link[tidylda]{tidylda}}.
#'   
#'   \code{summary} is the result of a call to \code{\link[tidylda]{summarize_topics}}
#'   
#'   \code{call} is the result of \code{\link[base]{match.call}} called at the top
#'     of \code{\link[tidylda]{tidylda}}
#'   
#'   \code{log_likelihood} is a \code{\link[tibble]{tibble}} whose columns are 
#'     the iteration and log likelihood at that iteration. This slot is only populated
#'     if \code{calc_likelihood = TRUE}
#'   
#'   \code{r2} is a numeric scalar resuting from a call to 
#'     \code{\link[textmineR]{CalcTopicModelR2}}. This slot only populated if
#'     \code{calc_r2 = TRUE}
#' @note
#'   In general, the arguments of this function should be what the user passed
#'   when calling \code{\link[tidylda]{tidylda}}. 
#'   
#'   \code{burnin} is used only to determine whether or not burn in iterations 
#'   were used when fitting the model. If \code{burnin > -1} then posteriors 
#'   are calculated using \code{lda$Cd_mean} and \code{lda$Cv_mean} respectively. 
#'   Otherwise, posteriors are calculated using \code{lda$Cd_mean} and 
#'   \code{lda$Cv_mean}.
#'   
#'   The class of \code{call} isn't checked. It's just passed through to the 
#'   object returned by this function. Might be useful if you are using this
#'   function for troubleshooting or something. 
new_tidylda <- function(lda, dtm, burnin, is_prediction = FALSE, 
                           alpha = NULL, beta = NULL, 
                           optimize_alpha = NULL, calc_r2 = NULL, 
                           calc_likelihood = NULL, call = NULL,
                           ...) {
  
  ### format theta ----
  if (burnin > -1) {
    
    theta <- t(t(lda$Cd_mean) + lda$alpha)
    
  } else {
    
    theta <- t(t(lda$Cd) + lda$alpha)
    
  }
  
  theta <- theta / rowSums(theta)
  
  theta[is.na(theta)] <- 0 # just in case of a numeric issue
  
  colnames(theta) <- seq_len(ncol(theta))
  
  rownames(theta) <- rownames(dtm)
  
  ### format phi and all the rest ----
  
  if (! is_prediction) {
    ### format posteriors correctly ----
    if (burnin > -1) { # if you used burnin iterations use Cd_mean etc.
      
      phi <- lda$Cv_mean + lda$beta
      
      
    } else { # if you didn't use burnin use standard counts (Cd etc.)
      
      phi <- lda$Cv + lda$beta
      
    }
    
    phi <- phi / rowSums(phi)
    
    phi[is.na(phi)] <- 0 # just in case of a numeric issue
    
    colnames(phi) <- colnames(dtm)
    
    rownames(phi) <- colnames(theta)
    
    
    ### collect the results ----
    
    # gamma
    gamma <- textmineR::CalcGamma(phi = phi, theta = theta, 
                                  p_docs = Matrix::rowSums(dtm))
    
    # beta
    colnames(lda$beta) <- colnames(phi)
    
    if (beta$beta_class == "scalar") {
      
      beta_out <- lda$beta[1, 1]
      
    } else if (beta$beta_class == "vector") {
      
      beta_out <- lda$beta[1, ]
      
    } else if (beta$beta_class == "matrix") {
      
      beta_out <- lda$beta
      
    } else { # this should be impossible, but science is hard and I am dumb.
      beta_out <- lda$beta
      
      message("something went wrong with 'beta'. This isn't your fault. Please 
            contact Tommy at jones.thos.w[at]gmail.com and tell him you got this
            error when you ran 'tidylda'.")
    }
    
    # alpha
    
    if (alpha$alpha_class == "scalar" & !optimize_alpha) {
      
      alpha_out <- lda$alpha[1]
      
    } else if (alpha$alpha_class == "vector" | optimize_alpha) {
      
      alpha_out <- lda$alpha
      
      names(alpha_out) <- rownames(phi)
      
    }
    
    # resulting object
    result <- list(phi = phi,
                   theta = theta,
                   gamma = gamma,
                   alpha = alpha_out,
                   beta = beta_out,
                   log_likelihood = as_tibble(data.frame(iteration = lda$log_likelihood[1,],
                                                         log_likelihood = lda$log_likelihood[2, ]))
    ) # add other things here if necessary
    
    class(result) <- "tidylda"
    
    ### calculate and add other things ---
    
    result$summary <- summarize_topics(phi = result$phi, theta = result$theta,
                                       dtm = dtm)
    
    # get arguments for auditiability
    # result$other_call_args <- list(iterations = iterations, 
    #                                burnin = burnin,
    #                                optimize_alpha = optimize_alpha)
    
    # goodness of fit
    if (calc_r2) {
      result$r2 <- textmineR::CalcTopicModelR2(dtm = dtm, 
                                               phi = result$phi, 
                                               theta = result$theta, ...)
    }
    
    # call
    result$call <- call
    
    # a little cleanup here
    if (! calc_likelihood) {
      result$log_likelihood <- NULL
    }
    
  }
  
  
  ### return the final result ----
  if (is_prediction) {
    
    return(theta)
    
  } else {
    
    return(result)
    
  }
  
}
