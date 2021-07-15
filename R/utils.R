################################################################################
# Functions in this file are internal to tidylda
################################################################################

#' Convert various things to a \code{dgCMatrix} to work with various functions
#' and methods
#' @keywords internal
#' @description
#'   Presently, \code{tidylda} makes heavy usage of the \code{dgCMatrix} class.
#'   However, a user may have created a DTM (or TCM) in one of several classes.
#'   Since data could be in several formats, this function converts them to a
#'   \code{dgCMatrix} before passing them along.
#' @param dtm the data you want to convert
#' @return an object of class \code{dgCMatrix}
convert_dtm <- function(dtm) {
  if (inherits(dtm, c("Matrix", "matrix"))) {
    out <- methods::as(dtm, "dgCMatrix", strict = TRUE)
  } else if (inherits(dtm, "simple_triplet_matrix")) {
    out <- Matrix::sparseMatrix(
      i = dtm$i,
      j = dtm$j,
      x = dtm$v,
      dims = c(dtm$nrow, dtm$ncol),
      dimnames = list(
        rownames = dtm$dimnames$Docs,
        colnames = dtm$dimnames$Terms
      )
    )
  } else if (inherits(dtm, "numeric")) {
    if (is.null(names(dtm))) {
      stop(
        "it looks like dtm (or new_data if you called 'predict') is a numeric ",
        "vector without names. Did you mean to pass a single document? If so, ",
        "it needs a names attribute to index tokens"
      )
    }
    
    vocab <- names(dtm)
    
    out <- Matrix::Matrix(dtm, nrow = 1, sparse = TRUE)
    
    colnames(out) <- vocab
    
    rownames(out) <- 1
  } else {
    stop(
      "dtm (or new_data if you called 'predict') cannot be converted to dgCMatrix. Supported classes are ",
      "c('Matrix', 'matrix', 'simple_triplet_matrix', 'dfm', 'DocumentTermMatrix'), ",
      "However, I see class(dtm) = ", class(dtm)
    )
  }
  
  out
}

#' Format \code{eta} For Input into \code{fit_lda_c}
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{eta} but the C++ Gibbs
#'   sampler in \code{\link[tidylda]{fit_lda_c}} only takes it one way. This function does the
#'   appropriate formatting. It also returns errors if the user input a malformatted
#'   \code{eta}.
#' @param eta the prior for words over topics. Can be a numeric scalar, numeric
#'   vector, or numeric matrix.
#' @param k the number of topics.
#' @param Nv the total size of the vocabulary as inherited from \code{ncol(dtm)}
#'   in \code{\link[tidylda]{tidylda}}.
#' @return
#'   Returns a list with two elements: \code{eta} and \code{eta_class}.
#'   \code{eta} is the post-formatted version of \code{eta} in the form of a
#'   \code{k} by \code{Nv} numeric matrix. \code{eta_class} is a character
#'   denoting whether or not the user-supplied \code{eta} was a "scalar",
#'   "vector", or "matrix".
format_eta <- function(eta, k, Nv) {
  if (!is.numeric(eta) | sum(is.na(eta)) > 0 | sum(eta == 0) == length(eta)) {
    stop("eta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
  }
  
  if (length(eta) == 1) { # if eta is a scalar
    
    eta <- matrix(eta, nrow = k, ncol = Nv)
    
    eta_class <- "scalar"
  } else if (is.vector(eta)) { # if eta is a vector
    
    if (length(eta) != Nv) { # if you didn't specify this vector right
      stop("eta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    }
    
    # otherwise let's carry on...
    # make eta a matrix to format for C++ funciton
    eta <- t(eta + matrix(0, nrow = length(eta), ncol = k))
    
    eta_class <- "vector"
  } else if (is.matrix(eta)) { # if eta is a matrix
    
    # check dims before moving on
    if (nrow(eta) != k | ncol(eta) != Nv) {
      stop(
        "If eta is a matrix, it must have the same number of rows as topics ",
        "and it must have the same number of columns (tokens) as your dtm. ",
        "But I see nrow(eta) = ", nrow(eta), " and k = ", k, ". I also see ",
        "ncol(eta) = ", ncol(eta), " but ncol(dtm) = ", Nv
      )
    }
    
    eta_class <- "matrix"
  }
  
  
  list(
    eta = eta,
    eta_class = eta_class
  )
}

#' Format \code{alpha} For Input into \code{fit_lda_c}
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{alpha} but the C++ Gibbs
#'   sampler in \code{\link[tidylda]{fit_lda_c}} only takes it one way. This function does the
#'   appropriate formatting. It also returns errors if the user input a malformatted
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
  if (!is.numeric(alpha) | sum(is.na(alpha)) > 0 | sum(alpha == 0) == length(alpha)) {
    stop("alpha must be a numeric scalar or vector of length 'k' with no missing 
          values and at least one non-zero value")
  }
  
  if (length(alpha) == 1 & is.numeric(alpha)) {
    alpha <- numeric(k) + alpha
    
    alpha_class <- "scalar"
  } else if (length(alpha) != k | !is.vector(alpha)) {
    stop("alpha must be a numeric scalar or numeric vector of length 'k'")
  } else {
    alpha_class <- "vector"
  }
  
  list(
    alpha = alpha,
    alpha_class = alpha_class
  )
}

#' Get Count Matrices from Beta or Theta (and Priors)
#' @keywords internal
#' @description
#'   This function is a core component of \code{\link[tidylda]{initialize_topic_counts}}.
#'   See details, below.
#' @param prob_matrix a numeric \code{beta} or \code{theta} matrix
#' @param prior_matrix a matrix of same dimension as \code{prob_matrix} whose 
#'   entries represent the relevant prior (\code{alpha} or \code{eta})
#' @param total_vector a vector of token counts of length \code{ncol(prob_matrix)}
#' @details 
#'   This function uses a probability matrix (theta or beta), its prior (alpha or
#'   eta, respectively), and a vector of counts to simulate what the the Cd or
#'   Cv matrix would be at the end of a Gibbs run that resulted in that probability
#'   matrix.
#'   
#'   For example, theta is calculated from a matrix of counts, Cd, and a prior,
#'   alpha. Specifically, the i,j entry of theta is given by
#'   
#'   \code{(Cd[i, j] + alpha[i, j]) / sum(Cd[, j] + alpha[, j])}
#'   
#'   Similarly, beta comes from
#'   
#'   \code{(Cv[i, j] + eta[i, j]) / sum(Cv[, j] + eta[, j])}
#'   
#'   (The above are written to be general with respect to alpha and eta being
#'   matrices. They could also be vectors or scalars.)
#'   
#'   So, this function uses the above formulas to try and reconstruct Cd or Cv
#'   from theta and alpha or beta and eta, respectively. As of this writing,
#'   this method is experimental. In the future, there will be a paper with
#'   more technical details cited here.
#'   
#'   The priors must be matrices for the purposes of the function. This is to
#'   support topic seeding and model updates. The former requires eta to be a 
#'   matrix. The latter may require eta to be a matrix. Here, alpha is also
#'   required to be a matrix for compatibility.
#'   
#'   All that said, for now \code{\link[tidylda]{initialize_topic_counts}} only
#'   uses this function to calculate Cd.
recover_counts_from_probs <- function(prob_matrix, prior_matrix, total_vector) {
  # prob_matrix is D X k
  # prior_matrix is D X k
  # total_vector is of length D
  
  # first, we have to get the denominator
  denom <- total_vector + (ncol(prior_matrix) * prior_matrix)
  
  # then, multiply probabilities by the denominator 
  count_matrix <- prob_matrix * denom # pointwise multiplication
  
  # subtract the prior to get what the counts *should* be
  count_matrix <- count_matrix - prior_matrix
  
  # reconcile the counts so that they're integers and line up to the right totals
  count_matrix <- apply(count_matrix, 1, function(x){
    
    tot <- sum(x)
    
    round_x <- round(x)
    
    remainder <- round(tot - sum(round_x))
    
    if (remainder == 0) {
      
      return(round_x)
      
    } else if (remainder > 0) { # we need to add some
      
      sample_prob <- x
      
      sample_prob[sample_prob < 0] <- 0
      
      # account for length of x
      # no need to sample 1 object, can do deterministically
      if (length(x) == 1) {
        idx <- 1
      } else {
        idx <- 
          try({
            sample(seq_along(x), remainder, replace = TRUE, prob = sample_prob)
          })
        
        if (inherits(x = idx, what = "try-error")) {
          stop("Something went wrong allocating counts.\n",
               "This is a low level error. Please contact the 'tidylda' maintainer.\n",
               "Error occured while adding some counts.\n",
               "remainder = ", remainder, "\n",
               "length(x) = ", length(x), "\n",
               "length(sample_prob) = ", length(sample_prob))
        }
        
      }
      
      
      round_x[idx] <- round_x[idx] + 1
      
      return(round_x)
      
    } else { # we need to take some away
      
      sample_prob <- x[round_x > 0]
      
      sample_prob[sample_prob < 0] <- 0
      
      sample_from <- seq_along(x)[round_x > 0]
      
      sample_size <- -1 * remainder
      
      # have to accommodate tokens that are out of bounds. 
      # if number of tokens is 1 then do it deterministically
      if (length(sample_from) == 1) {
        
        idx <- sample_from
        
      } else {
        
        idx <- 
          try({
            sample(x = sample_from, size = sample_size, replace = TRUE, prob = sample_prob)
          })
        
        if (inherits(x = idx, what = "try-error")) {
          stop("Something went wrong allocating counts.\n",
               "This is a low level error. Please contact the 'tidylda' maintainer.\n",
               "Error occured while subtracting some counts.\n",
               "remainder = ", remainder, "\n",
               "sample_size = ", sample_size, "\n",
               "length(sample_from) = ", length(sample_from), "\n",
               "length(sample_prob) = ", length(sample_prob))
        }
      }
      
      round_x[idx] <- round_x[idx] - 1
      
      return(round_x)
    }
    
  })
  
  count_matrix <- t(count_matrix)
  
  
  count_matrix
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
#'   initialization by sampling from Dirichlet distributions parameterized by
#'   priors \code{alpha} and \code{eta}.
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @param k the number of topics
#' @param alpha the numeric vector prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_alpha}}
#' @param eta the numeric matrix prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_eta}}
#' @param beta_initial if specified, a numeric matrix for the probability of tokens
#'   in topics. Must be specified for predictions or updates as called by
#'   \code{\link[tidylda]{predict.tidylda}} or \code{\link[tidylda]{refit.tidylda}}
#'   respectively.
#' @param theta_initial if specified, a numeric matrix for the probability of
#'   topics in documents. Must be specified for updates as called by
#'   \code{\link[tidylda]{refit.tidylda}}
#' @param freeze_topics if \code{TRUE} does not update counts of tokens in topics.
#'   This is \code{TRUE} for predictions.
#' @param threads number of parallel threads, currently unused
#' @param ... Additional arguments, currently unused
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
initialize_topic_counts <- function(
  dtm, 
  k, 
  alpha, 
  eta, 
  beta_initial = NULL,
  theta_initial = NULL, 
  freeze_topics = FALSE,
  threads = 1,
  ...
) {
  
  # check inputs
  
  if (! is.numeric(threads)) {
    stop("threads must be an integer 1 or greater")
  } else if (threads < 1) {
    stop("threads must be an integer 1 or greater")
  } else {
    threads = as.integer(threads) # ignore decimal inputs
  }
  
  # initialize beta if not already specified
  # this beta is used to sample topics for inital counts in the C++ function
  if (is.null(beta_initial)) {
    # beta_initial <- gtools::rdirichlet(n = k, alpha = eta)
    
    beta_initial <- apply(eta, 1, function(x) {
      gtools::rdirichlet(n = 1, alpha = x) + 
        .Machine$double.eps # avoid underflow
    })
    
    beta_initial <- t(beta_initial)
  }
  
  # initialize theta if not already specified
  # if not specified (e.g. if this is a new model) make a matrix by sampling
  # from alpha.
  if (is.null(theta_initial)) {
    theta_initial <- gtools::rdirichlet(n = nrow(dtm), alpha = alpha) + 
      .Machine$double.eps # avoid underflow
  }
  
  # initialize Cd by calling recover_counts_from_probs
  # we don't need to initialize Cv because we can use the probabilities in beta,
  # along with our sampled Cd to do a single Gibbs iteration to populate all three
  # of Cd, Ck, and Cv
  
  # format alpha to be a matrix to feed into recover_counts_from_probs
  alph <- matrix(0, nrow = ncol(theta_initial), ncol = nrow(theta_initial))
  
  alph <- alpha + alph
  
  alph <- t(alph)
  
  # get Cd itself
  # (note to future Tommy: consider a scalable version of this using future_map)
  Cd_start <- recover_counts_from_probs(
    prob_matrix = theta_initial,
    prior_matrix = alph,
    total_vector = Matrix::rowSums(dtm)
  )
  
  # Initialize objects with that single Gibbs iteration mentioned above
  # executed in parallel with RcppThread
  lexicon <- 
    create_lexicon(
      Cd_in = Cd_start,
      Beta_in = beta_initial,
      dtm_in = dtm, 
      alpha = alpha,
      freeze_topics = freeze_topics
    )
  
  lexicon
}

#' Summarize a topic model consistently across methods/functions
#' @keywords internal
#' @description
#'   Summarizes topics in a model. Called by \code{\link[tidylda]{tidylda}}
#'   and \code{\link[tidylda]{refit.tidylda}} and used to augment
#'   \code{\link[tidylda]{print.tidylda}}.
#' @param theta numeric matrix whose rows represent P(topic|document)
#' @param beta numeric matrix whose rows represent P(token|topic)
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @return
#'   Returns a \code{\link[tibble]{tibble}} with the following columns:
#'   \code{topic} is the integer row number of \code{beta}.
#'   \code{prevalence} is the frequency of each topic throughout the corpus it
#'     was trained on normalized so that it sums to 100.
#'   \code{coherence} makes a call to \code{\link[tidylda]{calc_prob_coherence}}
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
#'   \code{prevalence <- colSums(dtm) * t(beta) \%>\% colSums()}
#'
#'   \code{prevalence <- (prevalence * 100) \%>\% round(3)}
summarize_topics <- function(theta, beta, dtm) {
  
  # probabilistic coherence with default value for m
  coherence <- calc_prob_coherence(beta = beta, data = dtm)
  
  # prevalence of each topic, weighted by terms
  prevalence <- Matrix::rowSums(dtm) * theta
  
  prevalence <- colSums(prevalence) / sum(prevalence)
  
  prevalence <- round(prevalence * 100, 2)
  
  # top 3 terms
  top_terms <- apply(beta, 1, function(x) {
    names(x)[order(x, decreasing = TRUE)][1:3]
  })
  
  top_terms <- apply(top_terms, 2, function(x) {
    paste(c(x, "..."), collapse = ", ")
  })
  
  # combine into a summary
  summary <- data.frame(
    topic = as.numeric(rownames(beta)),
    prevalence = prevalence,
    coherence = coherence,
    top_terms = top_terms,
    stringsAsFactors = FALSE
  )
  
  summary <- tibble::as_tibble(summary)
  
  summary
}

#' Construct a new object of class \code{tidylda}
#' @keywords internal
#' @description
#'   Since all three of \code{\link[tidylda]{tidylda}},
#'   \code{\link[tidylda]{refit.tidylda}}, and
#'   \code{\link[tidylda]{predict.tidylda}} call \code{\link[tidylda]{fit_lda_c}},
#'   we need a way to format the resulting posteriors and other user-facing
#'   objects consistently. This function does that.
#' @param lda list output of \code{\link[tidylda]{fit_lda_c}}
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}
#' @param burnin integer number of burnin iterations.
#' @param is_prediction is this for a prediction (as opposed to initial fitting,
#'   or update)? Defaults to \code{FALSE}
#' @param alpha output of \code{\link[tidylda]{format_alpha}}
#' @param eta output of \code{\link[tidylda]{format_eta}}
#' @param optimize_alpha did you optimize \code{alpha} when making a call to
#'   \code{\link[tidylda]{fit_lda_c}}?  If \code{is_prediction = TRUE}, this
#'   argument is ignored.
#' @param calc_r2 did the user want to calculate R-squared when calculating the
#'   the model? If \code{is_prediction = TRUE}, this argument is ignored.
#' @param calc_likelihood did you calculate the log likelihood when making a call
#'   to \code{\link[tidylda]{fit_lda_c}}?  If \code{is_prediction = TRUE}, this
#'   argument is ignored.
#' @param call the result of calling \code{\link[base]{match.call}} at the top of
#'   \code{\link[tidylda]{tidylda}}.
#' @param threads number of parallel threads
#' @return
#'   Returns an S3 object of class \code{tidylda} with the following slots:
#'
#'   \code{beta} is a numeric matrix whose rows are the posterior estimates
#'     of P(token|topic)
#'
#'   \code{theta} is a numeric matrix  whose rows are the posterior estimates of
#'     P(topic|document)
#'
#'   \code{lambda} is a numeric matrix whose rows are the posterior estimates of
#'     P(topic|token), calculated using Bayes's rule.
#'     See \code{\link[tidylda]{calc_lambda}}.
#'
#'   \code{alpha} is the prior for topics over documents. If \code{optimize_alpha}
#'     is \code{FALSE}, \code{alpha} is what the user passed when calling
#'     \code{\link[tidylda]{tidylda}}. If \code{optimize_alpha} is \code{TRUE},
#'     \code{alpha} is a numeric vector returned in the \code{alpha} slot from a
#'     call to \code{\link[tidylda]{fit_lda_c}}.
#'
#'   \code{eta} is the prior for tokens over topics. This is what the user passed
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
#'   \code{r2} is a numeric scalar resulting from a call to
#'     \code{\link[mvrsquared]{calc_rsquared}}. This slot only populated if
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
new_tidylda <- function(
  lda, 
  dtm, 
  burnin, 
  is_prediction = FALSE,
  alpha = NULL, 
  eta = NULL,
  optimize_alpha = NULL, 
  calc_r2 = NULL,
  calc_likelihood = NULL, 
  call = NULL,
  threads
) {
  

  ### format theta ###
  if (burnin > -1) {
    theta <- t(lda$Cd_mean + lda$alpha) # t(t(lda$Cd_mean) + lda$alpha)
    Cd <- lda$Cd_mean
  } else {
    theta <- t(lda$Cd + lda$alpha) # t(t(lda$Cd) + lda$alpha)
    Cd <- lda$Cd
  }
  
  theta <- t(theta)
  
  theta <- theta / rowSums(theta)
  
  theta[is.na(theta)] <- 0 # just in case of a numeric issue
  
  colnames(theta) <- seq_len(ncol(theta))
  
  rownames(theta) <- rownames(dtm)
  
  ### format beta and all the rest ###
  
  if (!is_prediction) {
    ### format posteriors correctly ###
    if (burnin > -1) { # if you used burnin iterations use Cd_mean etc.
      beta <- lda$Cv_mean + lda$eta
      Cv <- lda$Cv_mean
    } else { # if you didn't use burnin use standard counts (Cd etc.)
      beta <- lda$Cv + lda$eta
      Cv <- lda$Cv
    }
    
    beta <- beta / rowSums(beta)
    
    beta[is.na(beta)] <- 0 # just in case of a numeric issue
    
    colnames(beta) <- colnames(dtm)
    
    rownames(beta) <- colnames(theta)
    
    
    ### collect the results ###
    
    # lambda
    lambda <- calc_lambda(
      beta = beta, theta = theta,
      p_docs = Matrix::rowSums(dtm)
    )
    
    # eta
    colnames(lda$eta) <- colnames(beta)
    
    if (eta$eta_class == "scalar") {
      eta_out <- lda$eta[1, 1]
    } else if (eta$eta_class == "vector") {
      eta_out <- lda$eta[1, ]
    } else if (eta$eta_class == "matrix") {
      eta_out <- lda$eta
    } else { # this should be impossible, but science is hard and I am dumb.
      eta_out <- lda$eta
      
      warning("something went wrong formatting eta. refit.tidylda and predict.tidylda might be affected")
    }
    
    # alpha
    
    if (alpha$alpha_class == "scalar" & !optimize_alpha) {
      alpha_out <- lda$alpha[1]
    } else if (alpha$alpha_class == "vector" | optimize_alpha) {
      alpha_out <- lda$alpha
      
      names(alpha_out) <- rownames(beta)
    } else { # this should be impossible, but science is hard and I am dumb.
      alpha_out <- lda$alpha
      
      warning("something went wrong formatting alpha. refit.tidylda and predict.tidylda might be affected")
    }
    
    # resulting object
    summary <- try(
      tryCatch(
      summarize_topics(beta = beta, theta = theta, dtm = dtm),
      error = function(err){
        err$message <- "summarize_topics failed. model$summary corrupted."
        stop(err)
      })
    )

    log_likelihood <- as_tibble(data.frame(
      iteration = lda$log_likelihood[1, ],
      log_likelihood = lda$log_likelihood[2, ]
    ))
    
    result <- list(
      beta = beta,
      theta = theta,
      lambda = lambda,
      alpha = alpha_out,
      eta = eta_out,
      summary = summary,
      call = call,
      log_likelihood = log_likelihood,
      counts = list(
        Cd = Cd,
        Cv = Cv
      )
    )
    
    class(result) <- "tidylda"
    
    ### calculate and add other things ###
    
    # goodness of fit
    if (calc_r2) {
      result$r2 <- try(
        tryCatch(
          calc_lda_r2(
            dtm = dtm,
            theta = theta,
            beta = beta,
            threads
          ),
          error = function(err){
            err$message <- "calc_r2 failed. R-squared corrupted."
            stop(err)
          } 
        )
      )
    }
    
    # a little cleanup here
    if (!calc_likelihood) {
      result$log_likelihood <- NULL
    }
  }
  
  
  ### return the final result ###
  if (is_prediction) {
    return(theta)
  } else {
    return(result)
  }
}

#' Calculate R-squared for a tidylda Model
#' @keywords internal
#' @description Formats inputs and hands off to mvrsquared::calc_rsquared
#' @param dtm must be of class dgCMatrix
#' @param theta a theta matrix
#' @param beta a beta matrix
#' @param threads number of parallel threads
calc_lda_r2 <- function(dtm, theta, beta, threads) {
  
  # weight rows of theta by document length
  x <- Matrix::rowSums(dtm) * theta
  
  # claculate r-squared
  r2 <- mvrsquared::calc_rsquared(
    y = dtm,
    yhat = list(x = x, w = beta),
    return_ss_only = FALSE,
    threads = threads
  )
  
  # had an issue with preserving names
  names(r2) <- NULL
  
  r2
}

#' Utility function to tidy a simple triplet matrix
#' @keywords internal
#' @param x Object with rownames and colnames
#' @param triplets A data frame or list of i, j, x
#' @param row_names rownames, if not gotten from rownames(x)
#' @param col_names colnames, if not gotten from colnames(x)
#' @note This function ported from \code{\link[tidytext]{tidytext}}, copyright
#'   2017 David Robinson and Julia Silge. Moved the function here for stability
#'   reasons, as it is internal to tidytext
tidy_triplet <- function(x, triplets, row_names = NULL, col_names = NULL) {
  row <- triplets$i
  if (!is.null(row_names)) {
    row <- row_names[row]
  } else if (!is.null(rownames(x))) {
    row <- rownames(x)[row]
  }
  col <- triplets$j
  if (!is.null(col_names)) {
    col <- col_names[col]
  } else if (!is.null(colnames(x))) {
    col <- colnames(x)[col]
  }
  
  ret <- tibble::tibble(row = row, column = col, value = triplets$x)
  ret
}

#' Create a tidy tibble for a dgCMatrix
#' @keywords internal
#' @description Create a tidy tibble for a dgCMatrix. Will probably be a PR to
#'   \link[tidytext]{tidytext} in the future
#' @param x must be of class dgCMatrix
#' @param ... Extra arguments, not used
#' @return Returns a triplet matrix with columns "document", "term", and "count"
tidy_dgcmatrix <- function(x, ...) {
  triplets <- Matrix::summary(methods::as(x, "dgTMatrix"))
  ret <- tidy_triplet(x, triplets)
  colnames(ret) <- c("document", "term", "count")
  ret
}

#' Calculate a matrix whose rows represent P(topic_i|tokens)
#' @keywords internal
#' @description 
#' Use Bayes' rule to get P(topic|token) from the estimated parameters of a
#' probabilistic topic model.This resulting "lambda" matrix can be used for
#' classifying new documents in a frequentist context and supports
#' \code{\link[tidylda]{augment}}.
#' @param theta a theta matrix
#' @param beta a beta matrix
#' @param p_docs A numeric vector of length \code{nrow(theta)} that is
#'   proportional to the number of terms in each document,  defaults to NULL.
#' @param correct Logical. Do you want to set NAs or NaNs in the final result to
#'   zero? Useful when hitting computational underflow. Defaults to \code{TRUE}.
#'   Set to \code{FALSE} for troubleshooting or diagnostics.
#' @return
#' Returns a \code{matrix} whose rows correspond to topics and whose columns
#' correspond to tokens. The i,j entry corresponds to P(topic_i|token_j)
calc_lambda <- function(beta, theta, p_docs = NULL, correct = TRUE){
  
  # set up constants
  D <- nrow(theta)
  K <- ncol(theta)
  V <- ncol(beta)
  
  # probability of each document (assumed to be equiprobable)
  if(is.null(p_docs)){
    p_d <- rep(1/nrow(theta), nrow(theta))
  }else{
    if(sum(is.na(p_docs)) > 0){
      warning("found missing values in p_docs. Setting them as 0.")
      p_docs[ is.na(p_docs) ] <- 0 
    }
    p_d <- p_docs / sum(p_docs)
  }
  
  # get the probability of each topic
  p_t <- p_d %*% theta
  
  # get the probability of each word from the model    
  p_w <- p_t %*% beta
  
  
  
  # get our result
  lambda <- matrix(0, ncol=ncol(p_t), nrow=ncol(p_t))
  diag(lambda) <- p_t
  
  lambda <- lambda %*% beta
  
  lambda <- t(apply(lambda, 1, function(x) x / p_w))
  
  rownames(lambda) <- rownames(beta)
  colnames(lambda) <- colnames(beta)
  
  # give us zeros instead of NAs when we have NA or NaN entries
  if (correct) {
    lambda[is.na(lambda)] <- 0 
  }
  
  return(lambda)
}

#' Probabilistic coherence of topics
#' @description Calculates the probabilistic coherence of a topic or topics. 
#' This approximates semantic coherence or human understandability of a topic.
#' @param beta A numeric matrix or a numeric vector. The vector, or rows of the 
#' matrix represent the numeric relationship between topic(s) and terms. For
#' example, this relationship may be p(word|topic) or p(topic|word).
#' @param data A document term matrix or term co-occurrence matrix. The preferred
#'   class is a \code{\link[Matrix]{dgCMatrix-class}}. However there is support
#'   for any \code{\link[Matrix]{Matrix-class}} object as well as several other
#'   commonly-used classes such as \code{\link[base]{matrix}},
#'   \code{\link[quanteda]{dfm}}, \code{\link[tm]{DocumentTermMatrix}}, and
#'   \code{\link[slam]{simple_triplet_matrix}}
#' @param m An integer for the number of words to be used in the calculation. 
#' Defaults to 5
#' @return Returns an object of class \code{numeric} corresponding to the 
#' probabilistic coherence of the input topic(s).
#' @details 
#'   For each pair of words {a, b} in the top M words in a topic, probabilistic
#'   coherence calculates P(b|a) - P(b), where {a} is more probable than {b} in
#'   the topic. For example, suppose the top 4 words in a topic are {a, b, c, d}.
#'   Then, we calculate 1. P(a|b) - P(b), P(a|c) - P(c), P(a|d) - P(d)
#'   2. P(b|c) - P(c), P(b|d) - P(d)
#'   3. P(c|d) - P(d)
#'   All 6 differences are averaged together.
#' @examples
#' # Load a pre-formatted dtm and topic model
#' data(nih_sample_topic_model)
#' data(nih_sample_dtm) 
#' 
#' calc_prob_coherence(beta = nih_sample_topic_model$phi, data = nih_sample_dtm, m = 5)
#' @export 
calc_prob_coherence <- function(beta, data, m = 5){
  # code below is ported almost verbatim from textmineR. Copied here to reduce
  # cross dependencies between textmineR and tidylda
  
  # beta is a numeric matrix or numeric vector?
  if( ! is.numeric(beta) ){
    stop("beta must be a numeric matrix whose rows index topics and columns\n",
         " index terms or beta must be a numeric vector whose entries index terms.")
  }
  # Ensure dtm is of class dgCMatrix
  dtm <- convert_dtm(dtm = data)
  
  # is m numeric? If it is not an integer, give a warning.
  if( ! is.numeric(m) | m < 1){
    stop("M must be an integer in 1:ncol(beta) or 1:length(beta)")
  }
  
  if(length(m) != 1){
    warning("m is a vector when scalar is expected. Taking only the first value")
    m <- m[ 1 ]
  }
  
  if(floor(m) != m){
    warning("m is expected to be an integer. floor(m) is being used.")
    m <- floor(m)
  }
  
  # # dtm has colnames?
  # if( is.null(colnames(dtm))){
  #   stop("dtm must have colnames")
  # }
  
  # Names of beta in colnames(dtm)
  if( ! is.matrix(beta) ){
    if(sum(names(beta)[ 1:m ] %in% colnames(dtm)) != length(1:m)){
      stop("vocabulary of beta (i.e., colnames(beta)) does not match vocabulary of data")
    }
  }else if(sum(colnames(beta)[ 1:m ] %in% colnames(dtm)) != length(1:m)){
    stop("vocabulary of beta (i.e., colnames(beta)) does not match vocabulary of data")
  }
  
  # Declare a function to get probabilistic coherence on one topic
  pcoh <- function(topic, dtm, m){
    terms <- names(topic)[order(topic, decreasing = TRUE)][1:m]
    dtm.t <- dtm[, terms]
    dtm.t[dtm.t > 0] <- 1
    count.mat <- Matrix::t(dtm.t) %*% dtm.t
    num.docs <- nrow(dtm)
    p.mat <- count.mat/num.docs
    # result <- sapply(1:(ncol(count.mat) - 1), function(x) {
    #   mean(p.mat[x, (x + 1):ncol(p.mat)]/p.mat[x, x] - Matrix::diag(p.mat)[(x + 
    #                                                                           1):ncol(p.mat)], na.rm = TRUE)
    # })
    # mean(result, na.rm = TRUE)
    result <- sapply(1:(ncol(count.mat) - 1), function(x) {
      p.mat[x, (x + 1):ncol(p.mat)]/p.mat[x, x] - 
        Matrix::diag(p.mat)[(x + 1):ncol(p.mat)]
    })
    mean(unlist(result), na.rm = TRUE) 
  }
  
  # if beta is a single topic vector get that one coherence
  if( ! is.matrix(beta) ){
    return(pcoh(topic = beta, dtm = dtm, m = m))
  }
  
  # Otherwise, do it for all the topics
  apply(beta, 1, function(x){
    pcoh(topic = x, dtm = dtm, m = m)
  })
}

