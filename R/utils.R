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

#' Format \code{beta} For Input into \code{fit_lda_c}
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{beta} but the C++ Gibbs
#'   sampler in \code{\link[tidylda]{fit_lda_c}} only takes it one way. This function does the
#'   appropriate formatting. It also returns errors if the user input a malformatted
#'   \code{beta}.
#' @param beta the prior for words over topics. Can be a numeric scalar, numeric
#'   vector, or numeric matrix.
#' @param k the number of topics.
#' @param Nv the total size of the vocabulary as inherited from \code{ncol(dtm)}
#'   in \code{\link[tidylda]{tidylda}}.
#' @return
#'   Returns a list with two elements: \code{beta} and \code{beta_class}.
#'   \code{beta} is the post-formatted version of \code{beta} in the form of a
#'   \code{k} by \code{Nv} numeric matrix. \code{beta_class} is a character
#'   denoting whether or not the user-supplied \code{beta} was a "scalar",
#'   "vector", or "matrix".
format_beta <- function(beta, k, Nv) {
  if (!is.numeric(beta) | sum(is.na(beta)) > 0 | sum(beta == 0) == length(beta)) {
    stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
  }

  if (length(beta) == 1) { # if beta is a scalar

    beta <- matrix(beta, nrow = k, ncol = Nv)

    beta_class <- "scalar"
  } else if (is.vector(beta)) { # if beta is a vector

    if (length(beta) != Nv) { # if you didn't specify this vector right
      stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    }

    # otherwise let's carry on...
    # make beta a matrix to format for C++ funciton
    beta <- t(beta + matrix(0, nrow = length(beta), ncol = k))

    beta_class <- "vector"
  } else if (is.matrix(beta)) { # if beta is a matrix

    # check dims before moving on
    if (nrow(beta) != k | ncol(beta) != Nv) {
      stop(
        "If beta is a matrix, it must have the same number of rows as topics ",
        "and it must have the same number of columns (tokens) as your dtm. ",
        "But I see nrow(beta) = ", nrow(beta), " and k = ", k, ". I also see ",
        "ncol(beta) = ", ncol(beta), " but ncol(dtm) = ", Nv
      )
    }

    beta_class <- "matrix"
  }


  list(
    beta = beta,
    beta_class = beta_class
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

#' Get Count Matrices from Phi or Theta (and Priors)
#' @keywords internal
#' @description
#'   This function is a core component of \code{\link[tidylda]{initialize_topic_counts}}.
#'   See details, below.
#' @param prob_matrix a numeric \code{phi} or \code{theta} matrix
#' @param prior_matrix a matrix of same dimension as \code{prob_matrix} whose 
#'   entries represent the relevant prior (\code{alpha} or \code{beta})
#' @param total_vector a vector of token counts of length \code{ncol(prob_matrix)}
#' @details 
#'   This function uses a probability matrix (theta or phi), its prior (alpha or
#'   beta, respectively), and a vector of counts to simulate what the the Cd or
#'   Cv matrix would be at the end of a Gibbs run that resulted in that probability
#'   matrix.
#'   
#'   For example, theta is calculated from a matrix of counts, Cd, and a prior,
#'   alpha. Specifically, the i,j entry of theta is given by
#'   
#'   \code{(Cd[i, j] + alpha[i, j]) / sum(Cd[, j] + alpha[, j])}
#'   
#'   Similarly, phi comes from
#'   
#'   \code{(Cv[i, j] + beta[i, j]) / sum(Cv[, j] + beta[, j])}
#'   
#'   (The above are written to be general with respect to alpha and beta being
#'   matrices. They could also be vectors or scalars.)
#'   
#'   So, this function uses the above formulas to try and reconstruct Cd or Cv
#'   from theta and alpha or phi and beta, respectively. As of this writing,
#'   this method is experimental. In the future, there will be a paper with
#'   more technical details cited here.
#'   
#'   The priors must be matrices for the purposes of the function. This is to
#'   support topic seeding and model updates. The former requires beta to be a 
#'   matrix. The latter may require beta to be a matrix. Here, alpha is also
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
            sample(seq_along(x), remainder, prob = sample_prob)
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
            sample(x = sample_from, size = sample_size, prob = sample_prob)
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
#'   priors \code{alpha} and \code{beta}.
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @param k the number of topics
#' @param alpha the numeric vector prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_alpha}}
#' @param beta the numeric matrix prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_beta}}
#' @param phi_initial if specified, a numeric matrix for the probability of tokens
#'   in topics. Must be specified for predictions or updates as called by
#'   \code{\link[tidylda]{predict.tidylda}} or \code{\link[tidylda]{refit.tidylda}}
#'   respectively.
#' @param theta_initial if specified, a numeric matrix for the probability of
#'   topics in documents. Must be specified for updates as called by
#'   \code{\link[tidylda]{refit.tidylda}}
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

    phi_initial <- apply(beta, 1, function(x) {
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

  # initialize Cd by calling recover_counts_from_probs
  # we don't need to initialize Cv because we can use the probabilities in phi,
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
  # if we have more than 3000 documents, do it in parallel with furrr::future_map

  batches <- seq(1, nrow(dtm), by = 3000)

  lexicon <- furrr::future_map(
    .x = batches,
    .f = function(b) {
      rows <- b:min(b + 2999, nrow(dtm))

      # if statement to handle single observations
      if (length(rows) == 1) {
        cd_tmp <- Cd_start

        dtm_tmp <- dtm
      } else {
        cd_tmp <- Cd_start[rows, ]

        dtm_tmp <- dtm[rows, ]
      }

      l <- create_lexicon(
        Cd = cd_tmp,
        Phi = phi_initial,
        dtm = dtm_tmp,
        alpha = alpha,
        freeze_topics = freeze_topics
      )
    }, ...
  )

  # combine
  Zd <- Reduce("c", lapply(lexicon, function(l) l$Zd))

  docs <- Reduce("c", lapply(lexicon, function(l) l$docs))

  Cv <- Reduce("+", lapply(lexicon, function(l) l$Cv))

  Ck <- Reduce("+", lapply(lexicon, function(l) l$Ck))

  Cd <- do.call(rbind, lapply(lexicon, function(l) l$Cd))

  out <- list(
    docs = docs,
    Zd = Zd,
    Cd = Cd,
    Cv = Cv,
    Ck = Ck
  )

  out
}

#' Summarize a topic model consistently across methods/functions
#' @keywords internal
#' @description
#'   Summarizes topics in a model. Called by \code{\link[tidylda]{tidylda}}
#'   and \code{\link[tidylda]{refit.tidylda}} and used to augment
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
summarize_topics <- function(theta, phi, dtm) {

  # probabilistic coherence with default M = 5
  coherence <- textmineR::CalcProbCoherence(phi = phi, dtm = dtm)

  # prevalence of each topic, weighted by terms
  prevalence <- Matrix::rowSums(dtm) * theta

  prevalence <- colSums(prevalence) / sum(prevalence)

  prevalence <- round(prevalence * 100, 2)

  # top 3 terms
  top_terms <- t(textmineR::GetTopTerms(phi, 3))

  top_terms <- apply(top_terms, 1, function(x) {
    paste(c(x, "..."), collapse = ", ")
  })

  # combine into a summary
  summary <- data.frame(
    topic = as.numeric(rownames(phi)),
    prevalence = prevalence,
    coherence = coherence,
    top_terms = top_terms,
    stringsAsFactors = FALSE
  )

  summary <- tibble::as_tibble(summary)

  summary
}

#' Format the outputs of \code{\link[tidylda]{fit_lda_c}} consistently
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
#' @param beta output of \code{\link[tidylda]{format_beta}}
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
#'   \code{r2} is a numeric scalar resulting from a call to
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
format_raw_lda_outputs <- function(lda, dtm, burnin, is_prediction = FALSE,
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

  if (!is_prediction) {
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
    gamma <- textmineR::CalcGamma(
      phi = phi, theta = theta,
      p_docs = Matrix::rowSums(dtm)
    )

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
    summary <- summarize_topics(phi = phi, theta = theta, dtm = dtm)

    log_likelihood <- as_tibble(data.frame(
      iteration = lda$log_likelihood[1, ],
      log_likelihood = lda$log_likelihood[2, ]
    ))

    result <- new_tidylda(
      phi = phi,
      theta = theta,
      gamma = gamma,
      alpha = alpha_out,
      beta = beta_out,
      summary = summary,
      call = call,
      log_likelihood = log_likelihood
    )

    ### calculate and add other things ---

    # goodness of fit
    if (calc_r2) {
      result$r2 <- calc_lda_r2(
        dtm = dtm,
        theta = theta,
        phi = phi,
        batch_size = 3000, # hard coded for now
        ...
      )
    }

    # a little cleanup here
    if (!calc_likelihood) {
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

#' Calculate R-squared for a tidylda Model
#' @keywords internal
#' @description Formats inputs and hands off to mvrsquared::calc_rsquared
#' @param dtm must be of class dgCMatrix
#' @param theta a theta matrix
#' @param phi a phi matrix
#' @param batch_size for parallel processing
#' @param ... other arguments passed to \code{\link[furrr]{future_map}}
calc_lda_r2 <- function(dtm, theta, phi, batch_size, ...) {
  # divide things into batches
  batches <- furrr::future_map(
    .x = seq(1, nrow(dtm), by = batch_size),
    .f = function(b){
      # rows to select on
      rows <- b:min(b + batch_size - 1, nrow(dtm))
      
      # rows of the dtm
      y <- dtm[rows, ]
      
      if (! inherits(y, "dgCMatrix")) {
        y <- methods::as(object = y, Class = "dgCMatrix")
      }
      
      # rows of theta multiplied by doc lengths in y
      x <- Matrix::rowSums(y) * theta[rows, ]
      
      if (! inherits(x, "matrix")) {
        x <- matrix(x, nrow = 1)
      }
      
      list(y = y, x = x)
    }, 
    ...
  )
  
  ybar <- Matrix::colMeans(dtm)
  
  # MAP: calculate sum of squares
  ss <- furrr::future_map(
    .x = batches,
    .f = function(batch) {
      mvrsquared::calc_rsquared(
        y = batch$y,
        yhat = list(x = batch$x, w = phi),
        ybar = ybar,
        return_ss_only = TRUE)
    },
    ...
  )
  
  # REDUCE: get SST and SSE by summation
  ss <- colSums(
    do.call(rbind, ss)
  )
  
  # return r-squared
  r2 <- 1 - ss["sse"] / ss["sst"]
  
  names(r2) <- NULL
  
  r2
}


