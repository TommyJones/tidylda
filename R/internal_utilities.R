# internal functions:

### Makes sure beta is formatted correctly for input to C functions for LDA ----

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

### Makes sure alpha is formatted correctly for input to C functions for LDA ----
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

### Initialize topic counts for gibbs sampling ----
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

### summarize a topic model consistently across methods/functions ----

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

### format the outputs of fit_lda_c consistently ----

format_raw_lda <- function(lda, dtm, burnin, is_prediction = FALSE, 
                           alpha = NULL, beta = NULL, 
                           optimize_alpha, calc_r2 = NULL, 
                           calc_likelihood = NULL, 
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
            error when you ran 'fit_tidylda'.")
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
                   log_likelihood = data.frame(iteration = lda$log_likelihood[1,],
                                               log_likelihood = lda$log_likelihood[2, ])
    ) # add other things here if necessary
    
    class(result) <- "tidylda_model"
    
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
