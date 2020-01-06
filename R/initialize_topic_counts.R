
# make this an internal function somehow...
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
  
  if (nrow(dtm) <= 3000) { # if a small corpus, do sequential
    
    # Cd_start <- sapply(X = Matrix::rowSums(dtm),
    #                    FUN = cd_sampler)
    
    Cd_start <- mapply(FUN = cd_sampler,
                       size = Matrix::rowSums(dtm),
                       prob = as.list(data.frame(t(theta_initial)))) 
    
    Cd_start <- t(Cd_start)
    
  } else { # otherwise do it in parallel
    
    batches <- seq(1, nrow(dtm), by = 3000)
    
    iterator <- textmineR::TmParallelApply(batches, function(b){
      
      rows <- b:min(b + 2999, nrow(dtm))
      
      if (length(rows) > 1) {
        
        size <- Matrix::rowSums(dtm[rows, ])
        
        prob <- as.list(data.frame(t(theta_initial[rows, ])))
        
      } else {
        size <- sum(dtm[rows, ])
        
        prob <- theta_initial[rows, ]
      }
      
      list(size = size, prob = prob)
      
    }, export = c("dtm", "theta_initial"),
    ...)
    
    Cd_start <- textmineR::TmParallelApply(X = iterator,
                                           FUN = function(x){
                                             out <- mapply(FUN = cd_sampler,
                                                           size = x$size,
                                                           prob = x$prob)
                                             
                                             t(out)
                                           },
                                           export = "cd_sampler",
                                           ...)
    
    Cd_start <- do.call(rbind, Cd_start)
    
  }
  
  # Initialize objects with that single Gibbs iteration mentioned above
  if(nrow(dtm) > 3000){  # if we have more than 3,000 docs, do it in parallel
    
    batches <- seq(1, nrow(dtm), by = 3000)

    lexicon <- textmineR::TmParallelApply(batches, function(b){
      
      rows <- b:min(b + 2999, nrow(dtm))
      
      l <- create_lexicon(Cd = Cd_start[rows,],
                          Phi = phi_initial,
                          dtm = dtm[rows,],
                          alpha = alpha,
                          freeze_topics = freeze_topics)
      
    }, 
    export = c("alpha", "Cd_start", "phi_initial", "dtm", "freeze_topics"),
    ...)
    
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

  }else{ # if we have 3,000 or fewer docs do it sequentially
    out <- create_lexicon(Cd = Cd_start,
                          Phi = phi_initial,
                          dtm = dtm,
                          alpha = alpha,
                          freeze_topics = freeze_topics)
  }
  
  out
  
}


