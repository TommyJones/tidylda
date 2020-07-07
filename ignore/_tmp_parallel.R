################################################################################
# This script experiments with parallel Gibbs algorithms in R.
# Syntax is C++-ish rather than optimized for R. Don't @ me.
# If successful here, I can port to C++
################################################################################

library(tidylda)

set.seed(90210)

# the below calls some internal tidylda functions
# on style, I'm going to use c++ style variables and be explicit about allocating
# the variables in a way you don't have to do in R.
# For example, if I say alpha = 0.1, my R functions know to make alpha = seq(0.1, k)

### Initial set up of variables ----
dtm <- textmineR::nih_sample_dtm

Nk <- 10 # number of topics

Nd <- nrow(dtm)

Nv <- ncol(dtm)

alpha <- seq(0.1, Nk)

# a matrix b/c of the transfer learning stuff I'm preparing to do
# FWIW textmineR's sampler takes beta as a matrix too
beta <- matrix(0.05, ncol = ncol(dtm), nrow = Nk) 

iterations <- 100

# not going to optimize alpha, calculate the likelihood, or do burnin iterations
# this will keep things simpler

### Allocate some counts to format data as it goes into the C++ gibbs sampler ----

counts <- 
  tidylda:::initialize_topic_counts(
    dtm = dtm, 
    k = Nk,
    alpha = alpha, 
    beta = beta,
    threads = 1
  )

# pull out the individual objects as in the C++ gibbs sampler

# docs is a list of integer vectors
# length(docs[[d]]) == sum(dtm[d, ])
# the integers correspond to the column indices of dtm, Cv, and beta
docs <- counts$docs

docs <- lapply(docs, function(x) x + 1) # r is 1 indexed

# Zd is very similar to docs, except that its indices correspond to topics
# its entries index the rows of Cv, beta, and Cd (after Cd transpose, below)
# length(Zd[[d]]) == sum(dtm[d, ])
Zd <- counts$Zd

Zd <- lapply(Zd, function(x) x + 1) # r is 1 indexed

# Ck is a count of the number of times each topic has been sampled for each
# instance of a token in the whole corpus
Ck <- as.numeric(counts$Ck) # as.numeric b/c Ck comes out as a single column matrix

# Cd is the count of topics sampled for each document
# columns are documents, rows are topics
Cd <- counts$Cd

# > dim(Cd)
# [1] 100  10

# Cv is the count of topics sampled for each token
Cv <- counts$Cv

# > dim(Cv)
# [1]   10 5210

### Start sampling! ----

# initialize variables used only in the sampler

sum_beta <- sum(beta[1,]) # strict LDA has beta as a vector not a matrix, so sum only a row

sum_alpha <- sum(alpha)

# initialize objects used for parallel sampling
threads <- 4

Cv_list <- vector(mode = "list", length = threads)

Ck_list <- vector(mode = "list", length = threads)

batch_size <- ncol(Cd) / threads # may fail if not evenly divided. Can fancy later

batch_indices <- seq(1, ncol(Cd), by = round(batch_size))

batch_indices <- lapply(batch_indices, function(x) x:min(x + batch_size - 1, ncol(Cd)))


for (t in 1:iterations) { # for each iteration
  
  # copy global Ck and Cv to their lists
  Cv_list <- lapply(Cv_list, function(x) Cv)
  
  Ck_list <- lapply(Ck_list, function(x) Ck)
  
  # loop over batches, documents, tokens
  for (j in 1:threads) { # this will be parallel 
    for (d in batch_indices[[j]]) {
      
      for (n in 1:length(docs[[d]])) { # for each token in that document
        
        # decrement counts from previous iteration where we hit this token
        # purpose to correctly calculate the probabilities
        Cd[Zd[[d]][n], d] <- Cd[Zd[[d]][n], d] - 1
        
        Cv_list[[j]][Zd[[d]][n], docs[[d]][n]] <- Cv_list[[j]][Zd[[d]][n], docs[[d]][n]] - 1
        
        Ck_list[[j]][Zd[[d]][n]] <- Ck_list[[j]][Zd[[d]][n]] - 1
        
        # calculate the probability of sampling a topic at this location
        qz <- numeric(Nk) # initialize topic probability vector
        
        for (k in 1:Nk) {
          qz[k] <- 
            (Cv_list[[j]][k, docs[[d]][n]] + beta[k, docs[[d]][n]]) / (Ck_list[[j]][k] + sum_beta) *
            (Cd[k, d] + alpha[k]) / (length(docs[[d]]) + sum_alpha - 1)
        }
        
        # sample a topic at random
        Zd[[d]][n] <- sample(1:Nk, 1, prob = qz)
        
        # increase the counts of Ck, Cd, Cv where we just sampled
        Cd[Zd[[d]][n], d] <- Cd[Zd[[d]][n], d] + 1
        
        Cv_list[[j]][Zd[[d]][n], docs[[d]][n]] <- Cv_list[[j]][Zd[[d]][n], docs[[d]][n]] + 1
        
        Ck_list[[j]][Zd[[d]][n]] <- Ck_list[[j]][Zd[[d]][n]] + 1
        
      } # end loop over token indices
      
    } # end loop over docs
    
  } # end loop over batches
  
  # reconcile Cd and Ck across cores
  Cv <- (Reduce("+", Cv_list) - threads * Cv) + Cv
  
  Ck <- (Reduce("+", Ck_list) - threads * Ck) + Ck
   
} # end iterations


# At this point, I more or less expel everything from C++ and do residual
# calculations to properly construct the posterior parameter estimates
