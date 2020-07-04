################################################################################
# This script is a simple LDA gibbs sampler in R
# I'm using it to help troubleshoot and optimize my C++ code
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

iterations <- 200

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

for (t in 1:iterations) { # for each iteration
  
  for (d in 1:Nd) { # for each document
    
    for (n in 1:length(docs[[d]])) { # for each token in that document
      
      # decrement counts from previous iteration where we hit this token
      # purpose to correctly calculate the probabilities
      Cd[Zd[[d]][n], d] <- Cd[Zd[[d]][n], d] - 1
      
      Cv[Zd[[d]][n], docs[[d]][n]] <- Cv[Zd[[d]][n], docs[[d]][n]] - 1
      
      Ck[Zd[[d]][n]] <- Ck[Zd[[d]][n]] - 1
      
      # calculate the probability of sampling a topic at this location
      qz <- numeric(Nk) # initialize topic probability vector
      
      for (k in 1:Nk) {
        qz[k] <- 
          (Cv[k, docs[[d]][n]] + beta[k, docs[[d]][n]]) / (Ck[k] + sum_beta) *
          (Cd[k, d] + alpha[k]) / (length(docs[[d]]) + sum_alpha - 1)
      }
      
      # sample a topic at random
      Zd[[d]][n] <- sample(1:Nk, 1, prob = qz)
      
      # increase the counts of Ck, Cd, Cv where we just sampled
      Cd[Zd[[d]][n], d] <- Cd[Zd[[d]][n], d] + 1
      
      Cv[Zd[[d]][n], docs[[d]][n]] <- Cv[Zd[[d]][n], docs[[d]][n]] + 1
      
      Ck[Zd[[d]][n]] <- Ck[Zd[[d]][n]] + 1
      
    } # end loop over token indices
    
  } # end loop over docs
  
} # end iterations


# At this point, I more or less expel everything from C++ and do residual
# calculations to properly construct the posterior parameter estimates
