################################################################################
# Experimenting with a different (better?) way to get counts for model updates
################################################################################


library(tidytext)
library(tidyverse)
library(tidylda)
library(Matrix)
### Initial set up ---
# load some documents
docs <- textmineR::nih_sample
# tokenize using tidytext's unnest_tokens
tidy_docs <- docs %>%
  select(APPLICATION_ID, ABSTRACT_TEXT) %>%
  unnest_tokens(output = word,
                input = ABSTRACT_TEXT,
                stopwords = stop_words$word,
                token = "ngrams",
                n_min = 1, n = 2) %>%
  count(APPLICATION_ID, word) %>%
  filter(n>1) #Filtering for words/bigrams per document, rather than per corpus
tidy_docs <- tidy_docs %>% # filter words that are just numbers
  filter(! stringr::str_detect(tidy_docs$word, "^[0-9]+$"))
# turn a tidy tbl into a sparse dgCMatrix
# note tidylda has support for several document term matrix formats
d <- tidy_docs %>%
  cast_sparse(APPLICATION_ID, word, n)
# let's split the documents into two groups to demonstrate predictions and updates
d1 <- d[1:50, ]
d2 <- d[51:nrow(d), ]
# make sure we have different vocabulary for each data set to simulate the "real world"
# where you get new tokens coming in over time
d1 <- d1[, colSums(d1) > 0]
d2 <- d2[, colSums(d2) > 0]
set.seed(123)
lda <- tidylda(
  dtm = d1,
  k = 10,
  iterations = 200,
  burnin = 175,
  alpha = 0.1, # also accepts vector inputs
  beta = 0.05, # also accepts vector or matrix inputs
  optimize_alpha = FALSE, # experimental
  calc_likelihood = TRUE,
  calc_r2 = TRUE, # see https://arxiv.org/abs/1911.11061
  return_data = FALSE
)


#### Begin getting counts for new transfer learning nonsense ----

# function to recover counts for both theta and phi
recover_counts <- function(prob_matrix, prior_matrix, total_vector) {
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
      
      idx <- sample(seq_along(x), remainder, prob = sample_prob)
      
      round_x[idx] <- round_x[idx] + 1
      
      return(round_x)
      
    } else { # we need to take some away
      
      idx <- sample(seq_along(x)[round_x > 0], -1 * remainder, prob = x[round_x > 0])
      
      round_x[idx] <- round_x[idx] - 1
      
      return(round_x)
    }
    
  })
  
  count_matrix <- t(count_matrix)
  
  
  count_matrix
}

# begin with dot product to get topic distributions for each document
theta_hat <- predict(lda, d1, method = "dot", no_common_tokens = "uniform")

alph <- tidylda:::format_alpha(lda$alpha, nrow(lda$phi))

alph <- matrix(1, ncol = nrow(theta_hat), nrow = ncol(theta_hat)) + alph$alpha

alph <- t(alph)

# get Cd
Cd <- recover_counts(
  prob_matrix = theta_hat,
  prior_matrix = alph,
  total_vector = Matrix::rowSums(d1)
  )

# Cd is done now get Ck 
Ck <- colSums(Cd)

# use Ck to get counts for Cv
# in future, you'll have to do this *after* you reconcile vocab and add uniform
# counts over new words
bet <- tidylda:::format_beta(beta = lda$beta, k = ncol(Cd), Nv = ncol(d1))

Cv <- recover_counts(
  prob_matrix = lda$phi,
  prior_matrix = bet$beta,
  total_vector = Ck
)


