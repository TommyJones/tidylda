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
  alpha = 1, # also accepts vector inputs
  beta = 0.05, # also accepts vector or matrix inputs
  optimize_alpha = FALSE, # experimental
  calc_likelihood = TRUE,
  calc_r2 = TRUE, # see https://arxiv.org/abs/1911.11061
  return_data = FALSE
)


#### Begin getting counts for new transfer learning nonsense ----

# begin with dot product to get topic distributions for each document
theta_hat <- predict(lda, d1, method = "dot", no_common_tokens = "uniform")

alph <- tidylda:::format_alpha(lda$alpha, nrow(lda$phi))

alph <- matrix(1, ncol = nrow(theta_hat), nrow = ncol(theta_hat)) + alph$alpha

alph <- t(alph)

# get Cd
Cd <- tidylda:::recover_counts_from_probs(
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

Cv <- tidylda:::recover_counts_from_probs(
  prob_matrix = lda$phi,
  prior_matrix = bet$beta,
  total_vector = Ck
)


