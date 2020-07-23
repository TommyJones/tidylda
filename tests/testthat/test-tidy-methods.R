context("tests of tidy methods for tidylda")

dtm <- textmineR::nih_sample_dtm

d1 <- dtm[1:50, ]

d2 <- dtm[51:100, ]

# make sure we have different vocabulary for each data set
d1 <- d1[, Matrix::colSums(d1) > 0]

d2 <- d2[, Matrix::colSums(d2) > 0]

lda <- tidylda(
  data =  d1,
  k = 4,
  iterations = 20, burnin = 10,
  alpha = 0.1, eta = 0.05,
  optimize_alpha = TRUE,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  return_data = FALSE
)

### test tidy methods ----

test_that("tidy.tidylda works as expected", {
  
  # tidy phi
  tidy_phi <- tidy(
    x = lda,
    matrix = "phi"
  )
  
  expect_named(tidy_phi, c("topic", "token", "phi"))
  
  expect_type(tidy_phi[[1]], "double")
  
  expect_type(tidy_phi[[2]], "character")
  
  expect_type(tidy_phi[[3]], "double")
  
  # log to tidy phi
  tidy_phi_log <- tidy(
    x = lda,
    matrix = "phi",
    log = TRUE
  )
  
  expect_named(tidy_phi_log, c("topic", "token", "log_phi"))
  
  expect_type(tidy_phi_log[[1]], "double")
  
  expect_type(tidy_phi_log[[2]], "character")
  
  expect_type(tidy_phi_log[[3]], "double")
  
  expect_equal(sum(colnames(lda$phi) %in% tidy_phi$token), length(colnames(lda$phi)))
  
  # tidy theta
  tidy_theta <- tidy(
    x = lda,
    matrix = "theta"
  )
  
  expect_named(tidy_theta, c("document", "topic", "theta"))
  
  expect_type(tidy_theta[[1]], "character")
  
  expect_type(tidy_theta[[2]], "double")
  
  expect_type(tidy_theta[[3]], "double")
  
  expect_equal(sum(colnames(lda$theta) %in% tidy_theta$topic), length(colnames(lda$theta)))
  
  # tidy lambda
  tidy_lambda <- tidy(
    x = lda,
    matrix = "lambda"
  )
  
  expect_named(tidy_lambda, c("topic", "token", "lambda"))
  
  expect_type(tidy_lambda[[1]], "double")
  
  expect_type(tidy_lambda[[2]], "character")
  
  expect_type(tidy_lambda[[3]], "double")
  
  expect_equal(sum(colnames(lda$lambda) %in% tidy_lambda$token), length(colnames(lda$lambda)))
  
})

test_that("tidy throws errors for malformed inputs", {
  expect_error(
    tidy(
      x = lda,
      matrix = "WRONG"
    )
  )
  
  expect_error(
    tidy(
      x = lda,
      matrix = 1
    )
  )
  
  expect_error(
    tidy(
      x = lda,
      matrix = "phi",
      log = "WRONG"
    )
  )
  
  # matrices the same as above
  expect_error(
    tidy_theta <- tidy(
      x = lda$theta,
      matrix = "WRONG"
    )
  )
  
  expect_error(
    tidy_theta <- tidy(
      x = lda$theta,
      matrix = "theta",
      log = "WRONG"
    )
  )
  
})

### tests for the glance method ----

test_that("glance.tidylda behaves nicely", {
  
  # well-formed call
  g <- glance(lda)
  
  expect_named(g, c(
    "num_topics", "num_documents", "num_tokens",
    "iterations", "burnin"
  ))
  
  expect_equal(g$num_topics, nrow(lda$phi))
  
  expect_equal(g$num_documents, nrow(lda$theta))
  
  expect_equal(g$num_tokens, ncol(lda$phi))
  
  expect_equal(g$iterations, lda$call$iterations)
  
  expect_equal(g$burnin, lda$call$burnin)

})

test_that("glance works with updated models", {
  l2 <- refit(lda, d2, iterations = 20)
  
  g <- glance(l2)
  
  expect_named(g, c(
    "num_topics", "num_documents", "num_tokens",
    "iterations", "burnin"
  ))
  
  expect_equal(g$num_topics, nrow(l2$phi))
  
  expect_equal(g$num_documents, nrow(l2$theta))
  
  expect_equal(g$num_tokens, ncol(l2$phi))
  
  expect_equal(g$iterations, l2$call$iterations)
  
  expect_equal(g$burnin, NA)
})

### tests for the glance method ----

# make a tidy tibble
# note this uses unigrams and bigrams to ensure that there isn't 
# 100% overlap in vocabulary between data and model
tidy_docs <- 
  textmineR::nih_sample[1:10, ] %>% 
  dplyr::select(APPLICATION_ID, ABSTRACT_TEXT) %>% 
  tidytext::unnest_tokens(output = word, 
                          input = ABSTRACT_TEXT,
                          stopwords = tidytext::stop_words$word,
                          token = "ngrams",
                          n_min = 1, n = 2) 

test_that("augment.tidylda behaves nicely", {
  
  td <- tidy_docs
  
  colnames(td)[1:2] <- c("document", "term")
  
  # tidy tibble input
  a <- augment(
    x = lda, 
    data = td,
    type = "class"
  )
  
  # right number of rows?
  expect_equal(nrow(a), nrow(tidy_docs))
  
  # correctly identified topics?
  expect_equal(
    a %>% 
      dplyr::filter(term %in% colnames(lda$phi)) %>%
      dplyr::select(topic) %>%
      dplyr::summarise(na_topic = sum(is.na(topic))) %>%
      as.numeric,
    0,
    label = "augment is.na(topic class)"
  )
  
  # sparse dtm input
  
  a <- augment(
    x = lda,
    data = d1,
    type = "prob"
  )
  
  expect_false(
    a %>% 
      dplyr::filter(term %in% colnames(lda$phi)) %>%
      dplyr::select(-c(document, term)) %>%
      colSums() %>%
      sum %>% 
      is.na,
    label = "augment is.na(topic probs)"
  )
  
  # probably should test other dtm types...
  
})


test_that("augment.tidylda reacts to malformed inputs correctly", {
  
  expect_error(
    augment(
      x = "wat",
      data = d1,
      type = "prob"
    )
  )
  
  expect_error(
    augment(
      x = lda,
      data = "wat",
      type = "prob"
    )
  )
  
  expect_error(
    augment(
      x = lda,
      data = d1,
      type = "wat"
    )
  )
  
})
