context("tests of posterior methods")

library(dplyr)

dtm <- nih_sample_dtm

d1 <- dtm[1:50, ]

# make sure we have different vocabulary for each data set
d1 <- d1[, Matrix::colSums(d1) > 0]

lda <- tidylda(
  data = d1,
  k = 4,
  iterations = 20, burnin = 10,
  alpha = 0.1, eta = 0.05,
  optimize_alpha = TRUE,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  return_data = FALSE,
  verbose = FALSE
)


test_that("posterior methods function with good inputs",{
  
  # sample from theta 1 doc
  g <- posterior(
    x = lda,
    matrix = "theta",
    which = 1,
    times = 10
  )
  
  expect_equal(
    g %>% 
      group_by(document, sample) %>% 
      summarise(tot = sum(theta)) %>% 
      .[["tot"]] %>% 
      sum(),
    10
  )
  
  # sample from theta many docs
  g <- posterior(
    x = lda,
    matrix = "theta",
    which = 1:3,
    times = 10
  )
  
  expect_equal(
    g %>% 
      group_by(document, sample) %>% 
      summarise(tot = sum(theta)) %>% 
      .[["tot"]] %>% 
      sum(),
    30
  )
  
  # sample from beta 1 doc
  g <- posterior(
    x = lda,
    matrix = "beta",
    which = 1,
    times = 10
  )
  
  expect_equal(
    g %>% 
      group_by(token, sample) %>% 
      summarise(tot = sum(beta)) %>% 
      .[["tot"]] %>% 
      sum(),
    10
  )
  
  # sample from beta many docs
  g <- posterior(
    x = lda,
    matrix = "beta",
    which = 1:3,
    times = 10
  )
  
  expect_equal(
    g %>% 
      group_by(token, sample) %>% 
      summarise(tot = sum(beta)) %>% 
      .[["tot"]] %>% 
      sum(),
    30
  )
  
  
})


test_that("posterior methods throw errors when they should",{
  
  # malformed matrix
  expect_error(
    g <- posterior(
      x = lda,
      matrix = "something",
      which = 1,
      times = 10
    ),
    label = "malformed matrix"
  )
  
  # which has NAs
  expect_error(
    g <- posterior(
      x = lda,
      matrix = "theta",
      which = c(1, NA),
      times = 10
    ),
    label = "which has NA"
  )
  
  # which has negatives
  expect_error(
    g <- posterior(
      x = lda,
      matrix = "theta",
      which = c(1, -2),
      times = 10
    ),
    label = "which has negative"
  )
  
  # times is too long
  expect_error(
    g <- posterior(
      x = lda,
      matrix = "theta",
      which = c(1),
      times = 10:11
    ),
    label = "times is too long"
  )
  
  # times is not numeric
  expect_error(
    g <- posterior(
      x = lda,
      matrix = "theta",
      which = c(1),
      times = "10:11"
    ),
    label = "times is not numeric"
  )
  
  # times is negative
  expect_error(
    g <- posterior(
      x = lda,
      matrix = "theta",
      which = c(1),
      times = -10
    ),
    label = "times is negative"
  )
  
})