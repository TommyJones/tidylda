context("tests of posterior methods")

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
  
  p <- posterior(lda)
  
  # appropriate class
  expect_s3_class(p, "tidylda_posterior")
  
  # dimensions of matrices
  expect_equal(ncol(p$theta_par), nrow(lda$theta))
  expect_equal(nrow(p$theta_par), ncol(lda$theta))
  expect_equal(ncol(p$beta_par), nrow(lda$beta))
  expect_equal(nrow(p$beta_par), ncol(lda$beta))
  
  # sample from theta 1 doc
  g <- generate(
    x = p,
    matrix = "theta",
    which = 1,
    times = 10
  )
  
  # sample from theta many docs
  g <- generate(
    x = p,
    matrix = "theta",
    which = 1:3,
    times = 10
  )
  
  # sample from beta 1 doc
  g <- generate(
    x = p,
    matrix = "beta",
    which = 1,
    times = 10
  )
  
  # sample from beta many docs
  g <- generate(
    x = p,
    matrix = "beta",
    which = 1:3,
    times = 10
  )
  
})


test_that("posterior methods throw errors when they should",{
  
  p <- posterior(lda)
  
  # malformed matrix
  expect_error(
    g <- generate(
      x = p,
      matrix = "something",
      which = 1,
      times = 10
    ),
    label = "malformed matrix"
  )
  
  # which has NAs
  expect_error(
    g <- generate(
      x = p,
      matrix = "theta",
      which = c(1, NA),
      times = 10
    ),
    label = "which has NA"
  )
  
  # which has negatives
  expect_error(
    g <- generate(
      x = p,
      matrix = "theta",
      which = c(1, -2),
      times = 10
    ),
    label = "which has negative"
  )
  
  # times is too long
  expect_error(
    g <- generate(
      x = p,
      matrix = "theta",
      which = c(1),
      times = 10:11
    ),
    label = "times is too long"
  )
  
  # times is not numeric
  expect_error(
    g <- generate(
      x = p,
      matrix = "theta",
      which = c(1),
      times = "10:11"
    ),
    label = "times is not numeric"
  )
  
  # times is negative
  expect_error(
    g <- generate(
      x = p,
      matrix = "theta",
      which = c(1),
      times = -10
    ),
    label = "times is negative"
  )
  
})