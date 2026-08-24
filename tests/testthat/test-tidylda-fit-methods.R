context("tidylda core tests")

### Define some common objects ----

dtm <- nih_sample_dtm

d1 <- dtm[1:50, ]

### Tests for initial fitting of topic models ----

test_that("can fit lda models without error", {

  # if any of the below throw an error, you've got a problem...

  # scalar priors without optimizing alpha
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  # make sure r2 is numeric since calc_r2 = TRUE
  expect_type(lda$r2, "double")
  
  # make sure r2 doesn't have a names element
  expect_null(names(lda$r2))

  # make sure that likelihood is correct since calc_likelihood = TRUE
  expect_s3_class(lda$log_likelihood, "tbl_df")
  
  # Three columns as of 0.1.0: the plug-in likelihood gained the collapsed
  # joint alongside it.
  expect_named(
    lda$log_likelihood, c("iteration", "log_likelihood", "log_joint")
  )
  
  # The likelihood is evaluated every likelihood_every-th iteration (roadmap
  # D11), so the old assertion -- nrow == tail(iteration, 1) + 1 -- no longer
  # holds; it assumed every iteration is recorded. This is NOT thinning: the
  # chain advances every iteration and every post-burnin iteration still
  # contributes to the count sums. Only the read-only diagnostic runs less often.
  #
  # What must still hold: iterations are 0-indexed, strictly increasing, spaced
  # by the interval, and the final iteration is always recorded so the curve
  # ends where the run does.
  expect_true(all(diff(lda$log_likelihood$iteration) > 0))
  expect_equal(lda$log_likelihood$iteration[1], 0)
  expect_equal(tail(lda$log_likelihood$iteration, 1), 20 - 1)

  # With the default interval of 10 over 20 iterations: 0, 10, and the forced
  # final 19.
  expect_equal(lda$log_likelihood$iteration, c(0, 10, 19))

  # An interval of 1 records every iteration, recovering the old expectation.
  lda_every <- tidylda(
    data = d1, k = 4, iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05, calc_likelihood = TRUE,
    verbose = FALSE, likelihood_every = 1
  )
  expect_equal(
    nrow(lda_every$log_likelihood),
    tail(lda_every$log_likelihood$iteration, 1) + 1
  )
  
  # while we're here... check dimensions and names of objects
  expect_s3_class(lda, "tidylda")

  expect_length(lda$alpha, 1)
  
  expect_length(lda$eta, 1)
  
  expect_equal(sum(dim(lda$beta) == c(4, ncol(d1))), 2)

  expect_equal(sum(dim(lda$beta) == dim(lda$lambda)), 2)

  expect_equal(sum(dim(lda$theta) == c(nrow(d1), nrow(lda$beta))), 2)

  expect_setequal(colnames(lda$beta), colnames(d1))

  expect_setequal(rownames(lda$beta), colnames(lda$theta))

  expect_setequal(rownames(lda$theta), rownames(d1))

  # scalar priors optimizing alpha
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    eta = 0.05,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_length(lda$alpha, 4)

  # vector priors
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20,
    burnin = 10,
    alpha = rep(0.1, 4),
    eta = rep(0.05, ncol(d1)),
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_length(lda$alpha, 4)

  expect_length(lda$eta, ncol(d1))

  # eta as matrix prior
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    eta = matrix(0.05, nrow = 4, ncol = ncol(d1)),
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_true(inherits(lda$eta, "matrix"))
})

test_that("sparse priors for eta don't cause underflow failures", {
  m <- tidylda(
    data = nih_sample_dtm,
    k = 10,
    iterations = 20,
    burnin = 15,
    alpha = 0.05,
    eta = 0.01,
    optimize_alpha = FALSE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(m, "tidylda")
})

test_that("errors hit for malformed parameters", {

  # k = 1 is bad
  expect_error(
    tidylda(
      data = d1,
      k = 1,
      iterations = 20, burnin = 10,
      alpha = 0.1, eta = 0.05,
      optimize_alpha = TRUE,
      calc_likelihood = FALSE,
      calc_r2 = FALSE,
      return_data = FALSE,
      verbose = FALSE
    ),
    regexp = "k must be 2 or greater"
  )

  # iterations not specified
  expect_error(
    tidylda(
      data = d1,
      k = 10,
      alpha = 0.1, eta = 0.05,
      optimize_alpha = TRUE,
      calc_likelihood = FALSE,
      calc_r2 = FALSE,
      return_data = FALSE,
      verbose = FALSE
    ),
    label = "iterations not specified"
  )
  

  # burnin >= iterations
  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 21,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  # non-numeric k
  expect_error(tidylda(
    data = d1,
    k = "4",
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  # iterations not specified
  expect_error(tidylda(
    data = d1,
    k = 4
  ))

  # non-logical logicals
  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = "FALSE",
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = "FALSE",
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = "FALSE",
    return_data = FALSE
  ))

  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = "FALSE"
  ))
  
  # These two used to assert that `threads` above nrow(dtm) is an error and that
  # fewer than 100 documents per thread warns of "a poor fit". Both described the
  # abandoned batched Gibbs implementation, where the partition genuinely changed
  # the model. The warpLDA engine seeds every work item from its own index
  # (D12), so the result does not depend on the thread count at all, and the word
  # pass parallelizes over the vocabulary rather than over documents. Neither
  # condition is a problem now, and neither should complain.
  expect_no_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    threads = nrow(d1) + 1,
    verbose = FALSE
  ))

  expect_no_warning(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    threads = 2,
    verbose = FALSE
  ))
  
  # data doesn't have column names
  d3 <- d1
  colnames(d3) <- NULL
  expect_error(
    tidylda(
      data = d3,
      k = 4,
      iterations = 20
    )
  )
  
})

# note as of this writing, not parallel, 
# but use of threads argument should not throw errors
test_that("parallelism works as expected", {
  suppressWarnings(
    lda <- tidylda(
      data = d1,
      k = 4,
      iterations = 20, burnin = 10,
      alpha = 0.1, eta = 0.05,
      optimize_alpha = FALSE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE,
      threads = 2,
      verbose = FALSE
    )
  )

  
  expect_s3_class(lda, "tidylda")
  
})
