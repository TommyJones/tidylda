context("tests of predict method for tidylda")

dtm <- nih_sample_dtm

d1 <- dtm[1:50, ]

d2 <- dtm[51:100, ]

# make sure we have different vocabulary for each data set
d1 <- d1[, Matrix::colSums(d1) > 0]

d2 <- d2[, Matrix::colSums(d2) > 0]

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

### Tests for predictions ----

test_that("can make predictions without error", {
  # one row gibbs with burnin
  p <- predict(
    object = lda, 
    new_data = d2[1, ], 
    method = "mh", 
    iterations = 20, 
    burnin = 10,
    verbose = FALSE
  )
  
  expect_equal(nrow(p), 1)
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # multi-row gibbs with burnin
  p <- predict(
    object = lda, 
    new_data = d2, 
    method = "mh", 
    iterations = 20, 
    burnin = 10,
    verbose = FALSE
  )
  
  expect_equal(nrow(p), nrow(d2))
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # single row dot method
  p <- predict(object = lda, new_data = d2[1, ], method = "dot")
  
  expect_equal(nrow(p), 1)
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # multi-row dot method
  p <- predict(object = lda, new_data = d2, method = "dot")
  
  expect_equal(nrow(p), nrow(d2))
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # multi row parallel 
  # (no longer parallel, but checks that threads arg doesn't cause an error)
  p <- predict(
    object = lda, 
    new_data = d2, 
    method = "mh", 
    iterations = 20, 
    burnin = 10,
    threads = 2,
    verbose = FALSE
  )  
  
  expect_true(inherits(p, "matrix"))
  
  # single row class with dot
  p <- predict(
    object = lda,
    new_data = d2[1, ],
    type = "class",
    method = "dot"
  )
  
  expect_equal(length(p), 1)
  
  expect_true(inherits(p, "integer"))
  
  # multi row class with dot
  p <- predict(
    object = lda,
    new_data = d2,
    type = "class",
    method = "dot"
  )
  
  expect_equal(length(p), nrow(d2))
  
  expect_true(inherits(p, "integer"))
  
  # single row distribution with gibbs
  p <- predict(
    object = lda, 
    new_data = d2[1, ],
    type = "distribution",
    method = "mh", 
    iterations = 20, 
    burnin = 10,
    times = 10,
    threads = 1,
    verbose = FALSE
  )  
  
  expect_true(inherits(p, "tbl"))
  
  # multi row distribution with dot
  p <- predict(
    object = lda,
    new_data = d2,
    type = "distribution",
    method = "dot",
    times = 10
  )
  expect_true(inherits(p, "tbl"))
  
})



test_that("malformed args in predict throw errors", {
  
  # threads > nrow(dtm)
  expect_message(
    predict(
      object = lda, 
      new_data = d2, 
      method = "mh", 
      iterations = 20, 
      burnin = 10,
      threads = nrow(d2) + 2,
      verbose = FALSE
    ), label = "threads > nrow(dtm)"
  )
  # no iterations specified
  expect_error(
    predict(object = lda, new_data = d2, method = "mh")
  )
  
  # burnin >= iterations
  expect_error(
    predict(object = lda, new_data = d2, method = "mh", iterations = 5, burnin = 6)
  )
  
  # incorrect method
  expect_error(
    predict(object = lda, new_data = d2, method = "oops")
  )
  
  # no overlap in vocabulary throws warning on "dot" by default
  nd <- numeric(10)
  
  names(nd) <- seq_along(nd) # numbers means no vocab overlap
  
  expect_warning(
    predict(object = lda, new_data = nd, method = "dot")
  )
  # no overlap in vocabulary doesn't throw warning on "dot" if specified
  expect_message(
    predict(object = lda, new_data = nd, method = "dot", no_common_tokens = "zero")
  )
  
  # no overlap in vocabulary sets every topic to 1/k and no message or warning
  p <- predict(object = lda, new_data = nd, method = "dot", no_common_tokens = "uniform")
  
  expect_equal(mean(p), 1 / nrow(lda$beta))
  
  # no_common_tokens has illegal value
  expect_error(
    predict(object = lda, new_data = nd, method = "dot", no_common_tokens = "WRONG!")
  )
  
  # type misspecified
  expect_error(
    predict(object = lda, new_data = d2, type = "blah", method = "dot")
  )
  
  # times misspecified while type = "distribution"
  expect_error(
    predict(
      object = lda, 
      new_data = d2, 
      type = "distribution", 
      method = "dot",
      times = NA
    )
  )
  
  expect_error(
    predict(
      object = lda, 
      new_data = d2, 
      type = "distribution", 
      method = "dot",
      times = "yeet"
    )
  )
  
  expect_error(
    predict(
      object = lda, 
      new_data = d2, 
      type = "distribution", 
      method = "dot",
      times = 0
    )
  )
  
})



test_that('method = "gibbs" still works but is deprecated', {
  # The sampler stopped being collapsed Gibbs in 0.1.0. The string is public API
  # so it keeps working, but it names something that no longer exists. Warned
  # once per session, so this test asserts the mapping rather than the warning --
  # which may already have fired elsewhere in the run.
  p_old <- predict(lda, d2[1:3, ], method = "gibbs", iterations = 20, burnin = 5,
                   verbose = FALSE)
  p_new <- predict(lda, d2[1:3, ], method = "mh", iterations = 20, burnin = 5,
                   verbose = FALSE)

  expect_equal(dim(p_old), dim(p_new))
  expect_equal(unname(rowSums(p_old)), rep(1, 3))

  # And "mh" is the default, so ordinary callers never trip the deprecation.
  expect_equal(formals(tidylda:::predict.tidylda)$method[[2]], "mh")
})
