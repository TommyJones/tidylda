context("tests of print method for tidylda")

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


### Tests for the print method ----

test_that("print.tidylda behaves as expected", {
  
  # no error
  capture.output(print(lda))
  
  # assignment creates a new object of class tidylda
  capture.output(m <- print(lda))
  
  expect_true("tidylda" %in% class(m))
  
  expect_named(m, names(lda))
  
  # can modify digits
  m2 <- tidylda(
    data = d1, 
    k = 5, 
    iterations = 20, 
    calc_r2 = TRUE, 
    verbose = FALSE
  )
  
  capture.output(print(m2, digits = 2))
})
