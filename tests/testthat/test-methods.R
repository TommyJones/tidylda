context("tests of methods for tidylda")

dtm <- textmineR::nih_sample_dtm

d1 <- dtm[1:50, ]

d2 <- dtm[51:100, ]

# make sure we have different vocabulary for each data set
d1 <- d1[, Matrix::colSums(d1) > 0]

d2 <- d2[, Matrix::colSums(d2) > 0]

lda <- tidylda(dtm = d1, 
               k = 4, 
               iterations = 20, burnin = 10,
               alpha = 0.1, beta = 0.05,
               optimize_alpha = TRUE,
               calc_likelihood = FALSE,
               calc_r2 = FALSE,
               return_data = FALSE)

### Tests for predictions ----

test_that("can make predictions without error",{
  # one row gibbs with burnin
  p <- predict(object = lda, newdata = d2[1,], method = "gibbs", iterations = 20, burnin = 10)
  
  expect_equal(nrow(p), 1)
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # multi-row gibbs with burnin
  p <- predict(object = lda, newdata = d2, method = "gibbs", iterations = 20, burnin = 10)
  
  expect_equal(nrow(p), nrow(d2))
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # single row dot method
  p <- predict(object = lda, newdata = d2[1,], method = "dot")
  
  expect_equal(nrow(p), 1)
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
  # multi-row dot method
  p <- predict(object = lda, newdata = d2, method = "dot")
  
  expect_equal(nrow(p), nrow(d2))
  
  expect_equal(ncol(p), ncol(lda$theta))
  
  expect_setequal(colnames(p), colnames(lda$theta))
  
})


### Tests for updates ----

test_that("can update models",{
  # continuation of the old model for another 20 iterations
  # matters because dtm lines up exactly with existing vocabulary etc.
  lda2 <- update(object = lda, 
                 dtm = d1,
                 additional_k = 0,
                 phi_as_prior = FALSE,
                 iterations = 20, 
                 burnin = 10,
                 optimize_alpha = TRUE,
                 calc_likelihood = FALSE,
                 calc_r2 = FALSE,
                 return_data = FALSE)
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d1))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta))
  
  expect_equal(sum(dim(lda2$phi) == dim(lda$phi)), 2)
  
  expect_equal(sum(dim(lda2$gamma) == dim(lda$gamma)), 2)
  
  
  # new data adding no extra topics no phi as prior
  lda2 <- update(object = lda, 
                 dtm = d2,
                 additional_k = 0,
                 phi_as_prior = FALSE,
                 iterations = 20, 
                 burnin = 10,
                 optimize_alpha = TRUE,
                 calc_likelihood = FALSE,
                 calc_r2 = FALSE,
                 return_data = FALSE)
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta))
  
  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))
  
  # 1 additonal topic and no phi as prior
  lda2 <- update(object = lda, 
                 dtm = d2,
                 additional_k = 1,
                 phi_as_prior = FALSE,
                 iterations = 20, 
                 burnin = 10,
                 optimize_alpha = TRUE,
                 calc_likelihood = FALSE,
                 calc_r2 = FALSE,
                 return_data = FALSE)
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 1)
  
  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))
  
  # 3 additional topics and no phi as prior
  lda2 <- update(object = lda, 
                 dtm = d2,
                 additional_k = 3,
                 phi_as_prior = FALSE,
                 iterations = 20, 
                 burnin = 10,
                 optimize_alpha = TRUE,
                 calc_likelihood = FALSE,
                 calc_r2 = FALSE,
                 return_data = FALSE)
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 3)
  
  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))
  
  # no additional topics and phi as prior
  lda2 <- update(object = lda, 
                 dtm = d2,
                 additional_k = 0,
                 phi_as_prior = TRUE,
                 iterations = 20, 
                 burnin = 10,
                 optimize_alpha = TRUE,
                 calc_likelihood = FALSE,
                 calc_r2 = FALSE,
                 return_data = FALSE)
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta))
  
  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))
  
  
  # 3 additonal topics and phi as prior
  lda2 <- update(object = lda, 
                 dtm = d2,
                 additional_k = 3,
                 phi_as_prior = TRUE,
                 iterations = 20, 
                 burnin = 10,
                 optimize_alpha = TRUE,
                 calc_likelihood = FALSE,
                 calc_r2 = FALSE,
                 return_data = FALSE)
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 3)
  
  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))
  
  
})


### Tests for the print method ----

test_that("print.tidylda behaves as expected",{
  
  # no error
  print(lda)
  
  # assignment creates a new object of class tidylda
  m <- print(lda)
  
  expect_true("tidylda" %in% class(m))
  
  expect_named(m, names(lda))
  
  # can modify digits
  m2 <- tidylda(dtm = d1, k = 5, iterations = 20, calc_r2 = TRUE)
  
  print(m2, digits = 2)
  
})