context("test utility functions")

library(tidytext, quietly = TRUE, verbose = FALSE)

docs <- nih_sample

tidy_docs <- unnest_tokens(
  tbl = docs[, c("APPLICATION_ID", "ABSTRACT_TEXT")],
  input = "ABSTRACT_TEXT",
  output = "word"
)

tidy_docs$count <- 1

triplet_dtm <- cast_dtm(tidy_docs,
  document = "APPLICATION_ID",
  term = "word",
  value = "count"
)

q_dfm <- cast_dfm(tidy_docs,
  document = "APPLICATION_ID",
  term = "word",
  value = "count"
)

mat <- as.matrix(q_dfm)

vec <- mat[1, ]

vec_nonames <- vec

names(vec_nonames) <- NULL

### tests for convert_dtm ----

test_that("convert_dtm can handle various inputs", {

  # simple triplet
  d <- convert_dtm(triplet_dtm)

  expect_true(inherits(d, "dgCMatrix"))

  expect_equal(nrow(d), triplet_dtm$nrow)

  expect_equal(ncol(d), triplet_dtm$ncol)

  expect_equivalent(colnames(d), triplet_dtm$dimnames$Terms)

  expect_equivalent(rownames(d), triplet_dtm$dimnames$Docs)


  # dfm
  d <- convert_dtm(q_dfm)

  expect_true(inherits(d, "dgCMatrix"))

  expect_equal(nrow(d), nrow(q_dfm))

  expect_equal(ncol(d), ncol(q_dfm))

  expect_equivalent(colnames(d), colnames(q_dfm))

  expect_equivalent(rownames(d), rownames(q_dfm))


  # dense matrix
  d <- convert_dtm(mat)

  expect_true(inherits(d, "dgCMatrix"))

  expect_equal(nrow(d), nrow(mat))

  expect_equal(ncol(d), ncol(mat))

  expect_equivalent(colnames(d), colnames(mat))

  expect_equivalent(rownames(d), rownames(mat))

  # vector with names
  d <- convert_dtm(vec)

  expect_true(inherits(d, "dgCMatrix"))

  expect_equal(nrow(d), 1)

  expect_equal(ncol(d), length(vec))

  expect_equivalent(colnames(d), names(vec))

  # vector without names
  expect_error(convert_dtm(vec_nonames))

  # not a supported class
  expect_error(convert_dtm(list(a = vec)))
})

### tests for format_eta and format_alpha ----

# since all variations (I think) that don't throw an error are tested in
# test-lda_core.R, this just tests bad inputs

test_that("format_eta chokes on bad inputs", {

  # eta non numeric
  expect_error(format_eta(eta = "WRONG!", k = 3, Nv = 10))

  # eta has na values
  expect_error(format_eta(eta = NA, k = 3, Nv = 10))

  # eta is zero
  expect_error(format_eta(eta = 0, k = 3, Nv = 10))

  # eta doesn't conform to vocabulary or topics
  expect_error(format_eta(eta = numeric(5) + 3, k = 3, Nv = 10))

  expect_error(format_eta(eta = matrix(1, nrow = 2, ncol = 10), k = 3, Nv = 10))

  # eta is a completely unsupported type
  expect_error(format_eta(eta = list(numeric(10) + 3), k = 3, Nv = 10))
})

test_that("format_alpha also chokes on bad inputs", {

  # alpha non numeric
  expect_error(format_alpha(alpha = "WRONG!", k = 3))

  # alpha has na values
  expect_error(format_alpha(alpha = NA, k = 3))

  # alpha is zero
  expect_error(format_alpha(alpha = 0, k = 3))

  # alpha doesn't conform to vocabulary or topics
  expect_error(format_alpha(alpha = numeric(5) + 3, k = 3))
})

test_that("tidy_dgcmatrix works as expected",{
  
  d <- convert_dtm(triplet_dtm)
  
  tmat <- tidy_dgcmatrix(d)
  
  expect_equal(sum(colnames(d) %in% tmat$term), ncol(d))
  
  expect_equal(nrow(tmat), sum(d > 0))
  
})

test_that("lambda works as expected",{
  dtm <- nih_sample_dtm
  
  d1 <- dtm[1:50, ]
  
  # make sure we have different vocabulary for each data set
  d1 <- d1[, Matrix::colSums(d1) > 0]
  
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20,
    verbose = FALSE
  )
  
  # proper function
  l <- 
    tidylda:::calc_lambda(
      beta = lda$beta,
      theta = lda$theta,
      p_docs = Matrix::rowSums(d1),
      correct = TRUE
    )
  
  expect_true(inherits(l, "matrix"))
  
  # p_docs is null
  l <- 
    tidylda:::calc_lambda(
      beta = lda$beta,
      theta = lda$theta,
      p_docs = NULL,
      correct = TRUE
    )
  
  expect_true(inherits(l, "matrix"))
  
  # p_docs contains NA values
  p <- Matrix::rowSums(d1)
  
  p[5] <- NA
  
  expect_warning(
    tidylda:::calc_lambda(
      beta = lda$beta,
      theta = lda$theta,
      p_docs = p,
      correct = TRUE
    )
  )
  
  
})
