context("test utility functions")

library(tidytext, quietly = TRUE, verbose = FALSE)

docs <- nih_sample

tidy_docs <- unnest_tokens(
  tbl = docs[, c("APPLICATION_ID", "ABSTRACT_TEXT")],
  input = "ABSTRACT_TEXT",
  output = "word"
)

tidy_docs$count <- 1



### tests for convert_dtm ----

test_that("convert_dtm can handle various inputs", {

  skip_if_not_installed("tm")
  skip_if_not_installed("quanteda")

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
  
  
  
  sparse_mat <- cast_sparse(tidy_docs,
                            row = "APPLICATION_ID",
                            column = "word",
                            value = "count"
  )
  
  mat <- as.matrix(q_dfm)
  
  vec <- mat[1, ]
  
  vec_nonames <- vec
  
  names(vec_nonames) <- NULL
  
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
  
  # sparse matrix from Matrix library
  d <- convert_dtm(sparse_mat)
  
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
  
  skip_if_not_installed("tm")
  
  triplet_dtm <- cast_dtm(tidy_docs,
                          document = "APPLICATION_ID",
                          term = "word",
                          value = "count"
  )
  
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


context("memory guards and allocation helpers")

test_that("check_result_size fires only above the ceiling", {
  # Well under the 1 GB default: silent.
  expect_silent(check_result_size(1e3, 4, "thing", "do something else"))

  # Over it: an error naming the size, the limit and the remedy.
  expect_error(
    check_result_size(1e9, 4, "thing", "do something else"),
    "thing would need about"
  )
  expect_error(check_result_size(1e9, 4, "thing", "do something else"), "GB limit")
  expect_error(check_result_size(1e9, 4, "thing", "try a slice"), "try a slice")
  expect_error(
    check_result_size(1e9, 4, "thing", "x"),
    "tidylda.max_result_size"
  )
})


test_that("check_result_size respects the option", {
  on.exit(options(tidylda.max_result_size = NULL), add = TRUE)

  # A request that is fine by default becomes an error under a tiny ceiling.
  expect_silent(check_result_size(1e5, 4, "thing", "x"))

  options(tidylda.max_result_size = 1000)
  expect_error(check_result_size(1e5, 4, "thing", "x"), "above the")

  # And raising it lets the same request through.
  options(tidylda.max_result_size = 1e12)
  expect_silent(check_result_size(1e5, 4, "thing", "x"))
})


test_that("session_memory_limit returns bytes or NA, never nonsense", {
  x <- session_memory_limit()

  expect_length(x, 1)
  expect_true(is.na(x) || (is.numeric(x) && x > 0))
})


test_that("posterior and tidy refuse impossible requests", {
  on.exit(options(tidylda.max_result_size = NULL), add = TRUE)

  d <- nih_sample_dtm[1:30, ]
  set.seed(1)
  m <- tidylda(d, k = 4, iterations = 15, burnin = 5,
               calc_r2 = FALSE, verbose = FALSE)

  # Both work normally.
  expect_s3_class(tidy(m, "beta"), "data.frame")
  expect_s3_class(posterior(m, "beta", which = 1, times = 3), "tbl_df")

  options(tidylda.max_result_size = 1000)

  expect_error(tidy(m, "beta"), "tidy\\(matrix")
  expect_error(posterior(m, "beta", which = 1, times = 3), "posterior\\(matrix")

  # The remedy each message suggests must actually work. tidy() dispatches on a
  # matrix, and slicing keeps the original topic labels.
  options(tidylda.max_result_size = NULL)
  sliced <- tidy(m$beta[3:4, ], "beta")
  expect_setequal(unique(sliced$topic), c(3, 4))
})
