library(tidyverse)
library(tidylda)
library(Matrix)
library(RSpectra)   # new dep; irlba is a drop-in

corpus_hypers <- function(dtm, alpha_0 = 1, k_max = 400, min_count = 10) {
  dtm <- dtm[, Matrix::colSums(dtm) >= min_count, drop = FALSE]
  n_d <- Matrix::rowSums(dtm)
  dtm <- dtm[n_d > 0, , drop = FALSE]; n_d <- n_d[n_d > 0]
  m   <- Matrix::colSums(dtm) / sum(dtm)
  D   <- nrow(dtm); V <- ncol(dtm); N <- sum(n_d)
  
  # X^2 = sum x^2 / (N_d m_v) - N, over stored nonzeros only
  tri <- Matrix::summary(as(dtm, "TsparseMatrix"))     # i, j, x
  x2  <- sum(tri$x^2 / (n_d[tri$i] * m[tri$j])) - N
  phi <- (x2 / (V - 1) - D) / N
  
  A <- Matrix::Diagonal(x = 1/sqrt(n_d)) %*% dtm %*% Matrix::Diagonal(x = 1/sqrt(m))
  u <- sqrt(n_d); v <- sqrt(m)
  d_obs <- RSpectra::svds(
    A      = function(x, a) as.numeric(a$A %*% x - a$u * sum(a$v * x)),
    Atrans = function(y, a) as.numeric(Matrix::crossprod(a$A, y) - a$v * sum(a$u * y)),
    k = min(k_max, D - 1, V - 1), dim = c(D, V),
    args = list(A = A, u = u, v = v)
  )$d
  
  edge  <- sqrt(D) + sqrt(V)
  k_hat <- sum(d_obs > edge) + 1
  Pi    <- (1 - 1/k_hat) / phi
  
  list(
    k_hat     = k_hat,
    k_screen  = (x2 - (V - 1) * D) / edge^2,
    phi       = phi,
    Pi        = Pi,                              # what's actually identified
    alpha_max = Pi - 1,                          # feasibility bound on alpha_0
    eta_0     = if (alpha_0 + 1 < Pi) Pi / (alpha_0 + 1) - 1 else NA_real_,
    feasible  = alpha_0 + 1 < Pi,
    edge      = edge,
    sv        = d_obs
  )
}


dat <- read_csv("../dissertation/data-raw/RePORTER_PRJABS_C_FY2014.csv")

dat <- dat |> 
  mutate(ABSTRACT_TEXT = iconv(ABSTRACT_TEXT, from = "UTF-8", to = "ASCII", sub = ""))

dtm <- textmineR::CreateDtm(
  doc_vec = dat$ABSTRACT_TEXT,
  doc_names = dat$APPLICATION_ID
)

dtm <- dtm[, colSums(dtm > 0) >= 5]

dtm <- dtm[rowSums(dtm) > 0, ]

hp <- corpus_hypers(
  dtm = dtm,
  alpha_k = 0.5,
  k_max = 600,
  min_count = 10
)
