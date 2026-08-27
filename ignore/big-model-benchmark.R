library(tidyverse)
library(tidylda)
library(Matrix)

dat <- read_csv("../dissertation/data-raw/RePORTER_PRJABS_C_FY2014.csv")

dat <- dat |> 
  mutate(ABSTRACT_TEXT = iconv(ABSTRACT_TEXT, from = "UTF-8", to = "ASCII", sub = ""))

dat <- dat |> 
  mutate(dat, sample1 = APPLICATION_ID %in% sample(dat$APPLICATION_ID, 20000))

dtm <- textmineR::CreateDtm(
  doc_vec = dat$ABSTRACT_TEXT,
  doc_names = dat$APPLICATION_ID
)

dtm1 <- dtm[dat$sample1, ]

dtm1 <- dtm1[, colSums(dtm1 > 0) >= 5]

dim(dtm1)

lda <- tidylda(
  data = dtm1,
  k = 300,
  iterations = 250,
  burnin = 200,
  eta = colSums(dtm1) / sum(dtm1) * 250,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  threads = parallel::detectCores() - 1,
  mh_steps = 4
)

lda

dtm2 <- dtm[ ! dat$sample1, ]

dtm2 <- dtm2[, colSums(dtm2 > 0) >= 5]

dim(dtm2)

preds <- predict(
  lda,
  new_data = dtm2,
  method = "mh",
  iterations = 250,
  burnin = 200,
  threads = parallel::detectCores() - 1,
  mh_steps = 4
)

lda_refit <- refit(
  lda,
  new_data = dtm2,
  iterations = 250,
  burnin = 200,
  prior_weight = 0.8,
  additional_k = 0,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  mh_steps = 4,
  threads = parallel::detectCores() - 1
)

lda_refit

system.time({
  lda2 <- tidylda(
    data = dtm2,
    k = 300,
    eta = colSums(dtm2) / sum(dtm2) * 250,
    iterations = 250,
    burnin = 200,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    mh_steps = 4,
    threads = parallel::detectCores() - 1
  )
})
