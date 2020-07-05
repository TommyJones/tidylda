library(tidyverse)


tcm <- textmineR::CreateTcm(
  doc_vec = c(textmineR::nih_sample$ABSTRACT_TEXT, textmineR::nih_sample$PROJECT_TITLE),
  skipgram_window = 10
)

tcm <- tcm + Matrix::t(tcm)

vocab <- Matrix::colSums(tcm)

vocab <- names(vocab[vocab > 5])
  
tcm <- tcm[vocab, vocab]

tidylda_benchmark <-
microbenchmark::microbenchmark(
  tidylda::tidylda(
    dtm = tcm,
    k = 10,
    iterations = 200,
    burnin = 175
  ),
  times = 10
)

textminer_benchmark <- 
microbenchmark::microbenchmark(
  textmineR::FitLdaModel(
    dtm = tcm,
    k = 10,
    iterations = 200,
    burnin = 175,
    cpus = 1
  ),
  times = 10
)


