# load libraries
library(tidyverse)
library(textmineR)
library(tidylda)

# load raw data
sbir <- read_csv("https://data.www.sbir.gov/awarddatapublic/award_data.csv")

colnames(sbir) <- 
  colnames(sbir) %>%
  tolower() %>%
  str_replace_all(" +", "_")

# add unique identifier
sbir$sbir_id <- 1:nrow(sbir)

# pull out text columns
sbir_text <- sbir %>%
  select(
    sbir_id,
    award_title,
    abstract
  ) %>%
  mutate(
    award_title = str_conv(award_title, "UTF-8"),
    abstract = str_conv(abstract, "UTF-8")
  )

# Create a TCM with 10-degree skipgrams
# Creating titles and abstracts as separate documents so that titles are handled
# in isolation when constructing skipgrams
sbir_tcm <- CreateTcm(
  doc_vec = c(sbir_text$award_title, sbir_text$abstract),
  skipgram_window = 10, # arbitrary but standard
  stopword_vec = stopwords::stopwords("en"),
  verbose = TRUE
)

# completing a step that textmineR should've done
sbir_tcm <- sbir_tcm + t(sbir_tcm) 

# vocabulary pruning
sbir_tf <- TermDocFreq(sbir_tcm) %>%
  as_tibble()

vocab_keep <- sbir_tf$term[sbir_tf$doc_freq > 20]

sbir_tcm <- sbir_tcm[vocab_keep, vocab_keep]

save(sbir_tcm, file = "ignore/large-sparse-mat.RData")

# train an LDA model off of it
sbir_embedding <- tidylda(
  data = sbir_tcm,
  k = 10, 
  iterations = 20,
  burnin = 17,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  verbose = TRUE
)

