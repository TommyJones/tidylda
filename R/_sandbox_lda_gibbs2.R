# expriment with informed start to LDA
# this, among other things, enables transfer learning with LDA

source("R/initialize_topic_counts.R")

source("R/fit_lda_model.R")

source("R/predict.lda_topic_model.R")

source("R/update.lda_topic_model.R")

Rcpp::sourceCpp("src/lda_gibbs2.cpp")

library(textmineR)

# declare some global variables I'll need
dtm <- textmineR::nih_sample_dtm

dtm_pred_update <- dtm[81:100, ]

dtm_pred_update <- dtm_pred_update[, Matrix::colSums(dtm_pred_update) > 0 ]

dtm <- dtm[1:80, ]

dtm <- dtm[, Matrix::colSums(dtm) > 0]

k <- 10

alpha <- rep(0.05, k)

beta <- Matrix::colSums(dtm) / sum(dtm) * 50

# test mah function

lda <- fit_lda_model(dtm = dtm, 
                     k = k, 
                     iterations = 100, burnin = 50,
                     alpha = alpha, beta = beta,
                     optimize_alpha = TRUE,
                     calc_likelihood = TRUE,
                     calc_coherence = TRUE,
                     calc_r2 = TRUE,
                     return_data = F)

### testing for predict function ----

# one row gibbs with burnin
p <- predict(object = lda, newdata = dtm_pred_update[1,], method = "gibbs", iterations = 100, burnin = 75)

# one row gibbs without burnin
p <- predict(object = lda, newdata = dtm_pred_update[1,], method = "gibbs", iterations = 100)

# multi-row gibbs with burnin
p <- predict(object = lda, newdata = dtm_pred_update, method = "gibbs", iterations = 100, burnin = 75)

# multi-row gibbs without burnin
p <- predict(object = lda, newdata = dtm_pred_update, method = "gibbs", iterations = 100)

# multi-row dot method
p <- predict(object = lda, newdata = dtm_pred_update, method = "dot")

# single row dot method
p <- predict(object = lda, newdata = dtm_pred_update[1,], method = "dot")

### test the update function ----

# continuation of the old model for another 200 iterations
lda2 <- update(object = lda, 
               dtm = dtm,
               additional_k = 0,
               phi_as_prior = FALSE,
               iterations = 200, 
               burnin = 175,
               optimize_alpha = TRUE,
               calc_likelihood = TRUE,
               calc_coherence = TRUE,
               calc_r2 = TRUE,
               return_data = FALSE)

# no additonal topics and no phi as prior
lda2 <- update(object = lda, 
               dtm = dtm_pred_update,
               additional_k = 0,
               phi_as_prior = FALSE,
               iterations = 100, 
               burnin = 75,
               optimize_alpha = TRUE,
               calc_likelihood = TRUE,
               calc_coherence = TRUE,
               calc_r2 = TRUE,
               return_data = FALSE)

# 1 additonal topic and no phi as prior
lda2 <- update(object = lda, 
               dtm = dtm_pred_update,
               additional_k = 1,
               phi_as_prior = FALSE,
               iterations = 100, 
               burnin = 75,
               optimize_alpha = TRUE,
               calc_likelihood = TRUE,
               calc_coherence = TRUE,
               calc_r2 = TRUE,
               return_data = FALSE)

# 3 additional topics and no phi as prior
lda2 <- update(object = lda, 
               dtm = dtm_pred_update,
               additional_k = 3,
               phi_as_prior = FALSE,
               iterations = 100, 
               burnin = 75,
               optimize_alpha = TRUE,
               calc_likelihood = TRUE,
               calc_coherence = TRUE,
               calc_r2 = TRUE,
               return_data = FALSE)

# no additional topics and phi as prior
lda2 <- update(object = lda, 
               dtm = dtm_pred_update,
               additional_k = 0,
               phi_as_prior = TRUE,
               iterations = 100, 
               burnin = 75,
               optimize_alpha = TRUE,
               calc_likelihood = TRUE,
               calc_coherence = TRUE,
               calc_r2 = TRUE,
               return_data = FALSE)

# 3 additonal topics and phi as prior
lda2 <- update(object = lda, 
               dtm = dtm_pred_update,
               additional_k = 3,
               phi_as_prior = TRUE,
               iterations = 100, 
               burnin = 75,
               optimize_alpha = TRUE,
               calc_likelihood = TRUE,
               calc_coherence = TRUE,
               calc_r2 = TRUE,
               return_data = FALSE)



