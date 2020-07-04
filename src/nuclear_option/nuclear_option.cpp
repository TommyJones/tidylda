// Functions to make a collapsed gibbs sampler for LDA

// Export this as a header for use in other packages
// [[Rcpp::interfaces(r, cpp)]] 

// #include "sample_int.h"
#include "matrix_conversions.h"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define ARMA_64BIT_WORD

#include <Rmath.h>
#include <R.h>
using namespace Rcpp;

#include <cmath>

// declare a sampling function
// will remove later
arma::uword samp_one(const arma::vec &pvec) {
  
  // check that all elements are positive
  if (arma::any(pvec < 0.0))
    Rcpp::stop("probabilities have to be positive");
  
  // generate a vector of indices, so that 0 represents the largest 
  // element, 1 the secon largest, and so on
  arma::uvec indx = arma::sort_index(pvec, "descend");
  
  // generate the q-vector, the vector of cumulative sums
  arma::vec qvec = arma::cumsum(arma::sort(pvec, "descend"));
  
  // draw randomly from (0,s)
  double u = arma::conv_to<double>::from(qvec.tail(1L));
  u *= arma::randu(); // note that this will not respect R seed! Must fix!
  
  // find interval into which u falls
  for (arma::uword k = 0; k < pvec.n_elem; ++k) {
    
    if (u <= qvec(k))
      return indx(k);
    
  }
  
  Rcpp::stop("couldn't find index (samp_one)");
  
}




//' Main C++ Gibbs sampler for Latent Dirichlet Allocationocation
//' @keywords internal
//' @description
//'   This is the C++ Gibbs sampler for LDA. "Abandon Allocation hope, ye who enter here."
//' @param Docs List with one element for each document and one entry for each token
//'   as formatted by \code{\link[tidylda]{initialize_topic_counts}}
//' @param Nk int number of topics
//' @param beta NumericMatrix for prior of tokens over topics
//' @param alpha NumericVector prior for topics over documents
//' @param Cd IntegerMatrix denoting counts of topics in documents
//' @param Cv IntegerMatrix denoting counts of tokens in topics
//' @param Ck IntegerVector denoting counts of topics across Allocation tokens
//' @param Zd List with one element for each document and one entry for each token
//'   as formatted by \code{\link[tidylda]{initialize_topic_counts}}
//' @param Phi NumericMatrix denoting probability of tokens in topics
//' @param iterations int number of gibbs iterations to run in total
//' @param burnin int number of burn in iterations
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
//' @param calc_likelihood bool do you want to calculate the log likelihood each
//'   iteration?
//' @param optimize_alpha bool do you want to optimize alpha each iteration?
//'
// [[Rcpp::export]]
Rcpp::List fit_lda_c(
    const std::vector<std::vector<std::size_t>>&  Docs, 
    const std::vector<std::vector<std::size_t>>&  Zd_in,
    const IntegerMatrix&                          Cd_in,
    const IntegerMatrix&                          Cv_in,
    const std::vector<std::size_t>&               Ck_in,
    const std::vector<double>                     alpha_in, 
    const NumericMatrix&                          beta_in, 
    const std::size_t&                            iterations, 
    const int&                                    burnin,
    const bool&                                   optimize_alpha, 
    const bool&                                   calc_likelihood,
    const NumericMatrix&                          Phi_in,
    const bool&                                   freeze_topics
) {
  
  // ***********************************************************************
  // Convert input matrices to 2-dimensional std::vector for speed
  // Also convert some other objects to avoid overwriting inputs
  // ***********************************************************************
  auto Zd = Zd_in;
  
  auto Cd = mat_to_vec(Cd_in);
  auto Cv = mat_to_vec(Cv_in, true);
  auto Ck = Ck_in;
  
  auto alpha = alpha_in;
  auto beta = mat_to_vec(beta_in, true);
  
  auto Phi = mat_to_vec(Phi_in, true);
  
  // ***********************************************************************
  // Variables and other set up for the main sampler
  // ***********************************************************************
  
  const std::size_t Nv = Cv[0].size();
  const std::size_t Nd = Cd.size();
  const std::size_t Nk = Cv.size();
  
  auto sum_tokens = std::accumulate(Ck.begin(), Ck.end(), 0);
  auto sum_alpha = std::accumulate(alpha.begin(), alpha.end(), 0.0); 
  auto sum_beta = std::accumulate(beta[0].begin(), beta[0].end(), 0.0);
  
  std::vector<double> qz(Nk); // placehodler for probability of topic
  
  arma::uword z; // placeholder for topic sampled
  
  // For aggregating samples post burn in
  std::vector<std::vector<std::size_t>> Cd_sum(Nd);
  std::vector<std::vector<double>> Cd_mean(Nd);
  
  for (std::size_t d = 0; d < Nd; d++) {
    for (std::size_t k = 0; k < Nk; k++) {
      Cd_sum[d].push_back(0);
      Cd_mean[d].push_back(0.0);
    }
  }
  
  std::vector<std::vector<std::size_t>> Cv_sum(Nk);
  std::vector<std::vector<double>> Cv_mean(Nk);
  
  for (std::size_t k = 0; k < Nk; k++) {
    for (std::size_t v = 0; v < Nv; v++) {
      Cv_sum[k].push_back(0);
      Cv_mean[k].push_back(0.0);
    }
  }
  
  // ***********************************************************************
  // Variables and other set up for the log likelihood function
  // ***********************************************************************
  std::vector<std::vector<double>> log_likelihood(iterations);
  
  double lgbeta(0.0); // if calc_likelihood, we need this term
  
  double lgalpha(0.0); // if calc_likelihood, we need this term
  
  if (calc_likelihood) { // if calc_likelihood, actuAllocationy populate this stuff
    
    for (std::size_t n = 0; n < Nv; n++) {
      lgbeta += lgamma(beta[0][n]);
    }
    
    lgbeta = (lgbeta - lgamma(sum_beta)) * Nk; 
    
    for (std::size_t k = 0; k < Nk; k++) {
      lgalpha += lgamma(alpha[k]);
    }
    
    lgalpha = (lgalpha - lgamma(sum_alpha)) * Nd;
  }
  
  
  // ***********************************************************************
  // Main Gibbs iterations
  // ***********************************************************************
  for (std::size_t t = 0; t < iterations; t++) { // for each iteration
    
    for (std::size_t d = 0; d < Nd; d++) { // for each document
      
      R_CheckUserInterrupt();
      
      // IntegerVector doc = Docs[d]; // placeholder for a document
      // IntegerVector zd = Zd[d]; // placeholder for doc-word-topic assigment
      auto doc = Docs[d]; // placeholder for a document
      auto zd = Zd[d]; // placeholder for doc-word-topic assigment
      
      for (std::size_t n = 0; n < doc.size(); n++) { // for each word in that document
        
        // handle things differently based on freeze_topics
        // note some copied code, but minimizes number of checks in this loop
        if (! freeze_topics) {
          
          // discount for the n-th word with topic z
          Cd[d][zd[n]]--; 
          Cv[zd[n]][doc[n]]--;
          Ck[zd[n]]--;
          
          // calculate probability vector
          for (std::size_t k = 0; k < Nk; k++) {
            qz[k] = (Cv[k][doc[n]] + beta[k][doc[n]]) / 
              (Ck[k] + sum_beta) *
              (Cd[d][k] + alpha[k]) / 
              (doc.size() + sum_alpha);
          }
          
          // sample topic
          z = samp_one(qz);
          
          // update counts based on sampled topic
          zd[n] = z;
          
          Cd[d][zd[n]]++; 
          Cv[zd[n]][doc[n]]++; 
          Ck[zd[n]]++; 
          
        } else {
          
          // discount for the n-th word with topic z
          Cd[d][zd[n]]--; 

          // calculate probability vector
          for (std::size_t k = 0; k < Nk; k++) {
            qz[k] = Phi[k][doc[n]] *
              (Cd[d][k] + alpha[k]) / 
              (doc.size() + sum_alpha);
          }
          
          // sample topic
          z = samp_one(qz);
          
          // update counts based on sampled topic
          zd[n] = z;
          
          Cd[d][zd[n]]++; 
          
        } // end if
      } // end loop over tokens
      
      Zd[d] = zd; // update vector of sampled topics
      
    } // end loop over documents
    
    // if using burnin, update sums
    if (burnin > -1 && t >= burnin) {
      
      for (std::size_t d = 0; d < Nd; d++) {
        for (std::size_t k = 0; k < Nk; k++) {
          Cd_sum[d][k] += Cd[d][k];
        }
      }
      
      if (! freeze_topics) {
        for (std::size_t k = 0; k < Nk; k++) {
          for (std::size_t v = 0; v < Nv; v++) {
            Cv_sum[k][v] += Cv[k][v];
          }
        }
      }
    }
    
    // if calculating log likelihood, do so 
    if (calc_likelihood && !freeze_topics) {
      
      std::vector<double> tmp(3);
      
      // get phi probability matrix @ this iteration
      std::vector<std::vector<double>> phi_prob(Nk);

      double denom(0.0);
      
      double lp_beta(0.0); // log probability of beta prior
      
      for (std::size_t k = 0; k < Nk; k++) {
        
        std::vector<double> tmp(Nv);
        
        phi_prob[k] = tmp;
        
        // get the denominator
        for (std::size_t v = 0; v < Nv; v++) {
          denom += Cv[k][v] + beta[k][v];
        }
        
        // get the probability
        for (std::size_t v = 0; v < Nv; v++) {
          phi_prob[k][v] = ((double)Cv[k][v] + beta[k][v]) / denom;
          
          lp_beta += (beta[k][v] - 1) * log(phi_prob[k][v]);
        }
      }
      
      lp_beta += lgbeta;
      
      // for each document, get the log probability of the words
      
      double lp_alpha(0.0); // log probability of alpha prior
      
      double lpd(0.0); // log probability of documents
      
      for (std::size_t d = 0; d < Nd; d++) {
        
        std::vector<double> theta_prob(Nk); // probability of each topic in this document
        
        auto doc = Docs[d];
        
        NumericVector lp(doc.size()); // log probability of each word under the model
        
        double denom(0.0);
        
        for (std::size_t k = 0; k < Nk; k++) {
          denom += (double)Cd[d][k] + alpha[k];
        }
        
        for (std::size_t k = 0; k < Nk; k++) {
          theta_prob[k] = ((double)Cd[d][k] + alpha[k]) / denom;
          
          lp_alpha += (alpha[k] - 1) * log(theta_prob[k]);
        }
        
        for (std::size_t n = 0; n < doc.size(); n++) {
          
          lp[n] = 0.0;
          
          for (std::size_t k = 0; k < Nk; k++) {
            
            lp[n] += theta_prob[k] * phi_prob[k][doc[n]];
            
          }
          
          lpd += log(lp[n]);
        }
      }
      
      lp_alpha += lgalpha;
      

      tmp[0] = static_cast<double>(t);
      
      tmp[1] = lpd; // log probability of whole corpus under the model w/o priors
      
      tmp[2] = lpd + lp_alpha + lp_beta; // second likelihood calculation here
      
      log_likelihood[t] = tmp;
      
    }
    // if optimizing alpha, do so
    if (optimize_alpha) {

      for (std::size_t k = 0; k < Nk; k++) {
        alpha[k] = (static_cast<double>(Ck[k]) / static_cast<double>(sum_tokens)) * sum_alpha;
      }
    } 
    
  } // end iterations
  
  
  // ***********************************************************************
  // Cleanup and return values
  // ***********************************************************************
  
  // If using burn in iterations, 
  // average Cd and Cv across post burn in iterations
  if (burnin > -1) { 
    const double diff(iterations - burnin);
    
    // average over chain after burnin
    for (std::size_t d = 0; d < Nd; d++) {
      for (std::size_t k = 0; k < Nk; k++) {
        Cd_mean[d][k] = static_cast<double>(Cd_sum[d][k]) / diff;
      }
    }
    
    if (! freeze_topics) {
      for (std::size_t v = 0; v < Nv; v++) {
        for (std::size_t k = 0; k < Nk; k++) {
          Cv_mean[k][v] = static_cast<double>(Cv_sum[k][v]) / diff;
        }
      }
    }
  }
  
  // Return the final list ***
  return Rcpp::List::create(
    _["Cd"]             = vec_to_mat(Cd, true),
    _["Cv"]             = vec_to_mat(Cv, true),
    _["Ck"]             = Ck,
    _["Cd_mean"]        = vec_to_mat(Cd_mean, true),
    _["Cv_mean"]        = vec_to_mat(Cv_mean, true),
    _["Cd_sum"]         = vec_to_mat(Cd_sum, true),
    _["Cv_sum"]         = vec_to_mat(Cv_sum, true),
    _["log_likelihood"] = vec_to_mat(log_likelihood),
    _["alpha"]          = alpha,
    _["beta"]           = beta_in // does not get modified, so don't waste compute converting beta
  );
}

/*** R
library(tidyverse)

library(testthat)

dtm <- textmineR::nih_sample_dtm

k <- 10

alpha <- rep(0.1, k)

beta <- matrix(0.05, nrow = k, ncol = ncol(dtm))

counts <- 
  tidylda:::initialize_topic_counts(
    dtm = dtm, 
    k = 10,
    alpha = rep(0.1, 10), 
    beta = matrix(0.05, nrow = 10, ncol = ncol(dtm)),
    threads = 1
  )

microbenchmark::microbenchmark({
  fit_lda_c(
    Docs = counts$docs,
    Zd_in = counts$Zd,
    beta_in = beta,
    alpha_in = alpha,
    Cd_in = counts$Cd,
    Cv_in = counts$Cv,
    Ck_in = counts$Ck,
    Phi = counts$Cv, # ignored
    iterations = 200,
    burnin = 175,
    freeze_topics = FALSE,
    calc_likelihood = FALSE,
    optimize_alpha = FALSE
  )},
  times = 20
)

### tests below to ensure correct computations
m <- fit_lda_c(
  Docs = counts$docs,
  Zd_in = counts$Zd,
  beta_in = beta,
  alpha_in = alpha,
  Cd_in = counts$Cd,
  Cv_in = counts$Cv,
  Ck_in = counts$Ck,
  Phi = counts$Cv, # ignored
  iterations = 200,
  burnin = 175,
  freeze_topics = FALSE,
  calc_likelihood = TRUE,
  optimize_alpha = TRUE
)

test_that("average coherence for Cv is greater than 0.1",{
  p <- m$Cv
  
  colnames(p) <- colnames(dtm)
  
  rownames(p) <- 1:k
  
  expect_true(mean(p) >= 0.1)
})

test_that("average coherence for Cv_mean is greater than 0.1",{
  p <- m$Cv_mean
  
  colnames(p) <- colnames(dtm)
  
  rownames(p) <- 1:k
  
  summary(textmineR::CalcProbCoherence(p, dtm))
})

test_that("checksums match expectation",{
  
  sum_tokens <- sum(dtm)
  
  expect_equal(sum(m$Cd), sum_tokens)
  
  expect_equal(sum(m$Cv), sum_tokens)
  
  expect_equal(sum(m$Cd_mean), sum_tokens)

  expect_equal(sum(m$Cv_mean), sum_tokens)
  
  
})


*/