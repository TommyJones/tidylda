// Functions to make a collapsed gibbs sampler for LDA

// Export this as a header for use in other packages
// [[Rcpp::interfaces(r, cpp)]] 

#include "parallel_gibbs_utils.h"
#include "sample_int.h"
#include "matrix_conversions.h"

#include <RcppArmadillo.h>
#define ARMA_64BIT_WORD

#include <RcppThread.h>

//' Make a lexicon for looping over in the gibbs sampler
//' @keywords internal
//' @description
//'   One run of the Gibbs sampler and other magic to initialize some objects.
//'   Works in concert with \code{\link[tidylda]{initialize_topic_counts}}.
//' @param Cd IntegerMatrix denoting counts of topics in documents
//' @param Phi NumericMatrix denoting probability of words in topics
//' @param dtm arma::sp_mat document term matrix
//' @param alpha NumericVector prior for topics over documents
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
// [[Rcpp::export]]
Rcpp::List create_lexicon(arma::imat&      Cd,
                          const arma::mat& Phi,
                          arma::sp_mat&    dtm,
                          const arma::vec  alpha,
                          const bool       freeze_topics,
                          const int        threads) {
  // ***************************************************************************
  // Initialize some variables
  // ***************************************************************************
  
  dtm = dtm.t(); // transpose dtm to take advantage of column major & parallel
  Cd  = Cd.t();  // transpose to put columns up
  
  auto                    sum_alpha = sum(alpha);
  std::vector<arma::imat> docs(dtm.n_cols);
  std::vector<arma::imat> Zd(dtm.n_cols);
  
  auto Nk = Cd.n_rows;
  
  // ***************************************************************************
  // Go through each document and split it into a lexicon and then sample a
  // topic for each token within that document
  // ***************************************************************************
  RcppThread::parallelFor(
    0,
    dtm.n_cols,
    [&](unsigned int d) {
      arma::vec qz(Nk);
      
      // arma::ivec       topic_index = Rcpp::seq_len(Nk) - 1;
      std::vector<int> topic_index(Nk);
      std::iota(std::begin(topic_index), std::end(topic_index), 0);
      
      // make a temporary vector to hold token indices
      auto nd(0);
      
      for (auto v = 0; v < dtm.n_rows; ++v) {
        nd += dtm(v, d);
      }
      
      arma::ivec doc(nd);
      arma::ivec zd(nd);
      arma::ivec z(1);
      
      // fill in with token indices
      std::size_t j = 0; // index of doc, advances when we have non-zero entries
      
      for (auto v = 0; v < dtm.n_rows; ++v) {
        if (dtm(v, d) > 0) { // if non-zero, add elements to doc
          
          // calculate probability of topics based on initially-sampled Phi and Cd
          for (auto k = 0; k < Nk; ++k) {
            qz[k] = log(Phi(k, v)) + log(Cd(k, d) + alpha[k]) - log(nd + sum_alpha - 1);
          }
          
          std::size_t idx(j + dtm(v, d)); // where to stop the loop below
          
          while (j < idx) {
            doc[j] = v;
            z      = lsamp_one(qz); // sample a topic here
            zd[j]  = z[0];
            j++;
          }
        }
      }
      
      // fill in docs[d] with the matrix we made
      docs[d] = doc;
      
      Zd[d] = zd;
      
      RcppThread::checkUserInterrupt();
    },
    threads);
  
  // ***************************************************************************
  // Calculate Cd, Cv, and Ck from the sampled topics
  // ***************************************************************************
  arma::imat Cd_out(Nk, dtm.n_cols);
  
  Cd_out.fill(0);
  
  arma::ivec Ck(Nk);
  
  Ck.fill(0);
  
  arma::imat Cv(Nk, dtm.n_rows);
  
  Cv.fill(0);
  
  for (auto d = 0; d < Zd.size(); ++d) {
    arma::ivec zd  = Zd[d];
    arma::ivec doc = docs[d];
    
    if (freeze_topics) {
      for (auto n = 0; n < zd.n_elem; ++n) {
        Cd_out(zd[n], d)++;
        Ck[zd[n]]++;
      }
      
    } else {
      for (auto n = 0; n < zd.n_elem; ++n) {
        Cd_out(zd[n], d)++;
        Ck[zd[n]]++;
        Cv(zd[n], doc[n])++;
      }
    }
  }
  
  // ***************************************************************************
  // Prepare output and expel it from this function
  // ***************************************************************************
  using Rcpp::_;
  return Rcpp::List::create(          //
    _["docs"] = Rcpp::wrap(docs),   //
    _["Zd"]   = Rcpp::wrap(Zd),     //
    _["Cd"]   = Rcpp::wrap(Cd_out), //
    _["Cv"]   = Rcpp::wrap(Cv),     //
    _["Ck"]   = Ck                  //
  );                                  //
}

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

//' Main C++ Gibbs sampler for Latent Dirichlet Allocation
//' @keywords internal
//' @description
//'   This is the C++ Gibbs sampler for LDA. "Abandon all hope, ye who enter here."
//' @param Docs List with one element for each document and one entry for each token
//'   as formatted by \code{\link[tidylda]{initialize_topic_counts}}
//' @param Zd_in List with one element for each document and one entry for each token
//'   as formatted by \code{\link[tidylda]{initialize_topic_counts}}
//' @param Cd_in IntegerMatrix denoting counts of topics in documents
//' @param Cv_in IntegerMatrix denoting counts of tokens in topics
//' @param Ck_in IntegerVector denoting counts of topics across all tokens
//' @param beta_in NumericMatrix for prior of tokens over topics
//' @param alpha_in NumericVector prior for topics over documents
//' @param iterations int number of gibbs iterations to run in total
//' @param burnin int number of burn in iterations
//' @param calc_likelihood bool do you want to calculate the log likelihood each
//'   iteration?
//' @param Phi_in NumericMatrix denoting probability of tokens in topics
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
//' @param optimize_alpha bool do you want to optimize alpha each iteration?
//' @param threads unsigned integer, how many parallel threads?
//' @param verbose bool do you want to print out a progress bar?
//' @details
//'   Arguments ending in \code{_in} are copied and their copies modified in
//'   some way by this function. In the case of \code{beta} and \code{Phi},
//'   the only modification is that they are converted from matrices to nested
//'   \code{std::vector} for speed, reliability, and thread safety. In the case
//'   of all others, they may be explicitly modified during training. 
// [[Rcpp::export]]
Rcpp::List fit_lda_c(
    const std::vector<std::vector<std::size_t>>&  Docs, 
    const std::vector<std::vector<std::size_t>>&  Zd_in,
    const IntegerMatrix&                          Cd_in,
    const IntegerMatrix&                          Cv_in,
    const std::vector<long>&                      Ck_in,
    const std::vector<double>                     alpha_in, 
    const NumericMatrix&                          beta_in, 
    const std::size_t&                            iterations, 
    const int&                                    burnin,
    const bool&                                   optimize_alpha, 
    const bool&                                   calc_likelihood,
    const NumericMatrix&                          Phi_in,
    const bool&                                   freeze_topics,
    const std::size_t&                            threads = 1,
    const bool&                                   verbose = false
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
  std::vector<std::vector<long>> Cd_sum(Nd);
  std::vector<std::vector<double>> Cd_mean(Nd);
  
  for (auto d = 0; d < Nd; d++) {
    for (auto k = 0; k < Nk; k++) {
      Cd_sum[d].push_back(0);
      Cd_mean[d].push_back(0.0);
    }
  }
  
  std::vector<std::vector<long>> Cv_sum(Nk);
  std::vector<std::vector<double>> Cv_mean(Nk);
  
  for (auto k = 0; k < Nk; k++) {
    for (auto v = 0; v < Nv; v++) {
      Cv_sum[k].push_back(0);
      Cv_mean[k].push_back(0.0);
    }
  }
  
  // ***********************************************************************
  // Set up variables for parallel sampling. 
  // ***********************************************************************
  
  // containers for thread specific Cv, Ck, and Phi
  // note: Phi does not get updated every iteration as its use means 
  // freeze_topics = true
  std::vector<std::vector<std::vector<long>>> Cv_batch(threads);
  std::vector<std::vector<long>> Ck_batch(threads);
  
  std::vector<std::vector<std::vector<double>>> Phi_batch(threads);
  for (auto j = 0; j < threads; j++) {
    Phi_batch[j] = Phi;
  }
  
  // container of document indices
  std::vector<std::vector<std::size_t>> batch_indices;
  batch_indices = allocate_batch_indices(threads, Nd);
  
  // ***********************************************************************
  // Variables and other set up for the log likelihood function
  // ***********************************************************************
  std::vector<std::vector<double>> log_likelihood(iterations);
  
  for (auto t = 0; t < iterations; t++) { // fill in empty likelihood
    std::vector<double> tmp(3);
    log_likelihood[t] = tmp;
  }
  
  double lgbeta(0.0); // if calc_likelihood, we need this term
  
  double lgalpha(0.0); // if calc_likelihood, we need this term
  
  if (calc_likelihood) { // if calc_likelihood, actuAllocationy populate this stuff
    
    for (auto n = 0; n < Nv; n++) {
      lgbeta += lgamma(beta[0][n]);
    }
    
    lgbeta = (lgbeta - lgamma(sum_beta)) * Nk; 
    
    for (auto k = 0; k < Nk; k++) {
      lgalpha += lgamma(alpha[k]);
    }
    
    lgalpha = (lgalpha - lgamma(sum_alpha)) * Nd;
  }
  
  
  // ***********************************************************************
  // Main Gibbs iterations
  // ***********************************************************************
  for (auto t = 0; t < iterations; t++) { // for each iteration
    
    // distribute copies of Ck and Cv to be local to each thread
    for (auto j = 0; j < threads; j++) {
      Cv_batch[j] = Cv;
      Ck_batch[j] = Ck;
    }
    
    RcppThread::parallelFor(
      0,
      threads,
      [&] (std::size_t j = 0) { // for each thread in parallel
      
      auto batch_idx = batch_indices[j];
      
      // do a loop over all documents in the thread with local Ck and Cv
      for (auto d = batch_idx[0]; d < batch_idx.size(); d++) { 
        
        RcppThread::checkUserInterrupt();
        
        // R_CheckUserInterrupt();
        
        auto doc = Docs[d]; // placeholder for a document
        auto zd = Zd[d]; // placeholder for doc-word-topic assigment
        
        for (auto n = 0; n < doc.size(); n++) { // for each word in that document
          
          // handle things differently based on freeze_topics
          // note some copied code, but minimizes number of checks in this loop
          if (! freeze_topics) {
            
            // discount for the n-th word with topic z
            Cd[d][zd[n]]--; 
            Cv_batch[j][zd[n]][doc[n]]--;
            Ck_batch[j][zd[n]]--;
            
            // calculate probability vector
            for (auto k = 0; k < Nk; k++) {
              qz[k] = (Cv_batch[j][k][doc[n]] + beta[k][doc[n]]) / 
                (Ck_batch[j][k] + sum_beta) *
                (Cd[d][k] + alpha[k]) / 
                (doc.size() + sum_alpha);
            }
            
            // sample topic
            z = samp_one(qz);
            
            // update counts based on sampled topic
            zd[n] = z;
            
            Cd[d][zd[n]]++; 
            Cv_batch[j][zd[n]][doc[n]]++; 
            Ck_batch[j][zd[n]]++; 
            
          } else {
            
            // discount for the n-th word with topic z
            Cd[d][zd[n]]--; 
            
            // calculate probability vector
            for (auto k = 0; k < Nk; k++) {
              qz[k] = Phi_batch[j][k][doc[n]] *
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
      
    }, threads); // end loop over threads
    
    // update global Ck and Cv using batch versions
    Ck = update_global_Ck(
      Ck,
      Ck_batch,
      threads
    );
    
    Cv = update_global_Cv(
      Cv,
      Cv_batch,
      threads
    );
    
    // if using burnin, update sums
    if (burnin > -1 && t >= burnin) {
      
      for (auto d = 0; d < Nd; d++) {
        for (auto k = 0; k < Nk; k++) {
          Cd_sum[d][k] += Cd[d][k];
        }
      }
      
      if (! freeze_topics) {
        for (auto k = 0; k < Nk; k++) {
          for (auto v = 0; v < Nv; v++) {
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
      
      for (auto k = 0; k < Nk; k++) {
        
        std::vector<double> tmp(Nv);
        
        phi_prob[k] = tmp;
        
        // get the denominator
        for (auto v = 0; v < Nv; v++) {
          denom += Cv[k][v] + beta[k][v];
        }
        
        // get the probability
        for (auto v = 0; v < Nv; v++) {
          phi_prob[k][v] = ((double)Cv[k][v] + beta[k][v]) / denom;
          
          lp_beta += (beta[k][v] - 1) * log(phi_prob[k][v]);
        }
      }
      
      lp_beta += lgbeta;
      
      // for each document, get the log probability of the words
      
      double lp_alpha(0.0); // log probability of alpha prior
      
      double lpd(0.0); // log probability of documents
      
      for (auto d = 0; d < Nd; d++) {
        
        std::vector<double> theta_prob(Nk); // probability of each topic in this document
        
        auto doc = Docs[d];
        
        NumericVector lp(doc.size()); // log probability of each word under the model
        
        double denom(0.0);
        
        for (auto k = 0; k < Nk; k++) {
          denom += (double)Cd[d][k] + alpha[k];
        }
        
        for (auto k = 0; k < Nk; k++) {
          theta_prob[k] = ((double)Cd[d][k] + alpha[k]) / denom;
          
          lp_alpha += (alpha[k] - 1) * log(theta_prob[k]);
        }
        
        for (auto n = 0; n < doc.size(); n++) {
          
          lp[n] = 0.0;
          
          for (auto k = 0; k < Nk; k++) {
            
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
      
      for (auto k = 0; k < Nk; k++) {
        alpha[k] = (static_cast<double>(Ck[k]) / static_cast<double>(sum_tokens)) * sum_alpha;
      }
    } 
    
    // progress bar
    if (verbose) {
      Rcout << "=";
      // every 100th iteration, add a new line
      if ((t + 1) % 100 == 0) {
        Rcout << std::endl;
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
    for (auto d = 0; d < Nd; d++) {
      for (auto k = 0; k < Nk; k++) {
        Cd_mean[d][k] = static_cast<double>(Cd_sum[d][k]) / diff;
      }
    }
    
    if (! freeze_topics) {
      for (auto v = 0; v < Nv; v++) {
        for (auto k = 0; k < Nk; k++) {
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
