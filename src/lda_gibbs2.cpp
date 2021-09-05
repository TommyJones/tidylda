// Functions to make a collapsed gibbs sampler for LDA

// Note that at one time, functions here used parallel for loops to do LDA
// in the spirit of Newman et al. 2009 (ish). That proved problematic for
// several reasons. So I've rolled this back to be sequential only.
// I've left a lot of the infrastructure here for parallel
// sampling (dividing data into batches, for example). That comes with computational
// costs, but allows easier re-parallelization down the road if desired.


#define ARMA_64BIT_WORD 1

#include "sample_int.h" // RcppArmadillo.h included here
#include "matrix_conversions.h"
#include "parallel_gibbs_utils.h"

#include <progress.hpp>
#include <progress_bar.hpp>



//' Make a lexicon for looping over in the gibbs sampler
//' @keywords internal
//' @description
//'   One run of the Gibbs sampler and other magic to initialize some objects.
//'   Works in concert with \code{\link[tidylda]{initialize_topic_counts}}.
//' @param Cd_in IntegerMatrix denoting counts of topics in documents
//' @param Beta_in NumericMatrix denoting probability of words in topics
//' @param dtm_in arma::sp_mat document term matrix
//' @param alpha NumericVector prior for topics over documents
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
//' @return Returns a list with five entries.
//' 
//'   \code{Docs} is a list of vectors. Each element is a document, and the contents
//'   are indices for tokens. Used as an iterator for the Gibbs sampler.
//'   
//'   \code{Zd} is a list of vectors, similar to Docs. However, its contents are topic
//'   assignments of each document/token pair. Used as an iterator for Gibbs
//'   sampling.
//'   
//'   \code{Cd} is a matrix counting the number of times each topic is sampled per
//'   document.
//'   
//'   \code{Cv} is a matrix counting the number of times each topic is sampled per token.
//'   
//'   \code{Ck} is a vector counting the total number of times each topic is sampled overall.
//'   
//'   \code{Cd}, \code{Cv}, and \code{Ck} are derivatives of \code{Zd}.
//' @details
//'   Arguments ending in \code{_in} are copied and their copies modified in
//'   some way by this function. In the case of \code{Cd_in} and \code{Beta_in},
//'   the only modification is that they are converted from matrices to nested
//'   \code{std::vector} for speed, reliability, and thread safety. \code{dtm_in}
//'   is transposed for speed when looping over columns. 
// [[Rcpp::export]]
Rcpp::List create_lexicon(
    const IntegerMatrix&       Cd_in,
    const NumericMatrix&       Beta_in,
    const arma::sp_mat&        dtm_in,
    const std::vector<double>& alpha,
    const bool&                freeze_topics
) {
  // ***************************************************************************
  // Initialize some variables
  // ***************************************************************************
  
  arma::sp_mat dtm = dtm_in.t(); // transpose dtm to take advantage of column major & parallel
  // Cd  = Cd.t();  // transpose to put columns up
  
  auto Cd = mat_to_vec(Cd_in, true);
  auto Beta = mat_to_vec(Beta_in, true);
  
  std::vector<std::vector<std::size_t>> Docs(dtm.n_cols);
  std::vector<std::vector<std::size_t>> Zd(dtm.n_cols);
  
  const std::size_t Nv = dtm.n_rows;
  const std::size_t Nd = Cd.size();
  const std::size_t Nk = Beta.size();
  
  auto sum_alpha = std::accumulate(alpha.begin(), alpha.end(), 0.0); 
  
  // ***************************************************************************
  // Go through each document and split it into a lexicon and then sample a
  // topic for each token within that document
  // ***************************************************************************
  for (auto d = 0; d < Nd; d++) {
    std::vector<double> qz(Nk); // placehodler for probability of topic
    
    arma::uword z; // placeholder for topic sampled
    
    // make a temporary vector to hold token indices
    auto nd(0);
    
    for (auto v = 0; v < Nv; v++) {
      nd += dtm(v, d);
    }
    
    std::vector<std::size_t> doc(nd);
    std::vector<std::size_t> zd(nd);
    
    // fill in with token indices
    std::size_t j = 0; // index of doc, advances when we have non-zero entries
    
    for (auto v = 0; v < Nv; v++) {
      if (dtm(v, d) > 0) { // if non-zero, add elements to doc
        
        // calculate probability of topics based on initially-sampled Beta and Cd
        for (auto k = 0; k < Nk; k++) {
          qz[k] = log(Beta[k][v]) + log(Cd[d][k] + alpha[k]) - log(nd + sum_alpha - 1);
        }
        
        std::size_t idx(j + dtm(v, d)); // where to stop the loop below
        
        while (j < idx) {
          doc[j] = v;
          z      = lsamp_one(qz); // sample a topic here
          zd[j]  = z;
          j++;
        }
      }
    }
    
    // fill in docs[d] with the matrix we made
    Docs[d] = doc;
    
    Zd[d] = zd;
    
    Rcpp::checkUserInterrupt();    
  }
  
  
  // ***************************************************************************
  // Calculate Cd, Cv, and Ck from the sampled topics
  // ***************************************************************************
  
  // initialize Ck
  std::vector<long> Ck(Nk);
  
  for (auto k = 0; k < Nk; k++) {
    Ck[k] = 0;
  }
  
  // initialize Cv; Unlike in fit_lda_c(), this Cv is token major, not topic major
  // this is so it's accessed efficiently below
  std::vector<std::vector<long>> Cv(Nv);
  
  for (auto v = 0; v < Nv; v++) {
    std::vector<long> cv_row(Nk);
    for (auto k = 0; k < Nk; k++) {
      cv_row[k] = 0;
    }
    Cv[v] = cv_row;
    
    Rcpp::checkUserInterrupt();
    
  }
  
  // loop over Zd and Docs to count topics sampled for each token
  // can't do this loop in parallel because of potential collisions in
  // Ck and Cv
  for (auto d = 0; d < Nd; d++) {
    
    std::vector<std::size_t> zd = Zd[d];
    std::vector<std::size_t> doc = Docs[d];
    
    // initialize vector of zeros
    std::vector<long> cd_row(Nk);
    
    for (auto k = 0; k < Nk; k++) {
      cd_row[k] = 0;
    }
    
    // count Cd, Cv, Ck based on whether or not we freeze topics
    if (freeze_topics) { // only count Cd
      
      for (auto n = 0; n < zd.size(); n++) {
        cd_row[zd[n]]++;
      }
      
      Cd[d] = cd_row;
      
    } else { // count Cd, Cv, Ck
      
      for (auto n = 0; n < zd.size(); n++) {
        cd_row[zd[n]]++;
        Cv[doc[n]][zd[n]]++;
        Ck[zd[n]]++;
      }
      
      Cd[d] = cd_row;
    }
    
    Rcpp::checkUserInterrupt();
  }
  
  // ***************************************************************************
  // Prepare output and expel it from this function
  // ***************************************************************************
  using Rcpp::_;
  return Rcpp::List::create(
    _["Docs"] = Docs,
    _["Zd"]   = Zd,
    _["Cd"]   = vec_to_mat(Cd, true),
    _["Cv"]   = vec_to_mat(Cv),
    _["Ck"]   = Ck
  );
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
//' @param eta_in NumericMatrix for prior of tokens over topics
//' @param alpha_in NumericVector prior for topics over documents
//' @param iterations int number of gibbs iterations to run in total
//' @param burnin int number of burn in iterations
//' @param calc_likelihood bool do you want to calculate the log likelihood each
//'   iteration?
//' @param Beta_in NumericMatrix denoting probability of tokens in topics
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
//' @param optimize_alpha bool do you want to optimize alpha each iteration?
//' @param threads unsigned integer, how many parallel threads?
//'        For now, nothing is actually parallel
//' @param verbose bool do you want to print out a progress bar?
//' @return Returns a list with the following entries.
//' 
//'   \code{Cd} is a matrix counting the number of times each topic is sampled per
//'   document.
//'   
//'   \code{Cv} is a matrix counting the number of times each topic is sampled per token.
//'   
//'   \code{Cd_mean} the same as \code{Cd} but values averaged across iterations
//'   greater than \code{burnin} iterations.
//'   
//'   \code{Cv_mean} the same as \code{Cv} but values averaged across iterations
//'   greater than \code{burnin} iterations.
//'   
//'   \code{Cd_sum} the same as \code{Cd} but values summed across iterations
//'   greater than \code{burnin} iterations.
//'   
//'   \code{Cv_sum} the same as \code{Cv} but values summed across iterations
//'   greater than \code{burnin} iterations.
//'   
//'   \code{log_likelihood} a matrix with one row indexing iterations and one
//'   row of the log likelihood for each iteration.
//'   
//'   \code{alpha} a vector of the document-topic prior
//'   
//'   \code{_eta} a matrix of the topic-token prior
//' @details
//'   Arguments ending in \code{_in} are copied and their copies modified in
//'   some way by this function. In the case of \code{eta_in} and \code{Beta_in},
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
    const NumericMatrix&                          eta_in, 
    const std::size_t&                            iterations, 
    const int&                                    burnin,
    const bool&                                   optimize_alpha, 
    const bool&                                   calc_likelihood,
    const NumericMatrix&                          Beta_in,
    const bool&                                   freeze_topics,
    const std::size_t&                            threads = 1,
    const bool&                                   verbose = true
) {
  
  // ***********************************************************************
  // Set up progress bar
  // ***********************************************************************
  Progress p(iterations, verbose);
  
  // ***********************************************************************
  // Convert input matrices to 2-dimensional std::vector for speed
  // Also convert some other objects to avoid overwriting inputs
  // ***********************************************************************
  auto Zd = Zd_in;
  
  auto Cd = mat_to_vec(Cd_in, true);
  auto Cv = mat_to_vec(Cv_in, true);
  auto Ck = Ck_in;
  
  auto alpha = alpha_in;
  auto eta = mat_to_vec(eta_in, true);
  
  auto Beta = mat_to_vec(Beta_in, true);
  
  // ***********************************************************************
  // Variables and other set up for the main sampler
  // ***********************************************************************
  
  const std::size_t Nv = Cv[0].size();
  const std::size_t Nd = Cd.size();
  const std::size_t Nk = Cv.size();
  
  auto sum_tokens = std::accumulate(Ck.begin(), Ck.end(), 0);
  auto sum_alpha = std::accumulate(alpha.begin(), alpha.end(), 0.0); 
  
  std::vector<double> sum_eta(eta.size());
  
  for (auto k = 0; k < Nk; k++) {
    sum_eta[k] = std::accumulate(eta[k].begin(), eta[k].end(), 0.0);
  }

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
  
  // containers for thread specific Cv, Ck, and Beta
  // note: Beta does not get updated every iteration as its use means 
  // freeze_topics = true
  std::vector<std::vector<std::vector<long>>> Cv_batch(threads);
  std::vector<std::vector<long>> Ck_batch(threads);
  
  std::vector<std::vector<std::vector<double>>> Beta_batch(threads);
  for (auto j = 0; j < threads; j++) {
    Beta_batch[j] = Beta;
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
  
  double lgeta(0.0); // if calc_likelihood, we need this term
  
  double lgalpha(0.0); // if calc_likelihood, we need this term
  
  if (calc_likelihood) { // if calc_likelihood, actually populate this stuff
    
    for (auto v = 0; v < Nv; v++) {
      lgeta += lgamma(eta[0][v]);
    }
    
    lgeta = (lgeta - lgamma(sum_eta[0])) * Nk; 
    
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
    
    for (auto j = 0; j < threads; j++) { // for each thread (used to be parallel)
      auto batch_idx = batch_indices[j];
      
      std::vector<double> qz(Nk); // placehodler for probability of topic
      
      arma::uword z; // placeholder for topic sampled
      
      
      // do a loop over all documents in the thread with local Ck and Cv
      for (auto d = batch_idx[0]; d < batch_idx.size(); d++) { 
        
        Rcpp::checkUserInterrupt();    
        
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
              qz[k] = (Cv_batch[j][k][doc[n]] + eta[k][doc[n]]) / 
                (Ck_batch[j][k] + sum_eta[k]) *
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
              qz[k] = Beta_batch[j][k][doc[n]] *
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
    } // end loop over threads
    

    // update global Ck and Cv using batch versions
    if (threads > 1) {
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
    } else {
      Ck = Ck_batch[0];
      
      Cv = Cv_batch[0];
    }
    
    
    
    // shuffle batch indices for next run
    // if (threads > 1) {
    //   shuffle_batch_indices(
    //     batch_indices,
    //     Nd
    //   );
    // }
    
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
      }// end if
      
    } // end if
    
    // if calculating log likelihood, do so 
    if (calc_likelihood && !freeze_topics) {
      
      std::vector<double> tmp(3);
      
      // get beta probability matrix @ this iteration
      std::vector<std::vector<double>> beta_prob(Nk);
      
      double denom(0.0);
      
      double lp_eta(0.0); // log probability of eta prior
      
      for (auto k = 0; k < Nk; k++) {
        
        std::vector<double> tmp(Nv);
        
        beta_prob[k] = tmp;
        
        // get the denominator
        for (auto v = 0; v < Nv; v++) {
          denom += Cv[k][v] + eta[k][v];
        }
        
        // get the probability
        for (auto v = 0; v < Nv; v++) {
          beta_prob[k][v] = ((double)Cv[k][v] + eta[k][v]) / denom;
          
          lp_eta += (eta[k][v] - 1) * log(beta_prob[k][v]);
        }
      }
      
      lp_eta += lgeta;
      
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
            
            lp[n] += theta_prob[k] * beta_prob[k][doc[n]];
            
          }
          
          lpd += log(lp[n]);
        }
      }
      
      lp_alpha += lgalpha;
      
      
      tmp[0] = static_cast<double>(t);
      
      tmp[1] = lpd; // log probability of whole corpus under the model w/o priors
      
      tmp[2] = lpd + lp_alpha + lp_eta; // second likelihood calculation here
      
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
      p.increment(); // update progress
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
      for (auto k = 0; k < Nk; k++) {
        for (auto v = 0; v < Nv; v++) {
          Cv_mean[k][v] = static_cast<double>(Cv_sum[k][v]) / diff;
        }
      }
    }
    
  } // end if
  
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
    _["eta"]           = eta_in // does not get modified, so don't waste compute converting eta
  );
}
