// Functions to make a collapsed gibbs sampler for LDA

#include "sample_int.h"
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
    [&Cd, &Phi, &dtm, &alpha, &sum_alpha, &docs, &Zd, &Nk](unsigned int d) {
      arma::vec qz(Nk);
      
      // arma::ivec       topic_index = Rcpp::seq_len(Nk) - 1;
      std::vector<int> topic_index(Nk);
      std::iota(std::begin(topic_index), std::end(topic_index), 0);
      
      // make a temporary vector to hold token indices
      auto nd(0);
      
      for (std::size_t v = 0; v < dtm.n_rows; ++v) {
        nd += dtm(v, d);
      }
      
      arma::ivec doc(nd);
      arma::ivec zd(nd);
      arma::ivec z(1);
      
      // fill in with token indices
      std::size_t j = 0; // index of doc, advances when we have non-zero entries
      
      for (std::size_t v = 0; v < dtm.n_rows; ++v) {
        if (dtm(v, d) > 0) { // if non-zero, add elements to doc
          
          // calculate probability of topics based on initially-sampled Phi and Cd
          for (std::size_t k = 0; k < Nk; ++k) {
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
  
  for (std::size_t d = 0; d < Zd.size(); ++d) {
    arma::ivec zd  = Zd[d];
    arma::ivec doc = docs[d];
    
    if (freeze_topics) {
      for (std::size_t n = 0; n < zd.n_elem; ++n) {
        Cd_out(zd[n], d)++;
        Ck[zd[n]]++;
      }
      
    } else {
      for (std::size_t n = 0; n < zd.n_elem; ++n) {
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
// Declare a bunch of voids to be called inside main sampling function
////////////////////////////////////////////////////////////////////////////////
// Functions down here are called inside of calc_lda_c()

// sample a new topic
void sample_topics(const std::vector<int>&    doc,
                   std::vector<int>&          zd,
                   const int                  d,
                   arma::uvec&                Ck,
                   arma::umat&                Cd,
                   arma::mat&                 Cv,
                   const bool                 freeze_topics,
                   const arma::mat&           Phi,
                   const arma::vec&           alpha,
                   const arma::mat&           beta,
                   const double               sum_alpha,
                   const double               sum_beta) {
  
  // initialize qz before loop and reset at end of each iteration
  std::vector<double> qz(Ck.n_elem, 1.0);
  
  
  // for each token instance in the document
  for (std::size_t n = 0; n < doc.size(); n++) {
    // discount counts from previous run ***
    Cd(zd[n], d)--;
    
    // update probabilities of each topic ***
    auto phi_kv(0.0);
    
    if (freeze_topics) {
      for (std::size_t k = 0; k < qz.size(); ++k) {
        // get the correct term depending on if we freeze topics or not
        // prevent branching inside loop by when `freeze_topics` condition
        phi_kv = Phi(k, doc[n]);
        qz[k]  = log(phi_kv) + log(Cd(k, d) + alpha[k]) - log(doc.size() + sum_alpha - 1);
      }
      
    } else {
      Cv(zd[n], doc[n])--;
      Ck[zd[n]]--;
      
      for (std::size_t k = 0; k < qz.size(); ++k) {
        phi_kv = log(Cv(k, doc[n]) + beta(k, doc[n])) - log(Ck[k] + sum_beta);
        qz[k]  = phi_kv + log(Cd(k, d) + alpha[k]) - log(doc.size() + sum_alpha - 1);
      }
    }
    
    // initialize z 
    arma::uvec z(1);
    
    // sample a topic ***
    z = lsamp_one(qz);
    
    // update counts ***
    Cd(z[0], d)++;
    
    if (!freeze_topics) {
      Cv(z[0], doc[n])++;
      Ck[z[0]]++;
    }
    
    // record this topic for this token/doc combination
    zd[n] = z[0];
    
    std::fill(std::begin(qz), std::end(qz), 1.0); // reset qz before next iteration
  }                                               // end loop over each token in doc
}

// self explanatory: calculates the (log) likelihood
void fcalc_likelihood(const std::size_t t,
                      const double      sum_beta,
                      const arma::uvec& Ck,
                      const arma::umat& Cd,
                      const arma::mat&  Cv,
                      const arma::vec&  alpha,
                      const arma::mat&  beta,
                      const double      lgalpha,
                      const double      lgbeta,
                      const double      lg_alpha_len,
                      arma::mat&        log_likelihood) {
  // calculate lg_beta_count1, lg_beta_count2, lg_alph_count for this iter
  // start by zeroing them out
  auto lg_beta_count1(0.0);
  auto lg_beta_count2(0.0);
  auto lg_alpha_count(0.0);
  
  for (std::size_t k = 0; k < Cd.n_rows; k++) {
    lg_beta_count1 += lgamma(sum_beta + Ck[k]);
  }
  
  for (std::size_t d = 0; d < Cd.n_cols; d++) {
    for (std::size_t k = 0; k < Cd.n_rows; k++) {
      lg_alpha_count += lgamma(alpha[k] + Cd(k, d));
    }
  }
  
  for (std::size_t v = 0; v < Cv.n_cols; v++) {
    for (std::size_t k = 0; k < Cv.n_rows; k++) {
      lg_beta_count2 += lgamma(beta(k, v) + Cv(k, v));
    }
  }
  
  lg_beta_count1 *= -1;
  log_likelihood(0, t) = t;
  log_likelihood(1, t) =
    lgalpha + lgbeta + lg_alpha_len + lg_alpha_count + lg_beta_count1 + lg_beta_count2;
}

// if user wants to optimize alpha, do that here.
// procedure likely to change similar to what Mimno does in Mallet
void foptimize_alpha(arma::vec&        alpha,
                     const arma::uvec& Ck,
                     const std::size_t sumtokens,
                     const double      sum_alpha,
                     const std::size_t Nk) {
  
  constexpr double denom = 2.0;
  for (std::size_t k = 0; k < Nk; ++k) {
    alpha[k] += ((Ck[k] / static_cast<double>(sumtokens) * sum_alpha) + alpha[k]) / denom;
  }
}

// Function aggregates counts across iterations after burnin iterations
void agg_counts_post_burnin(const bool        freeze_topics,
                            const arma::umat& Cd,
                            arma::umat&       Cd_sum,
                            const arma::mat&  Cv,
                            arma::mat&        Cv_sum) {
  
  if (freeze_topics) { // split these up to prevent branching inside loop
    
    for (std::size_t d = 0; d < Cd.n_cols; d++) { // consider parallelization
      for (std::size_t k = 0; k < Cd.n_rows; k++) {
        Cd_sum(k, d) += Cd(k, d);
      }
    }
    
    for (std::size_t v = 0; v < Cv.n_cols; v++) {
      for (std::size_t k = 0; k < Cv.n_rows; k++) {
        Cv_sum(k, v) += Cv(k, v);
      }
    }
    
  } else {
    for (std::size_t d = 0; d < Cd.n_cols; ++d) {
      for (std::size_t k = 0; k < Cd.n_rows; ++k) {
        Cd_sum(k, d) += Cd(k, d);
      }
    }
  }
}

//' Main C++ Gibbs sampler for Latent Dirichlet Allocation
//' @keywords internal
//' @description
//'   This is the C++ Gibbs sampler for LDA. "Abandon all hope, ye who enter here."
//' @param docs List with one element for each document and one entry for each token
//'   as formatted by \code{\link[tidylda]{initialize_topic_counts}}
//' @param Nk int number of topics
//' @param beta NumericMatrix for prior of tokens over topics
//' @param alpha NumericVector prior for topics over documents
//' @param Cd IntegerMatrix denoting counts of topics in documents
//' @param Cv IntegerMatrix denoting counts of tokens in topics
//' @param Ck IntegerVector denoting counts of topics across all tokens
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
Rcpp::List fit_lda_c(const std::vector<std::vector<int>>& docs,
                     const int                            Nk,
                     const arma::mat&                     beta,
                     arma::vec&                           alpha,
                     arma::umat&                          Cd,
                     arma::mat&                           Cv,
                     arma::uvec&                          Ck,
                     std::vector<std::vector<int>>&       Zd,
                     const arma::mat&                     Phi,
                     const int                            iterations,
                     const int                            burnin,
                     const bool                           freeze_topics,
                     const bool                           calc_likelihood,
                     const bool                           optimize_alpha) {
  
  // ***********************************************************************
  // Variables and other set up
  // ***********************************************************************
  
  const auto Nv = Cv.n_cols;
  const auto Nd = Cd.n_cols;
  
  // std::vector<double> k_alpha(alpha.size());
  // std::transform(std::begin(alpha),
  //                std::end(alpha),
  //                std::begin(k_alpha),
  //                [&Nk](double x) { return x * Nk; });
  
  arma::vec k_alpha = alpha * Nk;
  
  const arma::mat v_beta = beta * Nv;
  const auto sum_alpha = arma::sum(alpha);
  const auto sum_beta  = arma::sum(beta.row(1)); // needs updating for topic seeding
  const auto sumtokens = arma::sum(Ck);
  
  // related to burnin and averaging
  arma::mat Cv_sum(Nk, Nv);
  arma::mat Cv_mean(Nk, Nv);
  arma::umat Cd_sum(Nk, Nd);
  arma::mat Cd_mean(Nk, Nd);
  
  // related to the likelihood calculation
  arma::mat log_likelihood(2, iterations);
  
  auto lgbeta(0.0);       // calculated immediately below
  auto lgalpha(0.0);      // calculated immediately below
  auto lg_alpha_len(0.0); // calculated immediately below
  
  // indices for loops
  auto t(0);
  auto d(0);
  auto n(0);
  auto k(0);
  auto v(0);
  
  if (calc_likelihood &&
      !freeze_topics) { // if calc_likelihood, actually populate this stuff
      
      for (; n < Nv; ++n) {
        lgbeta += lgamma(beta[n]);
      }
      
      lgbeta = (lgbeta - lgamma(sum_beta)) * Nk; // rcpp sugar here
    
    for (; k < Nk; ++k) {
      lgalpha += lgamma(alpha[k]);
    }
    
    lgalpha = (lgalpha - lgamma(sum_alpha)) * Nd;
    
    for (d = 0; d < Nd; ++d) {
      lg_alpha_len += lgamma(sum_alpha + docs[d].size());
    }
    
    lg_alpha_len *= -1;
  }
  
  // ***********************************************************************
  // BEGIN ITERATIONS
  // ***********************************************************************
  
  for (; t < iterations; ++t) {
    // loop over documents
    for (; d < Nd; ++d) { // start loop over documents
      
      R_CheckUserInterrupt();
      
      auto doc = docs[d];
      auto zd = Zd[d];
      
      sample_topics(doc,
                    zd,
                    d,
                    Ck,
                    Cd,
                    Cv,
                    freeze_topics,
                    Phi,
                    alpha,
                    beta,
                    sum_alpha,
                    sum_beta);
      
    } // end loop over docs
    // calc likelihood ***
    if (calc_likelihood && !freeze_topics) {
      fcalc_likelihood(t,
                       sum_beta,
                       Ck,
                       Cd,
                       Cv,
                       alpha,
                       beta,
                       lgalpha,
                       lgbeta,
                       lg_alpha_len,
                       log_likelihood);
    }
    // optimize alpha ***
    if (optimize_alpha && !freeze_topics) {
      foptimize_alpha(alpha, Ck, sumtokens, sum_alpha, Nk);
    }
    
    // aggregate counts after burnin ***
    if (burnin > -1 && t >= burnin) {
      agg_counts_post_burnin(freeze_topics, Cd, Cd_sum, Cv, Cv_sum);
    }
    
  } // end iterations
  
  // ***********************************************************************
  // Cleanup and return list
  // ***********************************************************************
  
  // change sum over iterations to average over iterations ***
  if (burnin > -1) {
    const double diff(iterations - burnin);
    
    // average over chain after burnin
    for (; k < Nk; ++k) {
      for (; d < Nd; ++d) {
        Cd_mean(k, d) = Cd_sum(k, d) / diff;
      }
      
      for (; v < Nv; ++v) {
        Cv_mean(k, v) = Cv_sum(k, v) / diff;
      }
    }
  }
  
  // Return the final list ***
  using Rcpp::_;
  return Rcpp::List::create(                //
    _["Cd"]             = Cd,             //
    _["Cv"]             = Cv,             //
    _["Ck"]             = Ck,             //
    _["Cd_mean"]        = Cd_mean,        //
    _["Cv_mean"]        = Cv_mean,        //
    _["Cd_sum"]         = Cd_sum,        //
    _["Cv_sum"]         = Cv_sum,        //
    _["log_likelihood"] = log_likelihood, //
    _["alpha"]          = alpha,          //
    _["beta"]           = beta            //
  );                                        //
}
