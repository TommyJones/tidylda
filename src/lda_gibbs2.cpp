// Functions to make a collapsed gibbs sampler for LDA

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#define ARMA_64BIT_WORD

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "RcppThread.h"

#include <R.h>

#include <cmath>

#include <Rcpp.h>
using namespace Rcpp;

#include<vector>
using namespace std;


//' Make a lexicon for looping over in the gibbs sampler
//' @keywords internal
//' @description
//'   One run of the Gibbs sampler and other magic to initialize some objects.
//'   Works in concert with \code{\link[tidylda]{initialize_topic_counts}}.
//' @param Cd arma::imat denoting counts of topics in documents
//' @param Phi arma::mat denoting probability of words in topics
//' @param dtm arma::sp_mat document term matrix
//' @param alpha arma::vec prior for topics over documents
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
//[[Rcpp::export]]
List create_lexicon(
    arma::imat &Cd, 
    arma::mat &Phi, 
    arma::sp_mat &dtm,
    arma::vec alpha,
    bool freeze_topics,
    int threads
) {
  
  // ***************************************************************************
  // Initialize some variables
  // ***************************************************************************
  
  dtm = dtm.t(); // transpose dtm to take advantage of column major & parallel
  
  Cd = Cd.t(); // transpose to put columns up
  
  double sum_alpha = sum(alpha);
  
  std::vector<arma::ivec> docs(dtm.n_cols); 
  
  std::vector<arma::ivec> Zd(dtm.n_cols);
  
  unsigned int Nk = Cd.n_rows;
  
  // ***************************************************************************
  // Go through each document and split it into a lexicon and then sample a 
  // topic for each token within that document
  // ***************************************************************************
  RcppThread::parallelFor(
    0, 
    dtm.n_cols, 
    [&Cd,
     &Phi,
     &dtm,
     &alpha,
     &sum_alpha,
     &docs,
     &Zd,
     &Nk
    ] (unsigned int d) {
      
      arma::vec qz(Nk);
      
      arma::ivec topic_index = seq_len(Nk) - 1;
      
      // make a temporary vector to hold token indices
      int nd = 0;
      
      for (int v = 0; v < dtm.n_rows; v++) {
        nd += dtm(v, d);
      }
      
      arma::ivec doc(nd);
      
      arma::ivec zd(nd);
      
      arma::ivec z(1);
      
      // fill in with token indices
      unsigned int j = 0; // index of doc, advances when we have non-zero entries 
      
      for (unsigned int v = 0; v < dtm.n_rows; v++) {
        
        if (dtm(v, d) > 0) { // if non-zero, add elements to doc
          
          // calculate probability of topics based on initially-sampled Phi and Cd
          for (unsigned int k = 0; k < Nk; k++) {
            qz[k] = Phi(k, v) * ((double)Cd(k, d) + alpha[k]) / ((double)nd + sum_alpha - 1);
          }
          
          int idx = j + dtm(v, d); // where to stop the loop below
          
          while (j < idx) {
            
            doc[j] = v;
            
            z = RcppArmadillo::sample(topic_index, 1, false, qz);
            
            zd[j] = z[0];
            
            j += 1;
          }
          
        }
      }
      
      // fill in docs[d] with the matrix we made
      docs[d] = doc;
      
      Zd[d] = zd;
      
      RcppThread::checkUserInterrupt();
      
    }, 
    threads
  );
  
  // ***************************************************************************
  // Calculate Cd, Cv, and Ck from the sampled topics
  // ***************************************************************************
  arma::imat Cd_out(Nk, dtm.n_cols);
  
  Cd_out.fill(0);
  
  arma::ivec Ck(Nk);
  
  Ck.fill(0);
  
  arma::imat Cv(Nk, dtm.n_rows);
  
  Cv.fill(0);
  
  for (unsigned int d = 0; d < Zd.size(); d++) {
    arma::ivec zd = Zd[d]; 
    
    arma::ivec doc = docs[d];
    
    for (unsigned int n = 0; n < zd.n_elem; n++) {
      
      Cd_out(zd[n], d) += 1;
      
      Ck[zd[n]] += 1;
      
      if (! freeze_topics) {
        Cv(zd[n], doc[n]) += 1;
      }
      
    } 
    
  }
  
  
  // ***************************************************************************
  // Prepare output and expel it from this function
  // ***************************************************************************
  
  return List::create(
    Named("docs") = wrap(docs),
    Named("Zd") = wrap(Zd),
    Named("Cd") = wrap(Cd_out),
    Named("Cv") = wrap(Cv),
    Named("Ck") = wrap(Ck)
  );
  
}


////////////////////////////////////////////////////////////////////////////////
// Declare a bunch of voids to be called inside main sampling function
////////////////////////////////////////////////////////////////////////////////
// Functions down here are called inside of calc_lda_c()

// sample a new topic
void sample_topics(
    unsigned int &d,
    arma::ivec& doc,
    IntegerVector& zd,
    arma::ivec& Ck,
    arma::imat& Cd, 
    arma::mat& Cv,
    arma::ivec& topic_index,
    bool& freeze_topics,
    arma::mat& Phi,
    arma::vec& alpha,
    arma::mat& beta,
    double& sum_alpha,
    double& sum_beta,
    double& phi_kv
) {
  // initialize some variables
  arma::vec qz(topic_index.n_elem);
  
  qz.fill(1.0);
  
  arma::ivec z(1);
  
  // for each token instance in the document
  for (unsigned int n = 0; n < doc.n_elem; n++) {
    
    // discount counts from previous run ***
    Cd(zd[n], d) -= 1; 
    
    
    if (! freeze_topics) {
      Cv(zd[n], doc[n]) -= 1; 
      
      Ck[zd[n]] -= 1;
    }
    
    
    // update probabilities of each topic ***
    for (unsigned int k = 0; k < qz.n_elem; k++) {
      
      // get the correct term depending on if we freeze topics or not
      if (freeze_topics) {
        phi_kv = Phi(k, doc[n]);
      } else {
        phi_kv = ((double)Cv(k, doc[n]) + beta(k, doc[n])) /
          ((double)Ck[k] + sum_beta);
      }
      
      qz[k] =  phi_kv * ((double)Cd(k, d) + alpha[k]) / 
        ((double)doc.n_elem + sum_alpha - 1);
      
    }
    
    // sample a topic ***
    z = RcppArmadillo::sample(topic_index, 1, false, qz);
    
    // update counts ***
    Cd(z[0], d) += 1; 
    
    if (! freeze_topics) {
      
      Cv(z[0], doc[n]) += 1; 
      
      Ck[z[0]] += 1;
      
    }
    
    // record this topic for this token/doc combination
    zd[n] = z[0];
    
  } // end loop over each token in doc
  
}

// self explanatory: calculates the (log) likelihood
void fcalc_likelihood(
    double& lg_beta_count1,
    double& lg_beta_count2,
    double& lg_alpha_count,
    unsigned int& t,
    double& sum_beta,
    arma::ivec& Ck,
    arma::imat& Cd,
    arma::mat& Cv,
    arma::vec& alpha,
    arma::mat& beta,
    double& lgalpha,
    double& lgbeta,
    double& lg_alpha_len,
    arma::mat& log_likelihood
) {
  
  // calculate lg_beta_count1, lg_beta_count2, lg_alph_count for this iter
  // start by zeroing them out
  lg_beta_count1 = 0.0;
  lg_beta_count2 = 0.0;
  lg_alpha_count = 0.0;
  
  for (unsigned int k = 0; k < Ck.n_elem; k++) {
    
    lg_beta_count1 += lgamma(sum_beta + Ck[k]);
    
    for (unsigned int d = 0; d < Cd.n_cols; d++) {
      lg_alpha_count += lgamma(alpha[k] + Cd(k, d));
    }
    
    for (unsigned int v = 0; v < Cv.n_cols; v++) {
      lg_beta_count2 += lgamma(beta(k,v) + Cv(k, v));
    }
    
  }
  
  lg_beta_count1 *= -1;
  
  log_likelihood(0, t) = t;
  
  log_likelihood(1, t) = lgalpha + lgbeta + lg_alpha_len + lg_alpha_count + 
    lg_beta_count1 + lg_beta_count2;
  
}

// if user wants to optimize alpha, do that here.
// procedure likely to change similar to what Mimno does in Mallet
void foptimize_alpha(
    arma::vec& alpha, 
    arma::ivec& Ck,
    unsigned int& sumtokens,
    double& sum_alpha
) {
  
  arma::vec new_alpha(alpha.n_elem);
  
  new_alpha.fill(0.0);
  
  for (unsigned int k = 0; k < alpha.n_elem; k++) {
    
    new_alpha[k] += (double)Ck[k] / (double)sumtokens * (double)sum_alpha;
    
    new_alpha[k] += (new_alpha[k] + alpha[k]) / 2;
    
  }
  
  alpha = new_alpha;
  
}

// Function aggregates counts across iterations after burnin iterations
void agg_counts_post_burnin(
    bool& freeze_topics,
    arma::imat& Cd,
    arma::imat& Cd_sum,
    arma::mat& Cv,
    arma::mat& Cv_sum
) {
  
  for (unsigned int d = 0; d < Cd.n_cols; d++) { // consider parallelization here
    for (unsigned int k = 0; k < Cd.n_rows; k++) {
      
      Cd_sum(k, d) += Cd(k, d);
      
    }
  }
  
  if (! freeze_topics) {
    for (unsigned int v = 0; v < Cv.n_cols; v++) { // consider parallelization
      for (unsigned int k = 0; k < Cv.n_rows; k++) { 
        
        Cv_sum(k, v) += Cv(k, v);
        
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
//' @param beta arma::mat for prior of tokens over topics
//' @param alpha arma::vec prior for topics over documents
//' @param Cd arma::imat denoting counts of topics in documents
//' @param Cv arma::imat denoting counts of tokens in topics
//' @param Ck arma::ivec denoting counts of topics across all tokens
//' @param Zd List with one element for each document and one entry for each token
//'   as formatted by \code{\link[tidylda]{initialize_topic_counts}}
//' @param Phi arma::mat denoting probability of tokens in topics
//' @param iterations int number of gibbs iterations to run in total
//' @param burnin int number of burn in iterations
//' @param freeze_topics bool if making predictions, set to \code{TRUE}
//' @param calc_likelihood bool do you want to calculate the log likelihood each iteration?
//' @param optimize_alpha bool do you want to optimize alpha each iteration?
// [[Rcpp::export]]
List fit_lda_c(
    std::vector<arma::ivec> &docs,
    unsigned int &Nk,
    arma::mat &beta,
    arma::vec alpha,
    arma::imat Cd,
    arma::mat Cv,
    arma::ivec Ck,
    std::vector<IntegerVector> Zd,
    arma::mat &Phi,
    int &iterations,
    int &burnin,
    bool &freeze_topics,
    bool &calc_likelihood,
    bool &optimize_alpha
) {
  
  // ***********************************************************************
  // Variables and other set up
  // ***********************************************************************
  
  // set up some global variables
  unsigned int Nv = Cv.n_cols;
  
  unsigned int Nd = Cd.n_cols;
  
  arma::vec k_alpha = alpha * Nk;
  
  arma::mat v_beta = beta * Nv;
  
  double sum_alpha = arma::sum(alpha);
  
  double sum_beta = arma::sum(beta.col(0));
  
  unsigned int sumtokens = sum(Ck);
  
  double phi_kv(0.0);
  
  arma::ivec topic_index = seq_len(Nk) - 1;
  
  // variables for averaging post burn in
  arma::mat Cv_sum(Nk, Nv);
  
  arma::mat Cv_mean(Nk, Nv);
  
  arma::imat Cd_sum(Nk, Nd);
  
  arma::mat Cd_mean(Nk, Nd);
  
  // variables related to the likelihood calculation
  arma::mat log_likelihood(2, iterations);
  
  double lgbeta(0.0); // calculated immediately below
  
  double lgalpha(0.0); // calculated immediately below
  
  double lg_alpha_len(0.0); // calculated immediately below
  
  double lg_beta_count1(0.0); // calculated at the bottom of the iteration loop
  
  double lg_beta_count2(0.0); // calculated at the bottom of the iteration loop
  
  double lg_alpha_count(0.0); // calculated at the bottom of the iteration loop
  
  if (calc_likelihood && ! freeze_topics) { // if calc_likelihood, actually populate this stuff
    
    for (unsigned int n = 0; n < Nv; n++) {
      lgbeta += lgamma(beta[n]);
    }
    
    lgbeta = (lgbeta - lgamma(sum_beta)) * Nk; // rcpp sugar here
    
    for (unsigned int k = 0; k < Nk; k++) {
      lgalpha += lgamma(alpha[k]);
    }
    
    lgalpha = (lgalpha - lgamma(sum_alpha)) * Nd;
    
    for (unsigned int d = 0; d < Nd; d++) {
      arma::ivec doc = docs[d];
      
      lg_alpha_len += lgamma(sum_alpha + doc.n_elem);
    }
    
    lg_alpha_len *= -1;
  }
  
  
  
  // ***********************************************************************
  // BEGIN ITERATIONS
  // ***********************************************************************
  
  for (unsigned int t = 0; t < iterations; t++) {
    
    // loop over documents
    for (unsigned int d = 0; d < Nd; d++) { //start loop over documents
      
      R_CheckUserInterrupt();
      
      arma::ivec doc = docs[d];
      
      IntegerVector zd = Zd[d];
      
      sample_topics(
        d,
        doc,
        zd,
        Ck,
        Cd, 
        Cv,
        topic_index,
        freeze_topics,
        Phi,
        alpha,
        beta,
        sum_alpha,
        sum_beta,
        phi_kv
      );
      
    } // end loop over docs    
    // calc likelihood ***
    if (calc_likelihood && ! freeze_topics) {
      fcalc_likelihood(
        lg_beta_count1,
        lg_beta_count2,
        lg_alpha_count,
        t,
        sum_beta,
        Ck,
        Cd,
        Cv,
        alpha,
        beta,
        lgalpha,
        lgbeta,
        lg_alpha_len,
        log_likelihood
      );
    }
    // optimize alpha ***
    if (optimize_alpha && ! freeze_topics) {
      foptimize_alpha(
        alpha, 
        Ck,
        sumtokens,
        sum_alpha
      );  
    }
    
    // aggregate counts after burnin ***
    if (burnin > -1 && t >= burnin) {
      agg_counts_post_burnin(
        freeze_topics,
        Cd,
        Cd_sum,
        Cv,
        Cv_sum
      );
    }
    
  } // end iterations
  
  // ***********************************************************************
  // Cleanup and return list
  // ***********************************************************************
  
  // change sum over iterations to average over iterations ***
  
  if (burnin >-1) {
    
    double diff = iterations - burnin;
    
    // average over chain after burnin 
    for (unsigned int d = 0; d < Nd; d++) { // consider parallelization
      for (unsigned int k = 0; k < Nk; k++) {
        
        Cd_mean(k, d) = ((double)Cd_sum(k, d) / diff);
        
      }
    }
    
    for (unsigned int v = 0; v < Nv; v++) { // consider parallelization
      for (unsigned int k = 0; k < Nk; k++) {
        
        Cv_mean(k, v) = ((double)Cv_sum(k, v) / diff);
        
      }
    }
  }
  
  // Return the final list ***
  
  return List::create(
    Named("Cd") = Cd,
    Named("Cv") = Cv,
    Named("Ck") = Ck,
    Named("Cd_mean") = Cd_mean,
    Named("Cv_mean") = Cv_mean,
    Named("log_likelihood") = log_likelihood,
    Named("alpha") = alpha,
    Named("beta") = beta
  );  
}
