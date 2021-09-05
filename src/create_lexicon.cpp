// Create a formatted object for iterating a collapsed gibbs sampler

#define ARMA_64BIT_WORD 1

#include "matrix_conversions.h" 
#include "sample_int.h" // RcppArmadillo.h included here




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

