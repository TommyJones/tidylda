// This file is a header file for a replacement to RcppArmadillo::sample
// It was predominantly written by Barum Park. 
// Unmodified source: https://barumpark.com/blog/2019/Sampling-Integers/

#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rmath.h>

// a function to calculate the logarithm of the sum of two exponentials
inline double log_add_exp(const double x, const double y) {
  
  double d = x - y;
  if (d > 0.0)
    return x + std::log1p(std::exp(-d));
  if (d <= 0.0)
    return y + std::log1p(std::exp(d));
  
  return 1.0; // if you get here, you've messed up
  
}

// function to accumulate log-probabilities on the log-scale
template <typename T>
inline T log_accu_exp (const T& x) {
  
  // initialize with x
  T y(x);
  
  // accumulate probabilities on the log-scale
  typename T::iterator i = y.begin() + 1;
  for (; i < y.end(); i++) {
    *i = log_add_exp(*(i - 1), *i);
  }
  
  return y;
}

// Sample a single object from a vector of integers using log probabilities
arma::uword lsamp_one(const arma::vec &lpvec) {
  
  //check whether all elements are finite
  if (lpvec.has_inf())
    Rcpp::stop("log-probabilities have to be finite");
  
  if (lpvec.has_nan())
    Rcpp::stop("log-probability vector contains a missing value");
  
  // get indices of elements in descending order
  arma::uvec indx = arma::sort_index(lpvec, "descend");
  // logarithm of accumulated probabilities
  arma::vec lqvec = log_accu_exp(arma::vec(arma::sort(lpvec, "descend")));
  
  // sample log(unif(0,sum(exp(pvec))))
  double u = arma::conv_to<double>::from(lqvec.tail(1L));
  u -= R::rexp(1.0);
  
  // sample integer
  for (arma::uword k = 0; k < lqvec.n_elem; ++k) {
    if (u <= lqvec(k))
      return indx(k);
    
  }
  
  Rcpp::stop("couldn't find index (lsamp_one)");
  
}

// Sample a single object from a vector of integers using probabilities
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
  u *= arma::randu(); 
  
  // find interval into which u falls
  for (arma::uword k = 0; k < pvec.n_elem; ++k) {
    
    if (u <= qvec(k))
      return indx(k);
    
  }
  
  Rcpp::stop("couldn't find index (samp_one)");
  
}
