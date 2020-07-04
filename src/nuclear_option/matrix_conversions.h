// This header converts Rcpp matrices to std::vector<> and back again
#pragma once

#include <RcppArmadillo.h>
using namespace Rcpp;

// convert an NumericMatrix to two-dimensional vector
std::vector<std::vector<double>> mat_to_vec(
    const Rcpp::NumericMatrix& x,
    const bool&          by_rows
) {
  
  NumericMatrix tmp;
  
  if (by_rows) {
    tmp = Rcpp::transpose(x);
  } else {
    tmp = x;
  }
  
  std::vector<std::vector<double>> out(tmp.cols());
  
  for (std::size_t j = 0; j < tmp.cols(); j++) {
    std::vector<double> tmp_col(tmp.rows());
    
    for (std::size_t k = 0; k < tmp.rows(); k++){
      tmp_col[k] = tmp(k, j);
    }
    
    out[j] = tmp_col;
  }
  
  return out;
}

// convert an IntegerMatrix of only positive values to two-dimensional vector
std::vector<std::vector<std::size_t>> mat_to_vec(
    const Rcpp::IntegerMatrix& x,
    const bool&          by_rows
) {
  
  IntegerMatrix tmp;
  
  if (by_rows) {
    tmp = Rcpp::transpose(x);
  } else {
    tmp = x;
  }
  
  std::vector<std::vector<std::size_t>> out(tmp.cols());
  
  for (std::size_t j = 0; j < tmp.cols(); j++) {
    
    std::vector<std::size_t> tmp_col(tmp.rows());
    
    for (std::size_t k = 0; k < tmp.rows(); k++){
      tmp_col[k] = tmp(k, j);
    }
    
    out[j] = tmp_col;
  }
  
  return out;
}

// convert a std::vector to NumericMatrix
NumericMatrix vec_to_mat(
    const std::vector<std::vector<double>>& x,
    const bool&                             row_major
) {
  
  std::size_t n_cols = x.size();
  std::size_t n_rows = x[0].size();
  
  Rcpp::NumericMatrix out(n_rows, n_cols);
  
  for (std::size_t j = 0; j < n_cols; j++) {
    for (std::size_t k = 0; k < n_rows; k++) {
      out(k, j) = x[j][k];
    }
  }
  
  if (row_major) {
    out = Rcpp::transpose(out);
  }
  
  return out;
}

// convert a std::vector of only positive ints to IntegerMatrix
IntegerMatrix vec_to_mat(
    const std::vector<std::vector<std::size_t>>& x,
    const bool&                             row_major
) {
  
  std::size_t n_cols = x.size();
  std::size_t n_rows = x[0].size();
  
  Rcpp::IntegerMatrix out(n_rows, n_cols);
  
  for (std::size_t j = 0; j < n_cols; j++) {
    for (std::size_t k = 0; k < n_rows; k++) {
      out(k, j) = x[j][k];
    }
  }
  
  if (row_major) {
    out = Rcpp::transpose(out);
  }
  
  return out;
}

