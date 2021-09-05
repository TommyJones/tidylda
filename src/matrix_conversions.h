// This header converts Rcpp matrices to std::vector<> and back again
#pragma once

using namespace Rcpp;

// convert an NumericMatrix to two-dimensional vector
std::vector<std::vector<double>> mat_to_vec(
    const Rcpp::NumericMatrix& x,
    const bool&          by_rows = false
) {
  
  NumericMatrix tmp;
  
  if (by_rows) {
    tmp = Rcpp::transpose(x);
  } else {
    tmp = x;
  }
  
  std::vector<std::vector<double>> out(tmp.cols());
  
  for (auto j = 0; j < tmp.cols(); j++) {
    std::vector<double> tmp_col(tmp.rows());
    
    for (auto k = 0; k < tmp.rows(); k++){
      tmp_col[k] = tmp(k, j);
    }
    
    out[j] = tmp_col;
  }
  
  return out;
}

// convert an IntegerMatrix to two-dimensional vector
std::vector<std::vector<long>> mat_to_vec(
    const Rcpp::IntegerMatrix& x,
    const bool&          by_rows = false
) {
  
  IntegerMatrix tmp;
  
  if (by_rows) {
    tmp = Rcpp::transpose(x);
  } else {
    tmp = x;
  }
  
  std::vector<std::vector<long>> out(tmp.cols());
  
  for (auto j = 0; j < tmp.cols(); j++) {
    
    std::vector<long> tmp_col(tmp.rows());
    
    for (auto k = 0; k < tmp.rows(); k++){
      tmp_col[k] = tmp(k, j);
    }
    
    out[j] = tmp_col;
  }
  
  return out;
}

// convert a std::vector to NumericMatrix
NumericMatrix vec_to_mat(
    const std::vector<std::vector<double>>& x,
    const bool&                             row_major = false
) {
  
  std::size_t n_cols = x.size();
  std::size_t n_rows = x[0].size();
  
  Rcpp::NumericMatrix out(n_rows, n_cols);
  
  for (auto j = 0; j < n_cols; j++) {
    for (auto k = 0; k < n_rows; k++) {
      out(k, j) = x[j][k];
    }
  }
  
  if (row_major) {
    out = Rcpp::transpose(out);
  }
  
  return out;
}

// convert a std::vector of ints to IntegerMatrix
IntegerMatrix vec_to_mat(
    const std::vector<std::vector<long>>& x,
    const bool&                           row_major = false
) {
  
  std::size_t n_cols = x.size();
  std::size_t n_rows = x[0].size();
  
  Rcpp::IntegerMatrix out(n_rows, n_cols);
  
  for (auto j = 0; j < n_cols; j++) {
    for (auto k = 0; k < n_rows; k++) {
      out(k, j) = x[j][k];
    }
  }
  
  if (row_major) {
    out = Rcpp::transpose(out);
  }
  
  return out;
}

