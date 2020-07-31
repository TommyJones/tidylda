// This file is a header file for a replacement to RcppArmadillo::sample
// It is modified by code written by Barum Park, originally using Armadillo. 
// Unmodified source: https://barumpark.com/blog/2019/Sampling-Integers/#pragma once

#include <algorithm>
#include <numeric>
#include <cmath>

// make a categorical struct that maps a category to its probability
// function below makes a vector of categoricals
// Need to sort based on probabilities. This is more easily done in a struct
struct categorical
{
  int category;
  double probability;
  
  // declare constructor
  categorical(){};
  
  // define constructor
  categorical(
    int init_cat, 
    double init_prob
  ) {
    category = init_cat;
    probability = init_prob;
  }
};

// vector of categorical variables for sorting
typedef std::vector<categorical> categorical_variable;

// The function below does not use a Random number generator. Instead, you take
// a sample from a uniform(0,1) RV and pass it in the 'u' argument. This lets
// you use any RNG you'd like. 
int sample_one(
    const std::vector<double> &p, // vector of probabilities inexing categories 0 to k
    double                    &u  // sample from a U(0,1) RV
) {
  
  const std::size_t Nk = p.size();
  
  // check inputs
  if (std::any_of(p.begin(), p.end(), [](double j){ return j < 0.0;})) {
    throw std::runtime_error("Negative probabilities not allowed"); 
  }
  
  // construct categorical that combines index and probability
  // while we're at it, make sure probability is normalized
  categorical_variable cat_var;
  
  for (auto k = 0; k < Nk; k++) {
    cat_var.push_back(categorical(k, p[k]));
  }
  
  // sort our categorical from highest probability to lowest probability
  std::sort(cat_var.begin(), cat_var.end(), 
            [](categorical const &a, categorical const &b) { 
              return a.probability < b.probability; 
            }); 
  
  // generate a q vector of cumulate probabilities
  std::vector<double> q(Nk);
  
  q[0] = cat_var[0].probability;
  
  for (auto k = 1; k < Nk; k++) {
    q[k] = cat_var[k].probability + q[k - 1];
  }
  
  // compare to u and return the selected value
  u *= q.back();
  
  for (auto k = 0; k < Nk; k++) {
    if (u <= q[k])
      return cat_var[k].category;
  }
  
  throw std::runtime_error("couldn't find index (samp_one)");
}

// a function to calculate the logarithm of the sum of two exponentials
inline double log_add_exp(const double x, const double y) {
  
  double d = x - y;
  if (d > 0.0)
    return x + std::log1p(std::exp(-d));
  if (d <= 0.0)
    return y + std::log1p(std::exp(d));
  
  return 1.0; // if you get here, you've messed up
  
}

// function to sample a categorical variable using log probabilities
int log_sample_one(
    const std::vector<double> &log_p, // vector of log probabilities inexing categories 0 to k
    double                    &e      // sample from an exponential(1.0) RV
) {
  
  const std::size_t Nk = log_p.size();
  
  // check whether all elements are finite and numbers
  if (std::any_of(log_p.begin(), log_p.end(), [](double j){ 
    return std::isinf(j);
  })) {
    throw std::runtime_error("Infinite log probability detected!"); 
  }
  
  if (std::any_of(log_p.begin(), log_p.end(), [](double j){ 
    return std::isnan(j);
  })) {
    throw std::runtime_error("NaN log probability detected!"); 
  }
  
  // construct categorical that combines index and probability
  categorical_variable cat_var;
  
  for (auto k = 0; k < Nk; k++) {
    cat_var.push_back(categorical(k, log_p[k]));
  }
  
  // sort our categorical from highest probability to lowest probability
  std::sort(cat_var.begin(), cat_var.end(), 
            [](categorical const &a, categorical const &b) { 
              return a.probability < b.probability; 
            }); 
  
  // generate a q vector of cumulate probabilities
  std::vector<double> q(Nk);
  
  for (auto k = 0; k < Nk; k++) {
    q[k] = log_add_exp(q[k - 1], cat_var[k].probability);
  }

  // compare to u and return the selected value
  double u = q.back();
  
  u -= e; // u, because it's a uniform rv in log space
  
  for (auto k = 0; k < Nk; k++) {
    if (u <= q[k])
      return cat_var[k].category;
  }
  
  throw std::runtime_error("couldn't find index (log_sample_one)");
}



