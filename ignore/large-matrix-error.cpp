#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>

#include "../src/parallel_gibbs_utils.h"
#include "../src/sample_int.h"
#include "../src/matrix_conversions.h"

#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
arma::sp_mat transpose(arma::sp_mat x) {
  return x.t();
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(Matrix)

load("../ignore/large-sparse-mat.RData")

x <- transpose(sbir_tcm)
*/
