#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>


// [[Rcpp::export]]
arma::sp_mat transpose(arma::sp_mat x) {
  arma::sp_mat dtm = x.t(); 
  
  return dtm;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

# note: requires you to manually load ignore.large-sparse-mat.RData

my_transpose <- transpose(sbir_tcm)

*/
