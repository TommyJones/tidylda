# internal functions:

### Makes sure beta is formatted correctly for input to C functions for LDA ----

format_beta <- function(beta, k, Nv) {
  
  if (! is.numeric(beta) | sum(is.na(beta)) > 0 | sum(beta == 0) == length(beta))
    stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
  
  if (length(beta) == 1) { # if beta is a scalar
    
    beta <- matrix(beta, nrow = k, ncol = Nv)
    
    beta_class <- "scalar"
    
  } else if (is.vector(beta)){ # if beta is a vector
    
    if (length(beta) != Nv) # if you didn't specify this vector right
      stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    
    # otherwise let's carry on...
    # make beta a matrix to format for C++ funciton
    beta <- t(beta + matrix(0, nrow = length(beta), ncol = k))
    
    beta_class <- "vector"
    
  } else if (is.matrix(beta)) { # if beta is a matrix
    
    beta_class <- "matrix"
    
  } else { # if beta is of an un supported data type
    
    stop("beta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    
  }
  
  
  list(beta = beta,
       beta_class = beta_class)
}

### Makes sure alpha is formatted correctly for input to C functions for LDA ----
format_alpha <- function(alpha, k) {
  
  if (! is.numeric(alpha) | sum(is.na(alpha)) > 0 | sum(alpha == 0) == length(alpha))
    stop("alpha must be a numeric scalar or vector of length 'k' with no missing 
          values and at least one non-zero value")
  
  if (length(alpha) == 1 & is.numeric(alpha)) {
    alpha <- numeric(k) + alpha
    
    alpha_class <- "scalar"
    
  } else if (length(alpha) != k | ! is.vector(alpha)){
    
    stop("alpha must be a numeric scalar or numeric vector of length 'k'")
    
  } else {
    
    alpha_class <- "vector"
    
  }
  
  list(alpha = alpha,
       alpha_class = alpha_class)
  
}

