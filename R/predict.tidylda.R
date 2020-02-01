### Predict method for LDA objects
#' Get predictions from a Latent Dirichlet Allocation model
#' @description Obtains predictions of topics for new documents from a fitted LDA model
#' @param object a fitted object of class \code{tidylda}
#' @param newdata a DTM or TCM of class \code{dgCMatrix} or a numeric vector
#' @param method one of either "gibbs" or "dot". If "gibbs" Gibbs sampling is used
#'        and \code{iterations} must be specified.
#' @param iterations If \code{method = "gibbs"}, an integer number of iterations 
#'        for the Gibbs sampler to run. A future version may include automatic stopping criteria.
#' @param burnin If \code{method = "gibbs"}, an integer number of burnin iterations. 
#'        If \code{burnin} is greater than -1, the entries of the resulting "theta" matrix
#'        are an average over all iterations greater than \code{burnin}.
#'        Behavior is the same as documented in \code{\link[tidylda]{tidylda}}. 
#' @param ... Other arguments to be passed to \code{\link[furrr]{future_map}}
#' @return a "theta" matrix with one row per document and one column per topic
#' @examples
#' \dontrun{
#' # load some data
#' data(nih_sample_dtm, package = "textmineR")
#' 
#' # fit a model 
#' set.seed(12345)
#' 
#' m <- tidylda(dtm = nih_sample_dtm[1:20,], k = 5,
#'                    iterations = 200, burnin = 175)
#'
#' str(m)
#' 
#' # predict on held-out documents using gibbs sampling "fold in"
#' p1 <- predict(m, nih_sample_dtm[21:100,], method = "gibbs",
#'               iterations = 200, burnin = 175)
#' 
#' # predict on held-out documents using the dot product method
#' p2 <- predict(m, nih_sample_dtm[21:100,], method = "dot")
#'
#' # compare the methods
#' barplot(rbind(p1[1,],p2[1,]), beside = TRUE, col = c("red", "blue")) 
#' }
#' @export
predict.tidylda <- function(object, newdata, method = c("gibbs", "dot"), 
                            iterations = NULL, burnin = -1, ...) {
  
  ### Check inputs ----
  if (method[1] == "gibbs") {
    
    if (is.null(iterations)) {
      stop("when using method 'gibbs' iterations must be specified.")
    }
    
    if (burnin >= iterations) {
      stop("burnin must be less than iterations")
    }
    
  }
  
  if (sum(c("dgCMatrix", "numeric") %in% class(newdata)) < 1) {
    stop("newdata must be a matrix of class dgCMatrix or a numeric vector")
  }
  
  if (class(newdata) == "numeric") { # if newdata is a numeric vector, assumed to be 1 document
    if (is.null(names(newdata))) {
      stop("it looks like newdata is a numeric vector without names. 
           Did you mean to pass a single document?
           If so, it needs a names attribute to index tokens")
    }
    
    
    vocab <- names(newdata)
    
    newdata <- Matrix::Matrix(newdata, nrow = 1, sparse = TRUE)
    
    colnames(newdata) <- vocab
    
    rownames(newdata) <- 1
    
  } else { # assume you passed something that can be converted to dgCMatrix
    
    # Ensure dtm is of class dgCMatrix
    newdata <- convert_dtm(dtm = newdata)
    
    
  }
  
  if (sum(c("gibbs", "dot") %in% method) == 0) {
    stop("method must be one of 'gibbs' or 'dot'")
  }
  
  dtm_newdata <- newdata
  
  ### Align vocabulary ----
  # this is fancy because of how we do indexing in gibbs sampling
  vocab_original <- colnames(object$phi) # tokens in training set
  
  vocab_intersect <- intersect(vocab_original, colnames(dtm_newdata))
  
  vocab_add <- setdiff(vocab_original, vocab_intersect)
  
  add_mat <- Matrix::Matrix(0, nrow = nrow(dtm_newdata), ncol = length(vocab_add))
  
  colnames(add_mat) <- vocab_add
  
  dtm_newdata <- Matrix::cbind2(dtm_newdata, add_mat)
  
  if (nrow(dtm_newdata) == 1) {
    dtm_newdata <- Matrix::Matrix(dtm_newdata[,vocab_original], nrow = 1, sparse = TRUE)
    
    colnames(dtm_newdata) <- vocab_original
    
    rownames(dtm_newdata) <- 1
  } else {

    dtm_newdata <- dtm_newdata[, vocab_original]
    
  }
  
  ### Get predictions ----
  
  if (method[1] == "dot") { # dot product method
    
    result <- dtm_newdata[, vocab_original]
    
    # handle differently if one row
    if (nrow(dtm_newdata) == 1) {
      result <- result / sum(result)
    } else {
      result <- result / Matrix::rowSums(result)
    }
    
    result <- result %*% t(object$gamma[ ,vocab_original])
    result <- as.matrix(result)
    result[is.na(result)] <- 0
    
    rownames(result) <- rownames(dtm_newdata)
    colnames(result) <- rownames(object$phi)
    
    
  } else { # gibbs method
    # format inputs
    
    # get initial distribution with recursive call to "dot" method
    theta_initial <- predict.tidylda(object = object, newdata = newdata, method = "dot")
    
    # make sure priors are formatted correctly
    beta <- format_beta(object$beta, k = nrow(object$phi), Nv = ncol(dtm_newdata))
    
    alpha <- format_alpha(object$alpha, k = nrow(object$phi))
    
    # get initial counts
    lex <- initialize_topic_counts(dtm = dtm_newdata,
                                   k = nrow(object$phi),
                                   alpha = alpha$alpha,
                                   beta = beta$beta,
                                   phi_initial = object$phi,
                                   theta_initial = theta_initial,
                                   freeze_topics = TRUE,
                                   ...)
    
    # pass inputs to C function for prediciton
    lda <- fit_lda_c(docs = lex$docs,
                       Nk = nrow(object$phi),
                       beta = beta$beta,
                       alpha = alpha$alpha,
                       Cd = lex$Cd,
                       Cv = lex$Cv,
                       Ck = lex$Ck,
                       Zd = lex$Zd,
                       Phi = object$phi,
                       iterations = iterations,
                       burnin = burnin,
                       freeze_topics = TRUE,
                       calc_likelihood = FALSE,
                       optimize_alpha = FALSE)
    
    # format posterior prediction
    result <- format_raw_lda_outputs(lda = lda, dtm = dtm_newdata, 
                                     burnin = burnin, is_prediction = TRUE, ...)
    
  }
  
  # return result
  
  
  result
  
}
