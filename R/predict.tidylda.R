### Predict method for LDA objects
#' Get predictions from a Latent Dirichlet Allocation model
#' @description Obtains predictions of topics for new documents from a fitted LDA model
#' @param object a fitted object of class \code{tidylda}
#' @param new_data a DTM or TCM of class \code{dgCMatrix} or a numeric vector
#' @param method one of either "gibbs" or "dot". If "gibbs" Gibbs sampling is used
#'        and \code{iterations} must be specified.
#' @param iterations If \code{method = "gibbs"}, an integer number of iterations
#'        for the Gibbs sampler to run. A future version may include automatic stopping criteria.
#' @param burnin If \code{method = "gibbs"}, an integer number of burnin iterations.
#'        If \code{burnin} is greater than -1, the entries of the resulting "theta" matrix
#'        are an average over all iterations greater than \code{burnin}.
#'        Behavior is the same as documented in \code{\link[tidylda]{tidylda}}.
#' @param no_common_tokens behavior when encountering documents that have no tokens
#'        in common with the model. Options are "\code{default}", "\code{zero}",
#'        or "\code{uniform}". See 'details', below for explanation of behavior. 
#' @param threads Number of parallel threads, defaults to 1.
#' @param ... Additional arguments, currently unused
#' @return a "theta" matrix with one row per document and one column per topic
#' @details 
#'   If \code{predict.tidylda} encounters documents that have no tokens in common
#'   with the model in \code{object} it will engage in one of three behaviors based
#'   on the setting of \code{no_common_tokens}.
#'   
#'   \code{default} (the default) sets all topics to 0 for offending documents. This 
#'   enables continued computations downstream in a way that \code{NA} would not.
#'   However, if \code{no_common_tokens == "default"}, then \code{predict.tidylda}
#'   will emit a warning for every such document it encounters.
#'   
#'   \code{zero} has the same behavior as \code{default} but it emits a message
#'   instead of a warning.
#'   
#'   \code{uniform} sets all topics to 1/k for every topic for offending documents.
#'   it does not emit a warning or message.
#' @examples
#' \donttest{
#' # load some data
#' data(nih_sample_dtm, package = "textmineR")
#'
#' # fit a model
#' set.seed(12345)
#'
#' m <- tidylda(
#'   dtm = nih_sample_dtm[1:20, ], k = 5,
#'   iterations = 200, burnin = 175
#' )
#'
#' str(m)
#'
#' # predict on held-out documents using gibbs sampling "fold in"
#' p1 <- predict(m, nih_sample_dtm[21:100, ],
#'   method = "gibbs",
#'   iterations = 200, burnin = 175
#' )
#'
#' # predict on held-out documents using the dot product method
#' p2 <- predict(m, nih_sample_dtm[21:100, ], method = "dot")
#'
#' # compare the methods
#' barplot(rbind(p1[1, ], p2[1, ]), beside = TRUE, col = c("red", "blue"))
#' }
#' @export
predict.tidylda <- function(
  object, 
  new_data, 
  method = c("gibbs", "dot"),
  iterations = NULL, 
  burnin = -1, 
  no_common_tokens = c("default", "zero", "uniform"),
  threads = 1,
  ...
){
  
  ### Check inputs ----
  if (method[1] == "gibbs") {
    if (is.null(iterations)) {
      stop("when using method 'gibbs' iterations must be specified.")
    }

    if (burnin >= iterations) {
      stop("burnin must be less than iterations")
    }
  }

  # handle dtm
  new_data <- convert_dtm(dtm = new_data)

  if (sum(c("gibbs", "dot") %in% method) == 0) {
    stop("method must be one of 'gibbs' or 'dot'")
  }

  dtm_new_data <- new_data

  if (sum(no_common_tokens %in% c("default", "zero", "uniform")) <= 0) {
    stop(
      "no_common_tokens must be one of 'default', 'zero', or 'uniform'."
    )
  }
  
  ### Align vocabulary ----
  # this is fancy because of how we do indexing in gibbs sampling
  vocab_original <- colnames(object$phi) # tokens in training set

  vocab_intersect <- intersect(vocab_original, colnames(dtm_new_data))

  vocab_add <- setdiff(vocab_original, vocab_intersect)

  add_mat <- Matrix::Matrix(0, nrow = nrow(dtm_new_data), ncol = length(vocab_add))

  colnames(add_mat) <- vocab_add

  dtm_new_data <- Matrix::cbind2(dtm_new_data, add_mat)

  if (nrow(dtm_new_data) == 1) {
    dtm_new_data <- Matrix::Matrix(dtm_new_data[, vocab_original], nrow = 1, sparse = TRUE)

    colnames(dtm_new_data) <- vocab_original

    rownames(dtm_new_data) <- 1
  } else {
    dtm_new_data <- dtm_new_data[, vocab_original]
  }

  ### Get predictions ----

  if (method[1] == "dot") { # dot product method

    result <- dtm_new_data[, vocab_original]

    # handle differently if one row
    if (nrow(dtm_new_data) == 1) {
      result <- result / sum(result)
    } else {
      result <- result / Matrix::rowSums(result)
    }

    result <- result %*% t(object$gamma[, vocab_original])
    result <- as.matrix(result)
    
    repl <- is.na(result)
    
    bad_docs <- which(rowSums(repl) > 0)
    
    rownames(result) <- rownames(dtm_new_data)
    colnames(result) <- rownames(object$phi)
    
    # how do you want to handle empty documents?
    if (no_common_tokens[1] %in% c("default", "zero")) {
      if (length(bad_docs) > 0) {
        result[repl] <- 0
        if (no_common_tokens[1] == "default") {
          for (bad in bad_docs) {
            warning(
              "Document ", bad, " has no tokens in common with the model. ",
              "Setting predictions to 0 for all documents. To change this behavior ",
              "or silence this warning, change the value of 'no_common_tokens' in ",
              "the call to predict.tidylda."
            )
          }
        } else {
          for (bad in bad_docs) {
            message(
              "Document ", bad, " has no tokens in common with the model. ",
              "Setting predictions to 0 for all documents."
            )
          }
        } 
      } 
    } else { # means no_common_tokens == "uniform"
      result[bad_docs, ] <- 1 / ncol(object$theta)
    }
  } else { # gibbs method
    # format inputs

    # get initial distribution with recursive call to "dot" method
    theta_initial <- predict.tidylda(
      object = object, 
      new_data = new_data, 
      method = "dot", 
      no_common_tokens = "uniform"
    )

    # make sure priors are formatted correctly
    beta <- format_beta(object$beta, k = nrow(object$phi), Nv = ncol(dtm_new_data))

    alpha <- format_alpha(object$alpha, k = nrow(object$phi))

    # get initial counts
    lex <- initialize_topic_counts(
      dtm = dtm_new_data,
      k = nrow(object$phi),
      alpha = alpha$alpha,
      beta = beta$beta,
      phi_initial = object$phi,
      theta_initial = theta_initial,
      freeze_topics = TRUE,
      threads = threads,
    )

    # pass inputs to C function for prediciton
    lda <- fit_lda_c(
      docs = lex$docs,
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
      optimize_alpha = FALSE
    )

    # format posterior prediction
    result <- format_raw_lda_outputs(
      lda = lda, 
      dtm = dtm_new_data,
      burnin = burnin, 
      is_prediction = TRUE, 
      threads
    )
  }

  # return result


  result
}
