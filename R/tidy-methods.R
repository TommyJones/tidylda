#' Glance method for \code{tidylda} objects
#' @description
#'   \code{glance} constructs a single-row summary "glance" of a \code{tidylda}
#'   topic model.
#' @param x an object of class \code{tidylda}
#' @param ... other arguments passed to methods,currently not used
#' @return
#'   \code{glance} returns a one-row \code{\link[tibble]{tibble}} with the
#'   following columns:
#'
#'   \code{num_topics}: the number of topics in the model
#'   \code{num_documents}: the number of documents used for fitting
#'   \code{num_tokens}: the number of tokens covered by the model
#'   \code{iterations}: number of total Gibbs iterations run
#'   \code{burnin}: number of burn-in Gibbs iterations run
#' @examples
#' \donttest{
#' dtm <- textmineR::nih_sample_dtm
#'
#' lda <- tidylda(data = dtm, k = 10, iterations = 100, burnin = 75)
#'
#' glance(lda)
#' }
#' @export
glance.tidylda <- function(x, ...) {
  
  # get some objects to return
  if (inherits(x$call, "call")) {
    call_params <- as.list(x$call)
    
    if (!"burnin" %in% names(call_params)) {
      call_params$burnin <- NA
    }
  } else {
    call_params <- list(
      iterations = NA,
      burnin = NA
    )
  }
  
  out <- data.frame(
    num_topics = nrow(x$beta),
    num_documents = nrow(x$theta),
    num_tokens = ncol(x$beta),
    iterations = call_params$iterations,
    burnin = call_params$burnin,
    stringsAsFactors = FALSE
  )
  
  tibble::as_tibble(out)
}

#' Tidy a matrix from a \code{tidylda} topic model
#' @description
#' Tidy the result of a \code{tidylda} topic model
#' @param x an object of class \code{tidylda} or an individual \code{beta}, 
#'   \code{theta}, or \code{lambda} matrix.
#' @param matrix the matrix to tidy; one of \code{'beta'}, \code{'theta'}, or
#'   \code{'lambda'}
#' @param log do you want to have the result on a log scale? Defaults to \code{FALSE}
#' @param ... other arguments passed to methods,currently not used
#' @return
#'   Returns a \code{\link[tibble]{tibble}}.
#'
#'   If \code{matrix = "beta"} then the result is a table of one row per topic
#'   and token with the following columns: \code{topic}, \code{token}, \code{beta}
#'
#'   If \code{matrix = "theta"} then the result is a table of one row per document
#'   and topic with the following columns: \code{document}, \code{topic}, \code{theta}
#'
#'   If \code{matrix = "lambda"} then the result is a table of one row per topic
#'   and token with the following columns: \code{topic}, \code{token}, \code{lambda}
#' @note
#'   If \code{log = TRUE} then "log_" will be appended to the name of the third
#'   column of the resulting table. e.g "\code{beta}" becomes "\code{log_beta}".
#' @examples
#' \donttest{
#' dtm <- textmineR::nih_sample_dtm
#'
#' lda <- tidylda(data = dtm, k = 10, iterations = 100, burnin = 75)
#'
#' tidy_beta <- tidy(lda, matrix = "beta")
#'
#' tidy_theta <- tidy(lda, matrix = "theta")
#'
#' tidy_lambda <- tidy(lda, matrix = "lambda")
#' }
#' @export
tidy.tidylda <- function(x, matrix, log = FALSE, ...) {
  
  if (!inherits(matrix, "character") |
      !sum(c("beta", "theta", "lambda") %in% matrix) >= 1) {
    stop("matrix should be one of c('beta', 'theta', 'lambda')")
  }
  
  if (matrix == "beta") {
    out <- tidy.matrix(x = x$beta, matrix = matrix, log = log)
  } else if (matrix == "lambda") {
    out <- tidy.matrix(x = x$lambda, matrix = matrix, log = log)
  } else {
    out <- tidy.matrix(x = x$theta, matrix = matrix, log = log)
  }
  
  out
}

#' Tidy an individual matrix. Useful for predictions and called from tidy.tidylda
#' @describeIn tidy.tidylda Tidy an individual matrix.
#'   Useful for predictions and called from tidy.tidylda
#' @export
tidy.matrix <- function(x, matrix, log = FALSE, ...) {
  
  # check inputs
  if (!inherits(matrix, "character") |
      !sum(c("beta", "theta", "lambda") %in% matrix) >= 1) {
    stop("matrix should be one of c('beta', 'theta', 'lambda')")
  }
  
  if (!is.logical(log)) {
    stop("log must be logical.")
  }
  
  out <- as.data.frame(x, stringsAsFactors = FALSE)
  
  out$names_col <- rownames(x)
  
  out <- tidyr::pivot_longer(
    data = out, cols = setdiff(colnames(out), "names_col"),
    names_to = "index", values_to = "value"
  )
  
  if (matrix == "beta") {
    
    colnames(out) <- c("topic", "token", "beta")
    
    out$topic <- as.numeric(out$topic)
    
  } else if (matrix == "lambda") {
    
    colnames(out) <- c("topic", "token", "lambda")
    
    out$topic <- as.numeric(out$topic)
    
  } else { # meanse matrix  == theta
    
    colnames(out) <- c("document", "topic", "theta")
    
    out$topic <- as.numeric(stringr::str_replace_all(out$topic, "^X", ""))
    
  }
  
  if (log) {
    out[[3]] <- log(out[[3]])
    
    names(out)[3] <- paste0("log_", names(out)[3])
  }
  
  out
  
}


#' Augment method for \code{tidylda} objects
#' @description
#'   \code{augment} appends observation level model outputs.
#' @param x an object of class \code{tidylda}
#' @param data a tidy tibble containing one row per original document-token pair, 
#'   such as is returned by \link[tidytext]{tdm_tidiers} with column names
#'   c("document", "term") at a minimum.
#' @param type one of either "class" or "prob"
#' @param ... other arguments passed to methods,currently not used
#' @return
#'   \code{augment} returns a tidy tibble containing one row per document-token
#'   pair, with one or more columns appended, depending on the value of \code{type}.
#'   
#'   If \code{type = 'prob'}, then one column per topic is appended. Its value
#'   is P(topic | document, token).
#'   
#'   If \code{type = 'class'}, then the most-probable topic for each document-token
#'   pair is returned. If multiple topics are equally probable, then the topic
#'   with the smallest index is returned by default.
#' @details 
#'   The key statistic for \code{augment} is P(topic | document, token) =
#'   P(topic | token) * P(token | document). P(topic | token) are the entries
#'   of the 'lambda' matrix in the \code{\link[tidylda]{tidylda}} object passed
#'   with \code{x}. P(token | document) is taken to be the frequency of each
#'   token normalized within each document.
#' @export
augment.tidylda <- function(
  x,
  data,
  type = c("class", "prob"),
  ...
) {
  
  # check inputs
  if (sum(c("class", "prob") %in% type[1]) < 1) {
    stop("type must be one of 'class' or 'prob'")
  }
  
  # If they didn't pass a data.frame or if that data.frame doesn't have the right
  # columns, then try to make a dtm , get its relative frequencies, and re-create
  # a tibble with document-token pairs
  if (! inherits(data, "data.frame")) {
    
    dtm <- try(
      convert_dtm(data),
      silent = TRUE
    )
    
    if (inherits(dtm, "try-error")) { # if that fails exit
      stop(
        "data argument must either be a data.frame or tibble with 'document' and 'term' columns ",
        "or coercible to a dgCMatrix. Supported classes are ",
        "c('Matrix', 'matrix', 'simple_triplet_matrix', 'dfm', 'DocumentTermMatrix'), ",
        "However, I see class(data) = ", class(data)
      )
    } 
    
    # get each row to sum to one
    dtm <- dtm / Matrix::rowSums(dtm)
    
    # cast as a triplet matrix
    data <- tidy_dgcmatrix(dtm)
    
    
  }else {
    if (! all(c("document", "term") %in% colnames(data))) {
      stop("data is a data.frame but it does not have c('document' and 'term') colnames)")
    }
    
    # if a tidy tibble, need to get fraction of words in each document
    tmp <- 
      data %>% 
      dplyr::group_by(.data$document, .data$term) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::mutate(count = .data$n / sum(.data$n))
    
    data <- dplyr::left_join(
      data, 
      tmp,
      by = c("document" = "document", "term" = "term")
    )
    
  }
  
  tidy_lambda <- tidy(x, "lambda")
  
  tidy_lambda <- tidyr::pivot_wider(
    tidy_lambda, 
    names_from = .data$topic, 
    values_from = .data$lambda
  )
  
  result <- dplyr::left_join(
    data,
    tidy_lambda,
    by = c("term" = "token")
  )
  
  topic_names <- colnames(x$theta)
  
  result[, topic_names] <- 
    result[, topic_names] * result$count
  
  # return class or probs based on user input
  if (type[1] == "class") {
    tmp <- apply(result[, topic_names], 1,function(y) which.max(y)[1])
    
    result$topic <- tmp
    
    result <- result[, c("document", "term", "topic")]
  } else {
    result <- result[, c("document", "term", topic_names)]
  }
  
  tibble::as_tibble(result)  
}


