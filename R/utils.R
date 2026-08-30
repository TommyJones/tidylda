################################################################################
# Functions in this file are internal to tidylda
################################################################################

#' Convert various things to a \code{dgCMatrix} to work with various functions
#' and methods
#' @keywords internal
#' @description
#'   Presently, \code{tidylda} makes heavy usage of the \code{dgCMatrix} class.
#'   However, a user may have created a DTM (or TCM) in one of several classes.
#'   Since data could be in several formats, this function converts them to a
#'   \code{dgCMatrix} before passing them along.
#' @param dtm the data you want to convert
#' @return an object of class \code{dgCMatrix}
convert_dtm <- function(dtm) {
  if (inherits(dtm, "matrix")) { # regular R matrix
    
    out <- methods::as(dtm, "dgCMatrix", strict = TRUE)
    
  } else if (inherits(dtm, "Matrix")) { # Matrix class matrix 
    
    out <- methods::as(
      methods::as(
        methods::as(dtm, "dMatrix"), 
        "generalMatrix"
        ), 
      "CsparseMatrix"
      )
    
  }else if (inherits(dtm, "simple_triplet_matrix")) { # triplet matrix
    
    out <- Matrix::sparseMatrix(
      i = dtm$i,
      j = dtm$j,
      x = dtm$v,
      dims = c(dtm$nrow, dtm$ncol),
      dimnames = list(
        rownames = dtm$dimnames$Docs,
        colnames = dtm$dimnames$Terms
      )
    )
    
  } else if (inherits(dtm, "numeric")) { # numeric vector
    
    if (is.null(names(dtm))) {
      stop(
        "it looks like dtm (or new_data if you called 'predict') is a numeric ",
        "vector without names. Did you mean to pass a single document? If so, ",
        "it needs a names attribute to index tokens"
      )
    }
    
    vocab <- names(dtm)
    
    out <- Matrix::Matrix(dtm, nrow = 1, sparse = TRUE)
    
    colnames(out) <- vocab
    
    rownames(out) <- 1
    
  } else {
    
    stop(
      "dtm (or new_data if you called 'predict') cannot be converted to dgCMatrix. Supported classes are ",
      "c('Matrix', 'matrix', 'simple_triplet_matrix', 'dfm', 'DocumentTermMatrix'), ",
      "However, I see class(dtm) = ", class(dtm)
    )
    
  }
  
  out
}

#' Format \code{eta} for input into the sampler
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{eta} but the C++
#'   sampler in \code{\link[tidylda]{fit_lda_warp}} only takes it one way. This function does the
#'   appropriate formatting. It also returns errors if the user input a malformatted
#'   \code{eta}.
#' @param eta the prior for words over topics. Can be a numeric scalar, numeric
#'   vector, or numeric matrix.
#' @param k the number of topics.
#' @param Nv the total size of the vocabulary as inherited from \code{ncol(dtm)}
#'   in \code{\link[tidylda]{tidylda}}.
#' @return
#'   Returns a list with two elements: \code{eta} and \code{eta_class}.
#'   \code{eta} is the post-formatted version of \code{eta} in the form of a
#'   \code{k} by \code{Nv} numeric matrix. \code{eta_class} is a character
#'   denoting whether or not the user-supplied \code{eta} was a "scalar",
#'   "vector", or "matrix".
format_eta <- function(eta, k, Nv) {
  if (!is.numeric(eta) | sum(is.na(eta)) > 0 | sum(eta == 0) == length(eta)) {
    stop("eta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
  }
  
  if (length(eta) == 1) { # if eta is a scalar
    
    # D20: left as a scalar rather than materialized to k by Nv. At k = 500 and
    # Nv = 1e5 that matrix is 400 MB to hold one repeated number, and the engine
    # would then keep a 200 MB float copy of it. Consumers that genuinely need
    # the matrix call eta_matrix(); the sampler and most of the R side do not.
    eta <- as.numeric(eta)
    
    eta_class <- "scalar"
  } else if (is.vector(eta)) { # if eta is a vector
    
    if (length(eta) != Nv) { # if you didn't specify this vector right
      stop("eta must be a numeric scalar, a numeric vector of length 'ncol(dtm)', or
         a numeric matrix with 'k' rows and 'ncol(dtm)' columns with no missing 
         values and at least one non-zero value.")
    }
    
    # otherwise let's carry on...
    # make eta a matrix to format for C++ funciton
    eta <- t(eta + matrix(0, nrow = length(eta), ncol = k))
    
    eta_class <- "vector"
  } else if (is.matrix(eta)) { # if eta is a matrix
    
    # check dims before moving on
    if (nrow(eta) != k | ncol(eta) != Nv) {
      stop(
        "If eta is a matrix, it must have the same number of rows as topics ",
        "and it must have the same number of columns (tokens) as your dtm. ",
        "But I see nrow(eta) = ", nrow(eta), " and k = ", k, ". I also see ",
        "ncol(eta) = ", ncol(eta), " but ncol(dtm) = ", Nv
      )
    }
    
    eta_class <- "matrix"
  }
  
  
  list(
    eta = eta,
    eta_class = eta_class
  )
}

#' Materialize \code{eta} as a topics-by-tokens matrix
#' @keywords internal
#' @description
#'   Since D20, \code{\link[tidylda]{format_eta}} leaves a scalar prior as a
#'   scalar rather than expanding it to \code{k} by \code{Nv}. Call this where a
#'   full matrix is genuinely required. Most call sites do not need one --
#'   arithmetic against a scalar recycles correctly, and the sampler takes the
#'   scalar directly.
#' @param eta a list as returned by \code{\link[tidylda]{format_eta}}
#' @param k the number of topics
#' @param Nv the size of the vocabulary
#' @return a numeric matrix with \code{k} rows and \code{Nv} columns.
eta_matrix <- function(eta, k, Nv) {
  if (is.matrix(eta$eta)) {
    return(eta$eta)
  }
  matrix(eta$eta, nrow = k, ncol = Nv)
}

#' Row sums of \code{eta} without materializing it
#' @keywords internal
#' @description
#'   \code{rowSums(eta_matrix(...))} allocates a dense \code{k} by \code{Nv}
#'   matrix only to collapse it immediately. That is 8 GB at \code{k = 1000},
#'   \code{Nv = 1e6}. A scalar prior gives every row the same sum, so it needs
#'   no matrix at all; a matrix prior is summed directly, since
#'   \code{eta_matrix()} would have returned it unchanged anyway.
#'
#'   NOT BIT-IDENTICAL to the materialized form, and deliberately so.
#'   \code{rowSums()} accumulates in long double, so it returns
#'   1500.0000000000002274 where \code{Nv * eta} returns exactly 1500 --- a
#'   relative difference of 1.5e-16. The engine stores \code{eta} as
#'   \code{float} (D5), whose resolution is nine orders of magnitude coarser,
#'   so the difference is erased before the sampler sees it. The only trace is
#'   the 16th significant digit of a refitted model's \code{eta} slot.
#' @param eta a list as returned by \code{\link[tidylda]{format_eta}}
#' @param k the number of topics
#' @param Nv the size of the vocabulary
#' @return a numeric vector of length \code{k}.
eta_row_sums <- function(eta, k, Nv) {
  if (is.matrix(eta$eta)) {
    return(rowSums(eta$eta))
  }
  rep(Nv * eta$eta, k)
}

#' What this R session can plausibly allocate
#' @keywords internal
#' @description
#'   Best effort, and advisory only: the value is used to make an error message
#'   more informative, never to decide whether to raise one. A limit that moved
#'   with the machine would make the same script succeed on one box and fail on
#'   another, which is worse than a predictable ceiling.
#'
#'   \code{mem.maxVSize()} (base) is finite on macOS, where it is the cap behind
#'   "vector memory limit reached", and infinite elsewhere. On Linux the
#'   fallback is \code{MemAvailable} from \code{/proc/meminfo}, which reports
#'   the HOST inside a container and so may overstate what is really available.
#'   Windows gets neither and is simply omitted from the message.
#' @return a number of bytes, or \code{NA_real_} if nothing could be determined.
session_memory_limit <- function() {
  v <- tryCatch(mem.maxVSize(), error = function(e) Inf)

  if (is.finite(v)) {
    return(v * 1024^2) # mem.maxVSize() is in Mb
  }

  if (file.exists("/proc/meminfo")) {
    info <- tryCatch(
      readLines("/proc/meminfo", warn = FALSE),
      error = function(e) character(0)
    )

    avail <- grep("^MemAvailable:", info, value = TRUE)

    if (length(avail) == 1) {
      kb <- suppressWarnings(as.numeric(gsub("\\D", "", avail)))

      if (!is.na(kb)) {
        return(kb * 1024)
      }
    }
  }

  NA_real_
}

#' Refuse to build a result that cannot fit in memory
#' @keywords internal
#' @description
#'   \code{\link[tidylda]{posterior.tidylda}} and
#'   \code{\link[tidylda]{tidy.tidylda}} both return one row per cell of
#'   something that grows as \code{k * V}, so an innocuous-looking call can ask
#'   for an object of billions of rows. Previously such a call simply exhausted
#'   the session. This raises an error that says how large the result would be
#'   and what to do instead.
#'
#'   The ceiling is a fixed 1 GB, which
#'   \code{options(tidylda.max_result_size = <bytes>)} can raise. Fixed rather
#'   than derived from free memory, so that behavior is reproducible across
#'   machines; it is meant to catch a pathological request, not to track RAM.
#' @param n_rows,n_cols dimensions of the result that would be built
#' @param what character, the thing being built, for the message
#' @param suggestion character, a concrete alternative to offer the caller
#' @return \code{invisible(NULL)}, or an error.
check_result_size <- function(n_rows, n_cols, what, suggestion) {
  # Doubles are 8 bytes and character columns hold 8-byte pointers into R's
  # global string cache, so 8 per cell is a fair estimate either way.
  est <- as.numeric(n_rows) * as.numeric(n_cols) * 8

  max_size <- getOption("tidylda.max_result_size", 1024^3)

  if (est <= max_size) {
    return(invisible(NULL))
  }

  # Adaptive units, so a small ceiling does not report "0.0 GB above the 0.0 GB
  # limit".
  fmt <- function(x) {
    if (x >= 1024^3) {
      paste0(format(round(x / 1024^3, 1), nsmall = 1), " GB")
    } else if (x >= 1024^2) {
      paste0(round(x / 1024^2), " MB")
    } else {
      paste0(round(x / 1024), " KB")
    }
  }

  # Expressed as a multiple of a unit rather than a raw byte count, so it is
  # something a caller can reasonably retype.
  headroom <- if (est >= 1024^3) {
    paste0(ceiling(est * 1.5 / 1024^3), " * 1024^3")
  } else {
    paste0(ceiling(est * 1.5 / 1024^2), " * 1024^2")
  }

  limit <- session_memory_limit()

  stop(
    what, " would need about ", fmt(est), " (",
    format(n_rows, big.mark = ",", scientific = FALSE), " rows x ", n_cols,
    " columns), above the ", fmt(max_size), " limit.",
    if (!is.na(limit)) paste0(" This session looks able to allocate ", fmt(limit), "."),
    "\n  ", suggestion,
    "\n  To raise the limit: options(tidylda.max_result_size = ", headroom, ")",
    call. = FALSE
  )
}

#' Pad a document term matrix with empty columns for missing vocabulary
#' @keywords internal
#' @description
#'   Both \code{\link[tidylda]{refit.tidylda}} and
#'   \code{\link[tidylda]{predict.tidylda}} align a new DTM against a model's
#'   vocabulary by appending all-zero columns for terms the data lacks. This is
#'   that operation, in one place.
#'
#'   THE FILLER MUST BE SPARSE. Until 0.1.0 \code{refit()} built it with
#'   \code{matrix(0, ...)}, a dense allocation of \code{nrow(dtm)} by
#'   \code{length(add)} doubles --- 5.4 GB on a 48,508-document corpus with
#'   15,000 model-only terms, enough to exhaust a 32 GB session before sampling
#'   began. The waste was total, since \code{cbind()} of a sparse and a dense
#'   matrix returns a \code{dgCMatrix} regardless: the dense block existed only
#'   as an argument. \code{predict()} always did this correctly, which is why
#'   the two are now one function.
#' @param dtm a document term matrix of class \code{dgCMatrix}
#' @param add character vector of column names to append, possibly empty
#' @return \code{dtm} with one all-zero column per entry of \code{add},
#'   appended in order, with row and column names preserved.
pad_vocabulary <- function(dtm, add) {
  if (length(add) == 0) {
    return(dtm)
  }

  filler <- Matrix::Matrix(0, nrow = nrow(dtm), ncol = length(add))

  colnames(filler) <- add

  Matrix::cbind2(dtm, filler)
}

#' Format \code{alpha} for input into the sampler
#' @keywords internal
#' @description
#'   There are a bunch of ways users could format \code{alpha} but the C++
#'   sampler in \code{\link[tidylda]{fit_lda_warp}} only takes it one way. This function does the
#'   appropriate formatting. It also returns errors if the user input a malformatted
#'   \code{alpha}.
#' @param alpha the prior for topics over documents. Can be a numeric scalar or
#'   numeric vector.
#' @param k the number of topics.
#' @return
#'   Returns a list with two elements: \code{alpha} and \code{alpha_class}.
#'   \code{alpha} is the post-formatted version of \code{alpha} in the form of a
#'   \code{k}-length numeric vector. \code{alpha_class} is a character
#'   denoting whether or not the user-supplied \code{alpha} was a "scalar" or
#'   "vector".
format_alpha <- function(alpha, k) {
  if (!is.numeric(alpha) | sum(is.na(alpha)) > 0 | sum(alpha == 0) == length(alpha)) {
    stop("alpha must be a numeric scalar or vector of length 'k' with no missing 
          values and at least one non-zero value")
  }
  
  if (length(alpha) == 1 & is.numeric(alpha)) {
    alpha <- numeric(k) + alpha
    
    alpha_class <- "scalar"
  } else if (length(alpha) != k | !is.vector(alpha)) {
    stop("alpha must be a numeric scalar or numeric vector of length 'k'")
  } else {
    alpha_class <- "vector"
  }
  
  list(
    alpha = alpha,
    alpha_class = alpha_class
  )
}


#' Prepare the priors the sampler initializes from
#' @keywords internal
#' @description
#'   Implementing seeded (or guided) LDA models and transfer learning means that
#'   we can't initialize topics with a uniform-random start. This function
#'   prepares the two matrices the sampler needs in order to draw an informed
#'   starting assignment: \code{beta_initial}, giving P(token|topic), and
#'   \code{Cd_start}, the expected number of tokens each topic accounts for in
#'   each document. In the event that you aren't using fancy seeding or transfer
#'   learning, this makes a random initialization by sampling from Dirichlet
#'   distributions parameterized by priors \code{alpha} and \code{eta}.
#'
#'   The per-token work of building the token structure and sampling each
#'   token's starting topic happens inside
#'   \code{\link[tidylda]{fit_lda_warp}}, so nothing proportional to the token
#'   count crosses the R/C++ boundary.
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @param k the number of topics
#' @param alpha the numeric vector prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_alpha}}
#' @param eta the numeric matrix prior for topics over documents as formatted
#'   by \code{\link[tidylda]{format_eta}}
#' @param beta_initial if specified, a numeric matrix for the probability of tokens
#'   in topics. Must be specified for predictions or updates as called by
#'   \code{\link[tidylda]{predict.tidylda}} or \code{\link[tidylda]{refit.tidylda}}
#'   respectively.
#' @param theta_initial if specified, a numeric matrix for the probability of
#'   topics in documents. Must be specified for updates as called by
#'   \code{\link[tidylda]{refit.tidylda}}
#' @param freeze_topics if \code{TRUE} does not update counts of tokens in topics.
#'   This is \code{TRUE} for predictions.
#' @param threads number of parallel threads, currently unused
#' @param ... Additional arguments, currently unused
#' @return
#'   Returns a list with two elements, both of which the engine consumes
#'   directly:
#'
#'   \code{beta_initial} is a numeric matrix with one row per topic and one
#'   column per token, giving P(token|topic) to initialize from. Supplied by the
#'   caller for updates and predictions; sampled from \code{eta} otherwise.
#'
#'   \code{Cd_start} is a numeric matrix, documents by topics, holding
#'   \code{theta_initial * rowSums(dtm)} --- the expected number of tokens each
#'   topic accounts for in each document.
#'
#'   Together these define the informed initialization: the engine samples each
#'   token's starting topic from P(z) proportional to
#'   \code{beta_initial[k, v] * (Cd_start[d, k] + alpha[k])}, in log space. That
#'   per-token work happens in C++ (see \code{\link[tidylda]{fit_lda_warp}}), so
#'   nothing proportional to the token count crosses the R/C++ boundary.
initialize_topic_counts <- function(
  dtm, 
  k, 
  alpha, 
  eta, 
  beta_initial = NULL,
  theta_initial = NULL, 
  freeze_topics = FALSE,
  threads = 1,
  ...
) {
  
  # check inputs
  
  if (! is.numeric(threads)) {
    stop("threads must be an integer 1 or greater")
  } else if (threads < 1) {
    stop("threads must be an integer 1 or greater")
  } else {
    threads = as.integer(threads) # ignore decimal inputs
  }
  
  # initialize beta if not already specified
  # this beta is used to sample topics for inital counts in the C++ function
  if (is.null(beta_initial)) {
    # beta_initial <- gtools::rdirichlet(n = k, alpha = eta)
    
    # One rdirichlet call per topic, in topic order. Kept as a loop rather than a
    # single rdirichlet(n = k, ...) so the RNG is consumed exactly as before --
    # and written to accept a scalar eta without materializing k by Nv (D20).
    # Rows are drawn one at a time and written straight into their final
    # position. The previous version first built a list of all k rows --
    # lapply(seq_len(nrow(eta)), function(i) eta[i, ]) duplicates the entire
    # prior, +283 MB against a 191 MB matrix -- and then let vapply produce a
    # V by k result that had to be transposed, a second k by V copy. A scalar
    # prior was already fine, since rep(list(v), k) shares one vector.
    #
    # RNG ORDER IS LOAD-BEARING and unchanged: still exactly one
    # rdirichlet(n = 1, .) per topic, in topic order.
    beta_initial <- matrix(0, nrow = k, ncol = ncol(dtm))
    
    for (i in seq_len(k)) {
      eta_i <- if (is.matrix(eta)) eta[i, ] else rep(eta, ncol(dtm))
      
      beta_initial[i, ] <- as.numeric(gtools::rdirichlet(n = 1, alpha = eta_i)) +
        .Machine$double.eps # avoid underflow
    }
  }
  
  # initialize theta if not already specified
  # if not specified (e.g. if this is a new model) make a matrix by sampling
  # from alpha.
  if (is.null(theta_initial)) {
    theta_initial <- gtools::rdirichlet(n = nrow(dtm), alpha = alpha) + 
      .Machine$double.eps # avoid underflow
  }
  
  # Cd_start is the expected number of tokens each topic accounts for in each
  # document. Cv is not needed here: the engine derives it from the initial
  # assignment it samples.
  
  
  Cd_start <- theta_initial * Matrix::rowSums(dtm)
  
  # Phase 4 (D16): the per-token work -- building the token structure and
  # sampling each token's initial topic -- now happens inside the engine, in the
  # same walk of the DTM that the sampler uses. What used to come back here as
  # `Docs` and `Zd` (16 bytes per token, out to R and straight back in) is never
  # materialized. This function keeps only the parts that are O(KV + DK) rather
  # than per-token, and hands the engine its two starting matrices.
  list(
    beta_initial = beta_initial,
    Cd_start = Cd_start
  )
}

#' Summarize a topic model consistently across methods/functions
#' @keywords internal
#' @description
#'   Summarizes topics in a model. Called by \code{\link[tidylda]{tidylda}}
#'   and \code{\link[tidylda]{refit.tidylda}} and used to augment
#'   \code{\link[tidylda]{print.tidylda}}.
#' @param theta numeric matrix whose rows represent P(topic|document)
#' @param beta numeric matrix whose rows represent P(token|topic)
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}.
#' @return
#'   Returns a \code{\link[tibble]{tibble}} with the following columns:
#'   \code{topic} is the integer row number of \code{beta}.
#'   \code{prevalence} is the frequency of each topic throughout the corpus it
#'     was trained on normalized so that it sums to 100.
#'   \code{coherence} makes a call to \code{\link[tidylda]{calc_prob_coherence}}
#'     using the default 5 most-probable terms in each topic.
#'   \code{top_terms} displays the top 5 most-probable terms in each topic.
#' @note
#'   \code{prevalence} should be proportional to P(topic). It is calculated by
#'   weighting on document length. So, topics prevalent in longer documents get
#'   more weight than topics prevalent in shorter documents. It is calculated
#'   by
#'
#'   \code{prevalence <- rowSums(dtm) * theta \%>\% colSums()}
#'
#'   \code{prevalence <- (prevalence * 100) \%>\% round(3)}
#'
#'   An alternative calculation (not implemented here) might have been
#'
#'   \code{prevalence <- colSums(dtm) * t(beta) \%>\% colSums()}
#'
#'   \code{prevalence <- (prevalence * 100) \%>\% round(3)}
summarize_topics <- function(theta, beta, dtm) {
  
  # probabilistic coherence with default value for m
  if (nrow(dtm) == 1) {
    coherence <- rep(NA, nrow(beta))
  } else {
    coherence <- calc_prob_coherence(beta = beta, data = dtm)
  }
  
  # prevalence of each topic, weighted by terms
  #
  # crossprod, not `Matrix::rowSums(dtm) * theta` then colSums(). The latter
  # materializes a full D x K dense matrix --- 57 MB at D = 70k, K = 100 --- and
  # then immediately collapses it to a length-K vector. The matrix-vector
  # product computes the same thing with no intermediate.
  prevalence <- as.numeric(crossprod(theta, Matrix::rowSums(dtm)))

  prevalence <- prevalence / sum(prevalence)

  prevalence <- round(prevalence * 100, 2)
  
  # top 5 terms
  #
  # PARTIAL SORT. order() sorts all V entries of a row to keep five, which is
  # O(V log V) where selection is O(V): 4.06 s -> 0.82 s per topic at V = 1e6,
  # so roughly an hour off a k = 1000 model. sort(partial =) finds the 5th
  # largest without ordering anything below it; only the handful of values at
  # or above it are then sorted.
  #
  # Exact ties at the threshold yield more than five candidates, and the order
  # among tied values can differ from order()'s. These are display terms for a
  # printed summary, so that is acceptable --- but it is why this is not
  # guaranteed byte-for-byte identical to the previous implementation.
  n_top <- 5
  
  top_terms <- apply(beta, 1, function(x) {
    if (length(x) <= n_top) {
      return(names(x)[order(x, decreasing = TRUE)][seq_len(n_top)])
    }
    
    cut <- length(x) - n_top + 1
    
    keep <- which(x >= sort(x, partial = cut)[cut])
    
    names(x)[keep[order(x[keep], decreasing = TRUE)]][seq_len(n_top)]
  })
  
  top_terms <- apply(top_terms, 2, function(x) {
    paste(c(x, "..."), collapse = ", ")
  })
  
  # combine into a summary
  summary <- data.frame(
    topic = as.numeric(rownames(beta)),
    prevalence = prevalence,
    coherence = coherence,
    top_terms = top_terms,
    stringsAsFactors = FALSE
  )
  
  summary <- tibble::as_tibble(summary)
  
  summary
}

#' Read \code{counts$Cv} in the orientation this version expects
#' @keywords internal
#' @description
#'   D17 changed the exported word-topic counts from topics-by-words to
#'   words-by-topics, matching the orientation the engine holds them in. Models
#'   fitted by earlier versions carry the old shape, and both
#'   \code{\link[tidylda]{refit.tidylda}} and
#'   \code{\link[tidylda]{posterior.tidylda}} index this matrix directly.
#'
#'   A wrong orientation on the transfer-learning path would corrupt
#'   \eqn{\omega_k^{*(t)}} rather than fail loudly, so rather than trust a
#'   length mismatch to error, detect the shape against \code{beta} -- which is
#'   topics-by-words in every version -- and transpose an old object on read.
#' @param object a \code{tidylda} object
#' @return \code{object$counts$Cv} as a words-by-topics matrix.
counts_cv <- function(object) {
  cv <- object$counts$Cv

  if (is.null(cv)) {
    return(NULL)
  }

  k <- nrow(object$beta)

  # Words by topics already: ncol matches the topic count.
  if (ncol(cv) == k && nrow(cv) != k) {
    return(cv)
  }

  # Old topics-by-words layout.
  if (nrow(cv) == k && ncol(cv) != k) {
    return(Matrix::t(cv))
  }

  # Square, so shape cannot disambiguate: k == nrow(cv) == ncol(cv). Fall back on
  # the vocabulary labels instead.
  #
  # This is reachable only when the topic count exactly equals the vocabulary
  # size, which is near-impossible in real use and easy to hit in a toy example.
  # It was also broken until 0.1.0: Cv carried no dimnames, so the identity test
  # below could never succeed and every square Cv was transposed, corrupting a
  # correctly-oriented model rather than repairing a stale one. Models fitted by
  # 0.1.0 and later are labeled (see new_tidylda), so they take the first
  # branch; anything unlabeled is necessarily pre-0.1.0 and needs the
  # transpose.
  if (!is.null(colnames(object$beta)) && !is.null(rownames(cv)) &&
      identical(rownames(cv), colnames(object$beta))) {
    return(cv)
  }

  # Matrix::t() rather than t(): a 0.1.0 model stores Cv as a dgCMatrix, and
  # tidylda imports Matrix rather than attaching it, so a bare t() would dispatch
  # to t.default and error. It is correct for a base matrix too.
  Matrix::t(cv)
}

#' Construct a new object of class \code{tidylda}
#' @keywords internal
#' @description
#'   Since all three of \code{\link[tidylda]{tidylda}},
#'   \code{\link[tidylda]{refit.tidylda}}, and
#'   \code{\link[tidylda]{predict.tidylda}} call \code{\link[tidylda]{fit_lda_warp}},
#'   we need a way to format the resulting posteriors and other user-facing
#'   objects consistently. This function does that.
#' @param lda list output of \code{\link[tidylda]{fit_lda_warp}}
#' @param dtm a document term matrix or term co-occurrence matrix of class \code{dgCMatrix}
#' @param burnin integer number of burnin iterations.
#' @param is_prediction is this for a prediction (as opposed to initial fitting,
#'   or update)? Defaults to \code{FALSE}
#' @param alpha output of \code{\link[tidylda]{format_alpha}}
#' @param eta output of \code{\link[tidylda]{format_eta}}
#' @param optimize_alpha deprecated and ignored, retained so that callers
#'   passing it keep working. If \code{is_prediction = TRUE}, this argument is
#'   ignored.
#' @param calc_r2 did the user want to calculate R-squared when calculating the
#'   the model? If \code{is_prediction = TRUE}, this argument is ignored.
#' @param calc_likelihood did you calculate the log likelihood when making a call
#'   to \code{\link[tidylda]{fit_lda_warp}}?  If \code{is_prediction = TRUE}, this
#'   argument is ignored.
#' @param call the result of calling \code{\link[base]{match.call}} at the top of
#'   \code{\link[tidylda]{tidylda}}.
#' @param threads number of parallel threads
#' @return
#'   Returns an S3 object of class \code{tidylda} with the following slots:
#'
#'   \code{beta} is a numeric matrix whose rows are the posterior estimates
#'     of P(token|topic)
#'
#'   \code{theta} is a numeric matrix  whose rows are the posterior estimates of
#'     P(topic|document)
#'
#'   \code{lambda} is a numeric matrix whose rows are the posterior estimates of
#'     P(topic|token), calculated using Bayes's rule.
#'     See \code{\link[tidylda]{calc_lambda}}.
#'
#'   \code{alpha} is the prior for topics over documents. It is what the user
#'     passed when calling \code{\link[tidylda]{tidylda}}, formatted as a
#'     \code{k}-length numeric vector; the sampler does not modify it.
#'
#'   \code{eta} is the prior for tokens over topics. This is what the user passed
#'     when calling \code{\link[tidylda]{tidylda}}: a numeric scalar stays a
#'     scalar, and a matrix prior is a \code{k} by \code{ncol(dtm)} matrix.
#'
#'   \code{counts} is a list of two matrices holding the token-topic counts the
#'     sampler ended on. \code{Cd} is a dense matrix of documents by topics.
#'     \code{Cv} is a sparse \code{\link[Matrix]{dgCMatrix-class}} of tokens by
#'     topics, labeled with the model's vocabulary and topic names. Both are
#'     topics-in-columns, so \code{Cd} aligns with \code{theta} and \code{Cv}
#'     with \code{t(beta)}. If burn-in iterations were used these are averages
#'     over the post-burn-in iterations, and are therefore not integers.
#'
#'     \code{Cd} is deliberately dense: it is 38-81 percent nonzero, where a
#'     sparse form saves nothing and can cost 20 percent. \code{Cv} is 8-23
#'     percent nonzero and roughly 3 times smaller sparse.
#'
#'     NOTE: as of version 0.1.0 \code{Cv} is tokens by topics and sparse; it was
#'     topics by tokens and dense previously.
#'
#'   \code{summary} is the result of a call to \code{\link[tidylda]{summarize_topics}}
#'
#'   \code{call} is the result of \code{\link[base]{match.call}} called at the top
#'     of \code{\link[tidylda]{tidylda}}
#'
#'   \code{log_likelihood} is a \code{\link[tibble]{tibble}} with three columns,
#'     evaluated every \code{likelihood_every}-th iteration. \code{iteration} is
#'     the iteration number. \code{log_likelihood} is
#'     \eqn{P(tokens | \theta, \beta)}, the plug-in likelihood of the data under
#'     the current parameter estimates. \code{log_joint} is
#'     \eqn{P(tokens, topics | \alpha, \eta)}, the collapsed joint, with
#'     \code{theta} and \code{beta} integrated out. See
#'     \code{\link[tidylda]{tidylda}} for which to use when. This slot is only
#'     populated if \code{calc_likelihood = TRUE}
#'
#'   \code{r2} is a numeric scalar resulting from a call to
#'     \code{\link[mvrsquared]{calc_rsquared}}. This slot only populated if
#'     \code{calc_r2 = TRUE}
#' @note
#'   In general, the arguments of this function should be what the user passed
#'   when calling \code{\link[tidylda]{tidylda}}.
#'
#'   \code{burnin} is used only to determine whether or not burn in iterations
#'   were used when fitting the model. If \code{burnin > -1} then posteriors
#'   are calculated using \code{lda$Cd_mean} and \code{lda$Cv_mean} respectively.
#'   Otherwise, posteriors are calculated using \code{lda$Cd_mean} and
#'   \code{lda$Cv_mean}.
#'
#'   The class of \code{call} isn't checked. It's just passed through to the
#'   object returned by this function. Might be useful if you are using this
#'   function for troubleshooting or something.
new_tidylda <- function(
  lda, 
  dtm, 
  burnin, 
  is_prediction = FALSE,
  alpha = NULL, 
  eta = NULL,
  optimize_alpha = NULL, 
  calc_r2 = NULL,
  calc_likelihood = NULL, 
  call = NULL,
  threads
) {
  

  ### format theta ###
  # Cd is documents by topics and alpha is a K-vector, so alpha must be added
  # down the topic axis. Transposing first is what makes R's column-major
  # recycling land alpha[k] on topic k; adding it to the D x K matrix directly
  # recycles alpha diagonally instead, which silently corrupts theta for any
  # asymmetric alpha and is invisible for a symmetric one.
  if (burnin > -1) {
    Cd <- lda$Cd_mean
  } else {
    Cd <- lda$Cd
  }

  theta <- t(t(Cd) + lda$alpha)
  
  theta <- theta / rowSums(theta)
  
  # Guarded: is.na() on a D x K matrix allocates a D x K logical every call, to
  # repair a condition that essentially never arises. anyNA() short-circuits on
  # the first NA and scans nothing extra when there are none.
  if (anyNA(theta)) theta[is.na(theta)] <- 0 # just in case of a numeric issue
  
  colnames(theta) <- seq_len(ncol(theta))
  
  rownames(theta) <- rownames(dtm)
  
  ### format beta and all the rest ###
  
  if (!is_prediction) {
    ### format posteriors correctly ###
    # The engine returns Cv words-by-topics (D17). `beta` stays topics-by-words
    # for the public API, so the transpose happens here, once, rather than in
    # the engine on every run.
    if (burnin > -1) { # if you used burnin iterations use Cd_mean etc.
      Cv <- lda$Cv_mean
    } else { # if you didn't use burnin use standard counts (Cd etc.)
      Cv <- lda$Cv
    }

    # NAME THE AXES. Two reasons
    #
    # The obvious one: `beta` and `theta` are labeled, so `counts` should be
    # too. A user indexing Cv by token had no way to do it by name.
    #
    # The other: it is what makes a 0.1.0-or-later model self-identifying to
    # counts_cv(). That helper distinguishes the current words-by-topics layout
    # from the pre-0.1.0 topics-by-words one by shape, which cannot work when
    # k == ncol(dtm). Its tiebreaker compares rownames(Cv) against
    # colnames(beta) -- and before this, Cv carried no dimnames at all, so the
    # comparison could never succeed and the square case always transposed.
    #
    # The vocabulary here is the MODEL's, which is not always the data's: after
    # refit() adds vocabulary, ncol(dtm) < nrow(Cv). colnames(beta) is assigned
    # from colnames(dtm) further down for the same reason, and both are correct
    # because `dtm` at this point is the aligned matrix the sampler actually saw.
    rownames(Cv) <- colnames(dtm)
    colnames(Cv) <- colnames(theta)


    # A scalar prior arrives back as a 1 x 1 matrix (D20), so test length rather
    # than class; recycling then handles it.
    beta <- t(Cv) + if (length(lda$eta) == 1) as.numeric(lda$eta) else lda$eta
    
    beta <- beta / rowSums(beta)
    
    # Guarded for the same reason as theta above.
    if (anyNA(beta)) beta[is.na(beta)] <- 0 # just in case of a numeric issue
    
    colnames(beta) <- colnames(dtm)
    
    rownames(beta) <- colnames(theta)
    
    
    ### collect the results ###
    
    # lambda
    lambda <- calc_lambda(
      beta = beta, theta = theta,
      p_docs = Matrix::rowSums(dtm)
    )
    
    # eta
    if (eta$eta_class == "scalar") {
      eta_out <- as.numeric(lda$eta)[1]
    } else if (eta$eta_class == "vector") {
      eta_out <- lda$eta[1, ]
    } else if (eta$eta_class == "matrix") {
      colnames(lda$eta) <- colnames(beta)
      eta_out <- lda$eta
    } else { # this should be impossible, but science is hard and I am dumb.
      eta_out <- lda$eta
      
      warning("something went wrong formatting eta. refit.tidylda and predict.tidylda might be affected")
    }
    
    # alpha
    
    if (alpha$alpha_class == "scalar" & !optimize_alpha) {
      alpha_out <- lda$alpha[1]
    } else if (alpha$alpha_class == "vector" | optimize_alpha) {
      alpha_out <- lda$alpha
      
      names(alpha_out) <- rownames(beta)
    } else { # this should be impossible, but science is hard and I am dumb.
      alpha_out <- lda$alpha
      
      warning("something went wrong formatting alpha. refit.tidylda and predict.tidylda might be affected")
    }
    
    # resulting object
    summary <- try(
      tryCatch(
      summarize_topics(beta = beta, theta = theta, dtm = dtm),
      error = function(err){
        err$message <- "summarize_topics failed. model$summary corrupted."
        stop(err)
      })
    )

    # Two quantities, and they answer different questions. `log_likelihood` is
    # the plug-in P(data | parameters); `log_joint` is the collapsed joint
    # P(data, assignments | priors), with theta and beta integrated out. See
    # `?tidylda` for which to use when.
    log_likelihood <- as_tibble(data.frame(
      iteration = lda$log_likelihood[1, ],
      log_likelihood = lda$log_likelihood[2, ],
      log_joint = lda$log_likelihood[3, ]
    ))
    
    result <- list(
      beta = beta,
      theta = theta,
      lambda = lambda,
      alpha = alpha_out,
      eta = eta_out,
      summary = summary,
      call = call,
      log_likelihood = log_likelihood,
      counts = list(
        Cd = Cd,
        # SPARSE, AND ONLY THIS ONE. Cv is the largest slot in a fitted model at
        # scale --- V x K --- and it is genuinely sparse, because a word's tokens
        # can occupy at most min(n_w, K) topics and most words are rare. Measured
        # on nih_sample_dtm at k = 20 over 200 iterations: 7.7% nonzero without
        # burn-in, 23.2% with it.
        #
        # dgCMatrix carries a double value plus an int index per nonzero, so it
        # beats a dense matrix below 33% density for the integer counts and 67%
        # for the post-burn-in means. Cv clears both comfortably, shrinking by
        # roughly 4x and 3x.
        #
        # Cd stays DENSE on the same measurement: 38.4% and 81.3% nonzero, so
        # sparse storage would save nothing in the first case and cost 20% in the
        # second. D17 asked for both; only one of them is a good idea, and that
        # took measuring rather than reasoning.
        #
        # Converted HERE, at the point of storage, and not earlier: beta is built
        # from t(Cv) just above, and tidylda imports Matrix rather than attaching
        # it, so a bare t() on an S4 Matrix falls through to t.default and errors.
        # Keeping the conversion last means no internal arithmetic ever sees the
        # sparse form.
        Cv = Matrix::Matrix(Cv, sparse = TRUE)
      )
    )
    
    class(result) <- "tidylda"
    
    ### calculate and add other things ###
    
    # goodness of fit
    if (calc_r2) {
      result$r2 <- try(
        tryCatch(
          calc_lda_r2(
            dtm = dtm,
            theta = theta,
            beta = beta,
            threads
          ),
          error = function(err){
            err$message <- "calc_r2 failed. R-squared corrupted."
            stop(err)
          } 
        )
      )
    }
    
    # a little cleanup here
    if (!calc_likelihood) {
      result$log_likelihood <- NULL
    }
  }
  
  
  ### return the final result ###
  if (is_prediction) {
    return(theta)
  } else {
    return(result)
  }
}

#' Calculate R-squared for a tidylda Model
#' @keywords internal
#' @description Formats inputs and hands off to \link[mvrsquared]{calc_rsquared}
#' @param dtm must be of class dgCMatrix
#' @param theta a theta matrix
#' @param beta a beta matrix
#' @param threads number of parallel threads
#' @return Numeric scalar between negative infinity and 1
calc_lda_r2 <- function(dtm, theta, beta, threads) {
  
  # weight rows of theta by document length
  x <- Matrix::rowSums(dtm) * theta
  
  # calculate r-squared
  r2 <- mvrsquared::calc_rsquared(
    y = dtm,
    yhat = list(x = x, w = beta),
    return_ss_only = FALSE,
    threads = threads
  )
  
  # had an issue with preserving names
  names(r2) <- NULL
  
  r2
}

#' Utility function to tidy a simple triplet matrix
#' @keywords internal
#' @param x Object with rownames and colnames
#' @param triplets A data frame or list of i, j, x
#' @param row_names rownames, if not gotten from rownames(x)
#' @param col_names colnames, if not gotten from colnames(x)
#' @return returns a triplet matrix in the form of a data frame. The first
#'   column indexes rows. The second column indexes columns. The third column
#'   contains the i,j values.
#' @note This function ported from \code{\link[tidytext]{tidytext}}, copyright
#'   2017 David Robinson and Julia Silge. Moved the function here for stability
#'   reasons, as it is internal to tidytext
tidy_triplet <- function(x, triplets, row_names = NULL, col_names = NULL) {
  row <- triplets$i
  if (!is.null(row_names)) {
    row <- row_names[row]
  } else if (!is.null(rownames(x))) {
    row <- rownames(x)[row]
  }
  col <- triplets$j
  if (!is.null(col_names)) {
    col <- col_names[col]
  } else if (!is.null(colnames(x))) {
    col <- colnames(x)[col]
  }
  
  ret <- tibble::tibble(row = row, column = col, value = triplets$x)
  ret
}

#' Create a tidy tibble for a dgCMatrix
#' @keywords internal
#' @description Create a tidy tibble for a dgCMatrix. Will probably be a PR to
#'   \link[tidytext]{tidytext} in the future
#' @param x must be of class dgCMatrix
#' @param ... Extra arguments, not used
#' @return Returns a triplet matrix with columns "document", "term", and "count"
tidy_dgcmatrix <- function(x, ...) {
  triplets <- Matrix::summary(methods::as(x, "dgTMatrix"))
  ret <- tidy_triplet(x, triplets)
  colnames(ret) <- c("document", "term", "count")
  ret
}

#' Calculate a matrix whose rows represent P(topic_i|tokens)
#' @keywords internal
#' @description 
#' Use Bayes' rule to get P(topic|token) from the estimated parameters of a
#' probabilistic topic model.This resulting "lambda" matrix can be used for
#' classifying new documents in a frequentist context and supports
#' \code{\link[tidylda]{augment}}.
#' @param theta a theta matrix
#' @param beta a beta matrix
#' @param p_docs A numeric vector of length \code{nrow(theta)} that is
#'   proportional to the number of terms in each document,  defaults to NULL.
#' @param correct Logical. Do you want to set NAs or NaNs in the final result to
#'   zero? Useful when hitting computational underflow. Defaults to \code{TRUE}.
#'   Set to \code{FALSE} for troubleshooting or diagnostics.
#' @return
#' Returns a \code{matrix} whose rows correspond to topics and whose columns
#' correspond to tokens. The i,j entry corresponds to P(topic_i|token_j)
calc_lambda <- function(beta, theta, p_docs = NULL, correct = TRUE){
  
  # set up constants
  D <- nrow(theta)
  K <- ncol(theta)
  V <- ncol(beta)
  
  # probability of each document (assumed to be equiprobable)
  if(is.null(p_docs)){
    p_d <- rep(1/nrow(theta), nrow(theta))
  }else{
    if(sum(is.na(p_docs)) > 0){
      warning("found missing values in p_docs. Setting them as 0.")
      p_docs[ is.na(p_docs) ] <- 0 
    }
    p_d <- p_docs / sum(p_docs)
  }
  
  # get the probability of each topic
  p_t <- p_d %*% theta
  
  # get the probability of each word from the model    
  p_w <- p_t %*% beta
  
  
  
  # get our result
  #
  # This used to build a k by k matrix, set its diagonal to p_t, and multiply it
  # into beta -- an O(k^2 V) matmul whose off-diagonal terms are all zero, to do
  # what row scaling does in O(kV). beta is k by V, so a length-k vector recycles
  # down the columns and scales row i by p_t[i], which is what the diagonal
  # matmul computed. p_t arrives as a 1 by k matrix from `p_d %*% theta`, hence
  # as.numeric().
  lambda <- as.numeric(p_t) * beta
  
  # And this used to be t(apply(lambda, 1, function(x) x / p_w)), which builds a
  # V by k result and transposes it: two k by V copies to do one elementwise
  # division. p_w is indexed by token, so the division is along dimension 2.
  lambda <- sweep(lambda, 2, as.numeric(p_w), "/")
  
  rownames(lambda) <- rownames(beta)
  colnames(lambda) <- colnames(beta)
  
  # give us zeros instead of NAs when we have NA or NaN entries
  if (correct) {
    lambda[is.na(lambda)] <- 0 
  }
  
  return(lambda)
}

#' Probabilistic coherence of topics
#' @description Calculates the probabilistic coherence of a topic or topics. 
#' This approximates semantic coherence or human understandability of a topic.
#' @param beta A numeric matrix or a numeric vector. The vector, or rows of the 
#' matrix represent the numeric relationship between topic(s) and terms. For
#' example, this relationship may be p(word|topic) or p(topic|word).
#' @param data A document term matrix or term co-occurrence matrix. The preferred
#'   class is a \code{\link[Matrix]{dgCMatrix-class}}. However there is support
#'   for any \code{\link[Matrix]{Matrix-class}} object as well as several other
#'   commonly-used classes such as \code{\link[base]{matrix}},
#'   \code{\link[quanteda]{dfm}}, \code{\link[tm]{DocumentTermMatrix}}, and
#'   \code{\link[slam]{simple_triplet_matrix}}
#' @param m An integer for the number of words to be used in the calculation. 
#' Defaults to 5
#' @return Returns an object of class \code{numeric} corresponding to the 
#' probabilistic coherence of the input topic(s).
#' @details 
#'   For each pair of words \{a, b\} in the top M words in a topic, probabilistic
#'   coherence calculates P(b|a) - P(b), where \{a\} is more probable than \{b\} in
#'   the topic. For example, suppose the top 4 words in a topic are \{a, b, c, d\}.
#'   Then, we calculate 1. P(a|b) - P(b), P(a|c) - P(c), P(a|d) - P(d)
#'   2. P(b|c) - P(c), P(b|d) - P(d)
#'   3. P(c|d) - P(d)
#'   All 6 differences are averaged together.
#' @examples
#' # Load a pre-formatted dtm and topic model
#' data(nih_sample_dtm)
#' 
#' # fit a model
#' set.seed(12345)
#' model <- tidylda(
#'   data = nih_sample_dtm[1:20, ], k = 5,
#'   iterations = 100, burnin = 50
#' )
#' 
#' calc_prob_coherence(beta = model$beta, data = nih_sample_dtm, m = 5)
#' @export 
calc_prob_coherence <- function(beta, data, m = 5){
  # code below is ported almost verbatim from textmineR. Copied here to reduce
  # cross dependencies between textmineR and tidylda
  
  # beta is a numeric matrix or numeric vector?
  if( ! is.numeric(beta) ){
    stop("beta must be a numeric matrix whose rows index topics and columns\n",
         " index terms or beta must be a numeric vector whose entries index terms.")
  }
  # Ensure dtm is of class dgCMatrix
  dtm <- convert_dtm(dtm = data)
  
  # is m numeric? If it is not an integer, give a warning.
  if( ! is.numeric(m) | m < 1){
    stop("M must be an integer in 1:ncol(beta) or 1:length(beta)")
  }
  
  if(length(m) != 1){
    warning("m is a vector when scalar is expected. Taking only the first value")
    m <- m[ 1 ]
  }
  
  if(floor(m) != m){
    warning("m is expected to be an integer. floor(m) is being used.")
    m <- floor(m)
  }
  
  # # dtm has colnames?
  # if( is.null(colnames(dtm))){
  #   stop("dtm must have colnames")
  # }
  
  # Names of beta in colnames(dtm)
  if( ! is.matrix(beta) ){
    if(sum(names(beta)[ 1:m ] %in% colnames(dtm)) != length(1:m)){
      stop("vocabulary of beta (i.e., colnames(beta)) does not match vocabulary of data")
    }
  }else if(sum(colnames(beta)[ 1:m ] %in% colnames(dtm)) != length(1:m)){
    stop("vocabulary of beta (i.e., colnames(beta)) does not match vocabulary of data")
  }
  
  # Pick the top m terms of one topic, as COLUMN INDICES into the dtm.
  #
  # PARTIAL SORT. order() sorts all V entries to keep m, which is O(V log V)
  # where selection is O(V). Same reasoning as summarize_topics()'s top_terms,
  # and the same idiom.
  #
  # TIES ARE BROKEN BY TERM INDEX, which top_terms does not bother to do. There
  # it does not matter --- those are display terms. Here the selected terms feed
  # the returned coherence, so a different order among tied values would change
  # a number the user sees. Sorting the candidate set by (-value, index) makes
  # the choice independent of how the candidates were found, which is what keeps
  # this identical to the old order()-based selection on every input rather than
  # merely on inputs without ties.
  top_idx <- function(x, m) {
    if (length(x) <= m) {
      return(order(x, decreasing = TRUE)[seq_len(m)])
    }
    cut <- length(x) - m + 1L
    keep <- which(x >= sort(x, partial = cut)[cut])
    keep[order(-x[keep], keep)][seq_len(m)]
  }

  # Declare a function to get probabilistic coherence on one topic, given the
  # m x m matrix of co-occurrence probabilities for that topic's terms.
  pcoh <- function(p.mat){
    p.diag <- diag(p.mat)
    result <- sapply(1:(ncol(p.mat) - 1), function(x) {
      p.mat[x, (x + 1):ncol(p.mat)]/p.mat[x, x] -
        p.diag[(x + 1):ncol(p.mat)]
    })
    mean(unlist(result), na.rm = TRUE)
  }

  # ONE CROSSPRODUCT FOR THE WHOLE MODEL, not one per topic.
  #
  # This used to subset the dtm to a topic's m terms, binarize that subset, and
  # crossprod it --- k sparse column-subsets and k crossproducts to build what
  # is really one small co-occurrence problem. Topics share terms heavily, so
  # the union of all top-m sets is far smaller than k*m (236 distinct terms at
  # k = 100, 918 at k = 600 on 20 Newsgroups).
  #
  # So: select every topic's terms, take the union, binarize that submatrix
  # once, and take a single crossproduct. Each topic then reads its m x m block
  # out of the result. Measured 3.7x faster at k = 100 and 2.2x at k = 600, with
  # bit-identical output at k = 100, 300 and 600.
  #
  # Binarizing by overwriting @x is what makes it one pass: dtm[dtm > 0] <- 1
  # builds a logical sparse matrix and assigns through it, where a dgCMatrix's
  # stored values are exactly its nonzeros by construction.
  is_vec <- !is.matrix(beta)

  sel <- if (is_vec) {
    matrix(top_idx(beta, m), nrow = 1)
  } else {
    t(apply(beta, 1, top_idx, m = m))
  }

  terms_used <- sort(unique(as.vector(sel)))

  dtm_bin <- dtm[, terms_used, drop = FALSE]
  dtm_bin@x <- rep(1, length(dtm_bin@x))

  p_all <- Matrix::crossprod(dtm_bin) / nrow(dtm)

  # Positions of each topic's terms within the union, so the per-topic block is
  # a plain integer index rather than another name lookup.
  pos <- matrix(match(as.vector(sel), terms_used), nrow = nrow(sel))

  out <- vapply(
    seq_len(nrow(pos)),
    function(i) pcoh(as.matrix(p_all[pos[i, ], pos[i, ], drop = FALSE])),
    numeric(1)
  )

  if (is_vec) {
    return(out[[1]])
  }

  names(out) <- rownames(beta)

  out
}

