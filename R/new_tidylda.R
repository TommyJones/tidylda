#' Construct a new object of class \code{tidylda}
#' @keywords internal
#' @description
#'   This is a helper function for creating new \code{tidylda} objects.
#' @param phi a numeric matrix whose rows are the posterior estimates
#'   of P(token|topic)
#' @param theta a numeric matrix  whose rows are the posterior estimates of
#'   P(topic|document)
#' @param gamma a numeric matrix whose rows are the posterior estimates of
#'     P(topic|token)
#' @param alpha either a numeric scalar or numeric vector of length ncol(theta);
#'   the prior for topics over documents
#' @param beta either a numeric scalar, numeric vector of length ncol(phi), or
#'   numeric matrix with dimensions dim(phi); the prior for tokens over topics.
#' @param summary a \code{\link[tibble]{tibble}} with the following columns:
#'   \code{topic} (numeric), \code{prevalence} (numeric), \code{coherence} (numeric),
#'   \code{top_terms} (character)
#' @param call the function call used to create the model. See note, below.
#' @param ... Other elements you would like to populate in the \code{tidylda} object.
#' @return
#'   Returns an S3 object of class \code{tidylda} with the following slots:
#'
#'   \code{phi} is a numeric matrix whose rows are the posterior estimates
#'     of P(token|topic)
#'
#'   \code{theta} is a numeric matrix  whose rows are the posterior estimates of
#'     P(topic|document)
#'
#'   \code{gamma} is a numeric matrix whose rows are the posterior estimates of
#'     P(topic|token). This is calculated by making a call to
#'     \code{link[textmineR]{CalcGamma}} which uses Bayes's rule to calculate
#'     \code{gamma} from \code{phi}, \code{theta}, and P(document) (which is
#'     proportional to \code{Matrix::rowSums(dtm)}).
#'
#'   \code{alpha} is the prior for topics over documents. If \code{optimize_alpha}
#'     is \code{FALSE}, \code{alpha} is what the user passed when calling
#'     \code{\link[tidylda]{tidylda}}. If \code{optimize_alpha} is \code{TRUE},
#'     \code{alpha} is a numeric vector returned in the \code{alpha} slot from a
#'     call to \code{\link[tidylda]{fit_lda_c}}.
#'
#'   \code{beta} is the prior for tokens over topics. This is what the user passed
#'     when calling \code{\link[tidylda]{tidylda}}.
#'
#'   \code{summary} is the result of a call to \code{\link[tidylda]{summarize_topics}}
#'
#'   \code{call} is the result of \code{\link[base]{match.call}} called at the top
#'     of \code{\link[tidylda]{tidylda}}
#'
#'   \code{log_likelihood} is a \code{\link[tibble]{tibble}} whose columns are
#'     the iteration and log likelihood at that iteration. This slot is only populated
#'     if \code{calc_likelihood = TRUE}
#'
#'   \code{r2} is a numeric scalar resulting from a call to
#'     \code{\link[textmineR]{CalcTopicModelR2}}. This slot only populated if
#'     \code{calc_r2 = TRUE}
#' @note
#'   You can technically pass a \code{\link[base]{data.frame}} as the \code{summary}
#'   parameter. However, you may experience unintended side effects when using
#'   \code{\link[tidylda]{print.tidylda}}.
#'
#'   The class of \code{call} isn't checked. It's just passed through to the
#'   object returned by this function. Might be useful if you are using this
#'   function for troubleshooting or something.
new_tidylda <- function(phi, theta, gamma, alpha, beta, summary, call, ...) {

  ### formatting checks ----

  # class checks
  if (!inherits(phi, "matrix")) {
    stop("phi must be a numeric matrix")
  }

  if (!inherits(theta, "matrix")) {
    stop("theta must be a numeric matrix")
  }

  if (!inherits(gamma, "matrix")) {
    stop("gamma must be a numeric matrix")
  }

  if (!inherits(alpha, "numeric") & !inherits(alpha, "integer")) {
    stop("alpha must be a numeric scalar or numeric vector")
  }

  if (!inherits(beta, "numeric") & !inherits(beta, "integer") &
    !inherits(beta, "matrix")) {
    stop("beta must be a numeric scalar, numeric vector, or numeric matrix")
  }

  # compatibility of dimensionality
  if (ncol(theta) != nrow(phi)) {
    stop(
      "theta and phi matrices have incompatible dimensionality. ",
      "ncol(theta) must equal nrow(phi). I see ncol(theta) = ",
      ncol(theta), " and nrow(phi) = ", nrow(phi)
    )
  }

  if (nrow(phi) != nrow(gamma) | ncol(phi) != ncol(gamma)) {
    stop(
      "phi and gamma matrices must have the same dimensionality. ",
      "I see dim(phi) = c(", paste(dim(phi), collapse = ", "), ") ",
      "and dim(gamma) = c(", paste(dim(gamma), collapse = ", "), ")"
    )
  }

  if (length(alpha) > 1 & length(alpha) != ncol(theta)) {
    stop(
      "Since alpha is a vector, length(alpha) must equal ncol(theta). ",
      "I see length(alpha) = ", length(alpha), " But ncol(theta) = ",
      ncol(theta)
    )
  }

  if (is.matrix(beta)) {
    if (ncol(beta) != ncol(phi) | nrow(beta) != nrow(phi)) {
      stop(
        "Since beta is a matrix, it must have the same dimensionality as phi. ",
        "I see dim(beta) = c(", paste(dim(beta), collapse = ", "), ") but ",
        "dim(phi) = c(", paste(dim(phi), collapse = ", "), ")"
      )
    }
  }

  if (!is.matrix(beta) & length(beta) > 1 & length(beta) != ncol(phi)) {
    stop(
      "Since beta is a vector, length(beta) must equal ncol(phi). ",
      "I see length(beta) = ", length(beta), " But ncol(phi) = ",
      ncol(phi)
    )
  }

  # does summary conform to the minimum needed for the print method to work?
  if (!inherits(summary, "tbl") & !inherits(summary, "data.frame")) {
    stop(
      "class(summary) should be 'tbl' or 'data.frame'. However I see that ",
      "class(summary) = ", class(summary)
    )
  }

  if (sum(c("topic", "prevalence", "coherence", "top_terms") %in% names(summary)) < 4) {
    stop(
      "summary must have names c('topic', 'prevalence', 'coherence', 'top_terms'). ",
      "However, I see names(summary) = c('", paste(names(summary), collapse = "', '"),
      "')"
    )
  }


  ### construct and return ----
  result <- list(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = call,
    ...
  )

  class(result) <- "tidylda"

  result
}
