#' Starting values and lambda_max for the GLM path
#'
#' @description
#' Computes initial mean vector and the maximum penalty (`lambda_max`) used to
#' start regularization paths for Gaussian, Binomial, or Poisson models.
#'
#' @param x Numeric matrix (n x p) design matrix.
#' @param y Numeric response vector of length n (0/1 for binomial).
#' @param family Character string: "gaussian", "binomial", or "poisson".
#' @param intercept Logical; currently ignored.
#'
#' @return A list with:
#' \describe{
#'   \item{mu}{Length-n vector with all entries equal to `mean(y)`.}
#'   \item{lambda_max}{Scalar giving the maximal penalty.}
#' }
#'
#' @keywords internal
#' @noRd



get_start <- function(x, y, family, intercept) {
  nobs <- nrow(x); nvars <- ncol(x)

  # compute mu and null deviance
  # family = binomial() gives us warnings due to non-integer weights
  # to avoid, suppress warnings

  mu <- rep(mean(y), times = nobs)

  #nulldev <- dev_function(y, mu, family)

  ju <- rep(1, nvars)
  r <- y - mu

  # if some penalty factors are zero, we have to recompute mu

  # compute lambda max
  # Adjust calculations based on family
  if (family=="gaussian"){
    eta <- mu
    v <- 1
    m.e <- 1
  } else if (family == "binomial") {
    eta <- log(mu / (1 - mu)) # logit link for binomial
    v <- mu * (1 - mu) # binomial variance
    m.e <- mu * (1 - mu) # derivative of mu w.r.t. eta for binomial with logit link
  } else if (family == "poisson") {
    eta <- log(mu) # log link for poisson
    v <- mu # poisson variance
    m.e <- 1 # derivative of mu w.r.t. eta for poisson with log link
  }
  rv <- r / v * m.e
  if (inherits(x, "sparseMatrix")) {
    xm <- attr(x, "xm")
    xs <- attr(x, "xs")
    g <- abs((drop(t(rv) %*% x) - sum(rv) * xm) / xs)
  } else {
    g <- abs(drop(t(rv) %*% x))
  }
  g <- g * ju
  lambda_max <- max(g) / max(1, 1e-3) /nobs

  list(mu = mu, lambda_max = lambda_max) #nulldev = nulldev,
}
