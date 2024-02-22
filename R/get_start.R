dev_function <- function(y, mu, weights, family) {
  sum(family$dev.resids(y, mu,weights))
}

get_start <- function(x, y, family, intercept) {
  nobs <- nrow(x); nvars <- ncol(x)

  # compute mu and null deviance
  # family = binomial() gives us warnings due to non-integer weights
  # to avoid, suppress warnings
  weights = rep(1,nobs)

  mu <- rep(mean(y), times = nobs)

  nulldev <- dev_function(y, mu, weights/sum(weights), family)

  # if some penalty factors are zero, we have to recompute mu

  # compute lambda max
  ju <- rep(1, nvars)
  r <- y - mu
  eta <- family$linkfun(mu)
  v <- family$variance(mu)
  m.e <- family$mu.eta(eta)
  rv <- r / v * m.e
  if (inherits(x, "sparseMatrix")) {
    xm <- attr(x, "xm")
    xs <- attr(x, "xs")
    g <- abs((drop(t(rv) %*% x) - sum(rv) * xm) / xs)
  } else {
    g <- abs(drop(t(rv) %*% x))
  }
  g <- g * ju
  lambda_max <- max(g) / max(1, 1e-3)

  list(nulldev = nulldev, mu = mu, lambda_max = lambda_max)
}
