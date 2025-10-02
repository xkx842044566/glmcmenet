#' Indices of the minimum value (vector or matrix)
#'
#' @description
#' Returns the position of the smallest element of \code{x}. For a vector,
#' this is a single index. For a matrix, this is a length-2 integer vector
#' giving \code{(row, column)} of the first minimum (in column-major order).
#'
#' @param x A numeric vector or matrix.
#'
#' @details
#' Ties are broken by the first occurrence when \code{x} is vectorized in
#' column-major order (the same behavior as \code{order()} / \code{which.min()}).
#' For matrices, indices are 1-based \code{(row, col)}.
#'
#' @return
#' If \code{x} is a vector, a single integer index.
#'
#' If \code{x} is a matrix, an integer vector of length 2: \code{c(row, col)}.
#'
#' @seealso \code{\link[base]{which.min}}, \code{\link[base]{order}}
#' @keywords internal
#' @noRd


argmin <- function (x)
{
  vx = as.vector(x)
  imax = order(vx)[1]
  if (!is.matrix(x))
    imax
  else {
    d = dim(x)
    c1 = as.vector(outer(seq(d[1]), rep(1, d[2])))[imax]
    c2 = as.vector(outer(rep(1, d[1]), seq(d[2])))[imax]
    c(c1, c2)
  }
}
