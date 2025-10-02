#' Build the full ME + CME model matrix
#'
#' @description
#' Constructs the full design matrix containing main effects (MEs) and
#' conditional main effects (CMEs) from an input ME matrix. For each ordered
#' pair \eqn{(i, j)}, \eqn{i \ne j}, two CME columns are created:
#' \deqn{X_{i|j+} = X_i \cdot \mathbf{1}\{X_j \ge 0\}, \qquad
#'       X_{i|j-} = X_i \cdot \mathbf{1}\{X_j < 0\}.}
#' Columns are named as \code{"<MEi>|<MEj>+"} and \code{"<MEi>|<MEj>-"}.
#'
#' @param xme Numeric matrix of main effects with \eqn{n} rows and \eqn{p} columns.
#'   If \code{colnames(xme)} are missing, they are set to \code{V1, V2, ..., Vp}.
#'
#' @details
#' The output matrix has the first \eqn{p} columns equal to \code{xme}, followed by
#' \eqn{2\,p(p-1)} CME columns ordered by nested loops over \eqn{i=1,\dots,p} and
#' \eqn{j=1,\dots,p} with \eqn{i \ne j}, and within each pair the \code{"+", "-"}
#' columns in that order.
#'
#' The thresholding for CMEs uses \code{X_j >= 0} vs \code{X_j < 0}. In typical CME
#' workflows, \code{xme} is sign-coded (e.g., \code{-1/+1}) or centered so that
#' zero is a meaningful split.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{model.mtx} — numeric \eqn{n \times (p + 2\,p(p-1))} matrix combining
#'         MEs (first \eqn{p} columns) and CMEs (remaining columns). Column names are
#'         \code{xme} names for MEs and \code{"<MEi>|<MEj>+"}, \code{"<MEi>|<MEj>-"} for CMEs.
#'   \item \code{cme.mtx} — integer matrix of size \eqn{2\,p(p-1) \times 2}; each row
#'         stores the \code{(i, j)} indices (1-based) corresponding to the CME column at
#'         the same position in \code{model.mtx[, (p+1):ncol(model.mtx)]}. Rows appear in
#'         the same order as the CME columns (for each \eqn{(i,j)}: \code{+} then \code{-}).
#' }
#'
#' @references
#' Mak, S., & Wu, C. J. (2019). cmenet: A new method for bi-level variable selection
#' of conditional main effects. *Journal of the American Statistical Association*, 114(526), 844-856.
#'
#' @source
#' Function structure and notation follow the \pkg{cmenet} R package.
#'
#'
#' @seealso \code{\link{glmcmenet}}, \code{\link{predictcme}}
#' @export


full.model.mtx <- function (xme)
{
  memtx <- xme
  nn <- nrow(memtx)
  pp <- ncol(memtx)
  colnret <- rep(NA, pp + 2 * pp * (pp - 1))
  ret <- matrix(NA, nrow = nn, ncol = 2 * pp * (pp - 1) +
                  pp)
  ininames <- rep(NA, pp)
  if (is.null(colnames(memtx))) {
    for (i in 1:pp) {
      ininames[i] <- paste0("V", i)
    }
    colnames(memtx) <- ininames
  }
  ret[, 1:pp] <- memtx
  colnret[1:pp] <- colnames(memtx)
  cmenames <- rep(NA, 2 * pp * (pp - 1))
  count <- pp + 1
  intmtx <- matrix(NA, nrow = 2 * pp * (pp - 1), ncol = 2)
  for (i in 1:pp) {
    for (j in 1:pp) {
      if (i != j) {
        ret[, count] <- memtx[, i] * as.numeric(memtx[,
                                                      j] >= 0)
        cmenames[count - pp] <- paste0(colnames(memtx)[i],
                                       "|", colnames(memtx)[j], "+")
        intmtx[count - pp, 1] <- i
        intmtx[count - pp, 2] <- j
        count <- count + 1
        ret[, count] <- memtx[, i] * as.numeric(memtx[,
                                                      j] < 0)
        cmenames[count - pp] <- paste0(colnames(memtx)[i],
                                       "|", colnames(memtx)[j], "-")
        intmtx[count - pp, 1] <- i
        intmtx[count - pp, 2] <- j
        count <- count + 1
      }
    }
  }
  colnret[(pp + 1):ncol(ret)] <- cmenames
  colnames(ret) <- colnret
  return(list(model.mtx = ret, cme.mtx = intmtx))
}
