#' Compute CV loss over a tuning grid
#'
#' @description
#' Evaluates a loss/score for cross-validation across the full tuning grid
#' (columns × lambdas), using predictions \code{yhat} and a fitted CME object
#' to obtain model size (nonzeros). Supported measures are
#' \code{"deviance"}, \code{"class"}, and \code{"bic"}.
#'
#' @param fitobj A fitted CME model (e.g., from \code{\link{glmcmenet}}) that
#'   contains a 3D array \code{$coefficients} of shape
#'   \code{(p_total) × n_lambda_sib × n_lambda_cou}. The function uses this to
#'   count nonzeros per grid point.
#' @param y Numeric vector of length \code{n} with observed responses.
#'   For \code{family == "binomial"}, values should be 0/1.
#' @param yhat Array of predictions for all grid points; expected shape
#'   \code{n × n_lambda_sib × n_lambda_cou}. Typically produced by
#'   \code{\link{predictcme}} with \code{type = "response"} (for deviance/BIC)
#'   or \code{type = "link"} (if used differently).
#' @param n Integer number of observations (usually \code{length(y)}).
#' @param family Character; one of \code{"gaussian"}, \code{"binomial"}, or \code{"poisson"}.
#' @param type.measure Character; loss to compute. One of
#'   \code{"deviance"}, \code{"class"}, or \code{"bic"}.
#'
#' @details
#' The number of selected coefficients at each grid point is computed as
#' \code{selmat[i,j] = sum(beta[,i,j] != 0)} from \code{fitobj$coefficients}.
#'
#' \strong{Deviance:}
#' \itemize{
#'   \item Gaussian: \eqn{(y - \hat y)^2}.
#'   \item Binomial: \eqn{-2\,[y\log(\hat y) + (1-y)\log(1-\hat y)]}.
#'   \item Poisson: \eqn{2\,[y\log y - y + \hat y - y\log \hat y]}, with
#'         the convention \eqn{y\log y = 0} when \eqn{y = 0}.
#' }
#'
#' \strong{Classification error (binomial only):}
#' returns a logical array \code{(yhat < 0.5) == y}. For non-binomial families,
#' this falls back to squared error.
#'
#' \strong{BIC-like score:}
#' adds a penalty \eqn{\text{selmat}\cdot\log(n)/n} broadcast over samples to the
#' deviance for each family. Any \code{-Inf} created by logs is replaced by a large
#' finite value.
#'
#' @return
#' An array with the same dimensions as \code{yhat}
#' (\code{n × n_lambda_sib × n_lambda_cou}) containing the pointwise loss for
#' each observation and grid point. (The caller usually averages over the first
#' dimension.)
#'
#'
#' @seealso \code{\link{predictcme}}, \code{\link{cv.glmcmenet}}
#' @export


loss <- function(fitobj,y,yhat,n,family,type.measure=c("deviance","class","bic")){

  selmat <- apply(fitobj$coefficients,c(2,3),function(x){return(length(which(x!=0)))})

  if(type.measure=="deviance"){
    val <- array(NA, dim = dim(yhat))
    if (family=="gaussian") {
      val <- (y-yhat)^2
    } else if (family=="binomial") {
      val <- -2*(y*log(yhat)+(1-y)*log(1-yhat))
    }  else if (family=="poisson") {
      yly <- y*log(y)
      yly[y==0] <- 0
      val <- 2*(yly - y + yhat - y*log(yhat))
    }
  }else if(type.measure=="class"){
    val <- array(NA, dim = dim(yhat))
    if (family=="binomial") {
      val <- (yhat < 0.5) == y
    } else {
      val <- (y-yhat)^2
    }
  }else if(type.measure=="bic"){
    val <- array(NA, dim = dim(yhat))
    if (family=="gaussian") {
      val <- (y-yhat)^2+
        aperm(array(rep(selmat*log(n)/n, times=dim(yhat)[1]), dim(yhat)[c(2,3,1)]),c(3,1,2))
    } else if (family=="binomial") {
      val <- -2*(y*log(yhat)+(1-y)*log(1-yhat))+
        aperm(array(rep(selmat*log(n)/n, times=dim(yhat)[1]), dim(yhat)[c(2,3,1)]),c(3,1,2))
    }  else if (family=="poisson") {
      yly <- y*log(y)
      yly[y==0] <- 0
      val <- 2*(yly - y + yhat - y*log(yhat))+
        aperm(array(rep(selmat*log(n)/n, times=dim(yhat)[1]), dim(yhat)[c(2,3,1)]),c(3,1,2))
    }
    val <- replace(val, val == -Inf, 9999999)
  }
}
