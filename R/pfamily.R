#' Map linear predictor to mean on the response scale
#'
#' @description
#' Internal helper that converts a linear predictor \code{eta} to the mean
#' \eqn{\mu} under a GLM family:
#' identity (gaussian), logistic (binomial, with clipping at \eqn{\pm 16} for
#' numerical stability), and exponential (poisson).
#'
#' @param eta Numeric scalar linear predictor.
#' @param family Character scalar; one of \code{"gaussian"}, \code{"binomial"},
#'   or \code{"poisson"}.
#'
#' @return Numeric scalar \eqn{\mu} on the response scale.
#'
#' @details
#' This function is **not vectorized**; use from higher-level code that applies
#' it element-wise (e.g., via \code{sapply}/\code{vapply}), as done in
#' \code{predictcme()}.
#'
#' @keywords internal
#' @noRd


pfamily <- function (eta,family){
  if(family=="gaussian"){
    return(eta);
  }else if(family=="binomial"){
      if (eta > 16) {
        return(0.9999);
      } else if (eta < -16) {
        return(0.0001);
      } else {
      return(exp(eta)/(1+exp(eta)));
      }
    }else if(family=="poisson"){
      return(exp(eta));
    }
  }


