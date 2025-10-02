#' glmcmenet: Conditional Main-Effect Regularization for GLMs
#'
#'@description
#' \code{glmcmenet} performs variable selection of conditional main effects (CMEs)
#' via a adaptive bi-level penalization framework for generalized linear models.
#'
#'@references
#' Xie, K., & Xinwei, Deng. (2025). Adaptive Bi-Level Variable Selection of Conditional
#' Main Effects for Generalized Linear Models. *Technometrics*, to be appear.
#' Mak, S., & Wu, C. J. (2019). cmenet: A new method for bi-level variable selection
#' of conditional main effects. *Journal of the American Statistical Association*, 114(526), 844-856.
#'
#' @section C++ integration:
#' The package uses Rcpp/RcppArmadillo. Native routines are registered.
#'
#' @useDynLib glmcmenet, .registration = TRUE
#' @keywords internal
"_PACKAGE"
