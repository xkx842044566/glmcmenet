#' Cross-validated tuning for CME-regularized GLMs
#'
#' @description
#' Performs K-fold cross-validation to tune the CME model in two stages:
#' (1) a grid over \eqn{(\tau, \gamma)} then (2) a grid over \eqn{(\lambda_{\mathrm{sib}}, \lambda_{\mathrm{cou}})}.
#' Optionally uses warm starts from external solvers
#' (lasso / adaptive lasso / elastic net / MCP) to set an initial active set.
#'
#' @param xme Numeric matrix (n × p\_me) of main-effect columns.
#' @param xcme Numeric matrix (n × p\_cme) of CME columns.
#' @param y Numeric response vector of length \code{n}.
#' @param family Character; one of \code{"gaussian"}, \code{"binomial"}, or \code{"poisson"}.
#' @param nfolds Integer number of CV folds (default \code{10}). Must be \eqn{\ge} 3.
#' @param var.names Optional character vector of variable names (length \code{ncol(xme)+ncol(xcme)})
#'   used when reporting selected features.
#' @param nlambda.sib, nlambda.cou Integers: grid sizes for the sibling and cousin
#'   penalty sequences in stage (2) (defaults \code{20} each).
#' @param lambda.min.ratio Smallest \eqn{\lambda} as a fraction of the largest \eqn{\lambda}
#'   in each sequence (default \code{1e-6}).
#' @param ngamma,Integer Number of \code{gamma} values in stage (1) (default \code{20}).
#' @param max.gamma Largest \code{gamma} considered (default \code{150}).
#' @param ntau Integer number of \code{tau} values (default \code{20}).
#' @param max.tau Largest \code{tau} (default \code{0.2}).
#' @param tau.min.ratio Smallest \code{tau} as a fraction of \code{max.tau} (default \code{0.001}).
#' @param it.max,Integer Max iterations for the final refit on all data (default \code{250}).
#' @param it.max.cv Integer Max iterations used inside each CV fit (default \code{25}).
#' @param type.measure Character; CV loss to minimize. One of
#'   \code{"deviance"}, \code{"class"} (misclassification for binomial), or \code{"bic"}.
#' @param warm.str Character warm-start strategy for the initial active set, one of
#'   \code{"lasso"}, \code{"adaptive_lasso"}, \code{"elastic"}, \code{"ncvreg"}, or \code{"NULL"}.
#'   If not \code{"NULL"}, external solvers (glmnet/ncvreg) are used to preselect variables.
#' @param elastic_alpha Optional numeric \eqn{\alpha \in [0,1]} for the elastic-net warm start.
#' @param penalty.factor Numeric vector of per-coefficient weights (length \code{p\_me+p\_cme});
#'   modulates penalties (default all ones).
#' @param group.penalty Numeric vector of length \code{2*p\_me} with group weights for
#'   siblings and cousins (ordered \code{c(m_sib, m_cou)}). Default all ones.
#' @param screen_ind Logical; enable screening rules inside fits (default \code{FALSE}).
#'
#' @details
#' Stage (1) builds a grid for \code{tau} (from \code{max.tau} down to \code{max.tau * tau.min.ratio})
#' and a grid for \code{gamma} (from \code{max.gamma} down to a data-driven lower bound),
#' runs K-fold CV for each pair, and chooses the pair minimizing the mean CV loss.
#' Stage (2) then fixes \code{(tau, gamma)} at those values and cross-validates over
#' \code{lambda.sib × lambda.cou}.
#'
#' Warm starts (when \code{warm.str != "NULL"}) set \code{act.vec} by fitting an external model:
#' lasso / adaptive lasso / elastic net via \pkg{glmnet}, or MCP via \pkg{ncvreg}.
#'
#' @return
#' An object of class \code{"cv.glmcme"} (a list) with components:
#' \itemize{
#'   \item \code{params} — named vector with selected \code{lambda.sib}, \code{lambda.cou},
#'         \code{gamma}, and \code{tau}.
#'   \item \code{cvm.gt} — matrix of mean CV loss over the \code{(tau, gamma)} grid.
#'   \item \code{cvm.lambda} — matrix of mean CV loss over the \code{(lambda.sib, lambda.cou)} grid.
#'   \item \code{cme.fit} — full-data fit from \code{glmcmenet()} using the selected parameters.
#'   \item \code{select.idx}, \code{select.names} — indices and names of selected variables.
#'   \item \code{lambda.sib}, \code{lambda.cou}, \code{gamma}, \code{tau}, \code{penalty.factor} — grids/inputs used.
#'   \item \code{y}, \code{family} — echoed inputs.
#' }
#'
#' @references
#' Mak, S., & Wu, C. J. (2019). cmenet: A new method for bi-level variable selection
#' of conditional main effects. *Journal of the American Statistical Association*, 114(526), 844-856.
#'
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' n <- 50 #number of observations
#' p <- 20 #number of main effects

## Simulate model matrix for MEs and CMEs
#' set.seed(1)
#' rho <- 0 #correlation
#' ones <- matrix(1,p,p)
#' covmtx <- rho*ones+(1-rho)*diag(p)
#' latmtx <- mvrnorm(n,p,mu=rep(0,p),Sigma=covmtx) #equicorrelated cov. matrix
#' memtx <- (latmtx>=0)-(latmtx<0) #simulate model matrix for MEs
#' model.mtx <- full.model.mtx(memtx)$model.mtx #generate model matrix for MEs and CMEs
#' glist <- grouplist(model.mtx)
#' ## Set true model and generate response
#' num.act <- 2 # two siblings active
#' num.grp <- 4 # ... within four active groups
#' ind <- c()
#' for (ii in 1:num.grp){
#'  eff <- sample(seq(2*(p-1)),num.act)
#'  ind <- c(ind, p + eff + (ii-1)*(2*(p-1)))
#'}
#' colnames(model.mtx)[ind] # active CMEs

#' des.mtx <- model.mtx[,ind]
#' inter <- 0 #intercept
#' betatrue <- rep(1, length(ind))
#' xb <- inter + des.mtx %*% betatrue
#' y  <- rbinom(nrow(des.mtx), 1, 1 / (1 + exp(-xb)))

#' xme <- model.mtx[,1:p]
#' xcme <- model.mtx[,(p+1):ncol(model.mtx)]
#'
#' # weights from ridge fit (recompute with current y)
#' cv.ridge <- glmnet::cv.glmnet(cbind(xme, xcme), y, family = "binomial", alpha = 0, standardize = FALSE)
#' coefs <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
#' w  <- 1 / (abs(coefs) + 1 / n)       # element-wise
#' w[!is.finite(w)] <- 9.999e8
#' mg <- sapply(glist, function(idx) 1 / (sum(abs(coefs[idx])) + 1 / n))
#' mg[!is.finite(mg)] <- 9.999e8
#'
#' cvfit <- cv.glmcmenet(
#'   xme, xcme, y, family = "binomial", nfolds = 10,
#'   var.names=colnames(model.mtx), type.measure = "deviance",
#'   warm.str = "elastic", elastic_alpha=0.25,
#'   group.penalty = mg, penalty.factor = w
#' )
#' cvfit$params
#' cvfit$select.names
#'
#' # Fit at selected parameters:
#' fit <- cvfit$cme.fit
#' }
#'
#' @seealso \code{\link{glmcmenet}}, \code{\link{predictcme}}
#' @importFrom glmnet cv.glmnet
#' @importFrom ncvreg cv.ncvreg
#' @importFrom stats glm coef
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export


cv.glmcmenet <- function (xme, xcme, y, family = c("gaussian","binomial", "poisson"), nfolds = 10, var.names = NULL, nlambda.sib = 20,
                          nlambda.cou = 20, lambda.min.ratio = 1e-06, ngamma = 20,
                          max.gamma = 150, ntau = 20, max.tau = 0.2, tau.min.ratio = 0.001,
                          it.max = 250, it.max.cv = 25, type.measure=c("deviance","class","bic"),
                          warm.str = c("lasso","adaptive_lasso","elastic","ncvreg","NULL"), elastic_alpha = NULL,
                          penalty.factor=rep(1,ncol(xme) + ncol(xcme)), group.penalty=rep(1,2*ncol(xme)),
                          screen_ind=F)
{
  pme <- ncol(xme)
  pcme <- ncol(xcme)
  n <- nrow(xme)
  xmat <- cbind(xme, xcme)
  min.tau <- max.tau * tau.min.ratio
  min.gamma <- max(max(apply(xmat,2,function(x){8*nrow(xmat)/sum(x^2)}))+ 0.001,1/(0.125 - min.tau) + 0.001)
  act.vec <- rep(-1, ncol(xme) + ncol(xcme))
  if (warm.str == "lasso") {
    cvlas <- cv.glmnet(cbind(xme, xcme), y,family = family,alpha=1,type.measure = "deviance")
    lasfit <- cvlas$glmnet.fit
    lasind <- which(lasfit$beta[, which(cvlas$lambda ==cvlas$lambda.min)] != 0)
    act.vec[lasind] <- 1
  } else if (warm.str == "adaptive_lasso") {
    cv.ridge <- cv.glmnet(cbind(xme,xcme),y, family=family, alpha=0)
    w3 <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)
                       [, 1][2:(ncol(cbind(xme,xcme))+1)] ))^1 ## Using gamma = 1
    w3[w3[,1] == Inf] <- 999999999
    cvaplas <- cv.glmnet(cbind(xme,xcme),y,family=family, alpha=1, type.measure="deviance", penalty.factor=w3)
    aplasfit <- cvaplas$glmnet.fit
    aplasind <- which(aplasfit$beta[, which(cvaplas$lambda ==cvaplas$lambda.min)] != 0)[-1]
    act.vec[aplasind] <- 1
  }else if (warm.str == "elastic") {
    cv.ela <- cv.glmnet(cbind(xme,xcme),y, family=family, alpha=elastic_alpha, type.measure="deviance")#cv.elastic$finalModel
    fitela <- cv.ela$glmnet.fit
    elaind <- which(fitela$beta[,which(cv.ela$lambda==cv.ela$lambda.min)]!=0) #fitela$beta@i+1
    act.vec[elaind] <- 1
  } else if (warm.str == "ncvreg") {
    cvncv <- cv.ncvreg(cbind(xme,xcme),y,family = family,penalty="MCP")
    ncvfit <-cvncv$fit
    ncvind <- which(ncvfit$beta[,which(cvncv$lambda==cvncv$lambda.min)]!=0)[-1]-1
    act.vec[ncvind] <- 1
  } else if (warm.str == "NULL") {
    act.vec <- rep(-1, ncol(xme) + ncol(xcme))
  }
  start_val <- get_start(cbind(xme, xcme), y,family)
  max.lambda <- start_val$lambda_max
  lambda.sib <- exp(seq(from = log(max.lambda), to = log(max.lambda *
                                                           lambda.min.ratio), length = nlambda.sib))
  lambda.cou <- exp(seq(from = log(max.lambda), to = log(max.lambda *
                                                           lambda.min.ratio), length = nlambda.cou))
  gamma_vec = exp(seq(from = log(max.gamma), to = log(min.gamma),
                      length = ngamma - 1))
  gamma_vec = c(9.9e+35, gamma_vec)
  tau_vec = rev(exp(seq(from = log(max.tau), to = log(min.tau),
                        length = ntau)))

  # For each replicate ...

  parms1.min <- c(median(lambda.sib), median(lambda.cou))

  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }

  if (family=="binomial") {
    ind1 <- which(y==1)
    ind0 <- which(y==0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    fold1 <- 1:n1 %% nfolds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1==0] <- nfolds
    fold0[fold0==0] <- nfolds
    foldid <- double(n)
    foldid[y==1] <- sample(fold1)
    foldid[y==0] <- sample(fold0)
  }else {
    #foldid <- sample(1:n %% nfolds)
    foldid = sample(rep(seq(nfolds), length = n))
    foldid[foldid==0] <- nfolds
  }

  ## (1) Tuning gammas and tau:
  predmat = array(NA, c(n, length(tau_vec), length(gamma_vec)))

  #Perform N-fold CV
  print(paste0("Tuning gamma & tau: "))
  pb <- txtProgressBar(style = 3, width = floor(getOption("width")/2))
  for (i in seq(nfolds)) {
    setTxtProgressBar(pb, i/nfolds)
    which = (foldid == i)
    fitobj <- glmcmenet(xme = xme[!which, , drop = F], xcme = xcme[!which,, drop = F], y = y[!which],  family=family,
                        lambda.sib = parms1.min[1],lambda.cou = parms1.min[2], gamma = gamma_vec,
                        tau = tau_vec, act.vec = act.vec, penalty.factor=penalty.factor, group.penalty=group.penalty,
                        max.lambda = max.lambda, it.max = it.max.cv, screen_ind=F)
    xtest <- xmat[which, , drop = F]
    yhat <- predictcme(fitobj, xtest, type="response")
    predmat[which, , ] <- loss(fitobj,y[which],yhat,n=n,family=family,type.measure=type.measure)

  }
  cat("\n")


  #Compute sum of cv deviation
  cvm.gt <- t(apply(predmat, c(2, 3), mean))
  whichmin = argmin(cvm.gt) #select tau and gamma
  ind2.min = whichmin
  parms2.min = c(gamma_vec[whichmin[1]], tau_vec[whichmin[2]])

  ## (2) Tuning lambdas:
  lambda.sib <- exp(seq(from = log(max.lambda), to = log(max.lambda *
                                                           lambda.min.ratio), length = nlambda.sib))
  lambda.cou <- exp(seq(from = log(max.lambda), to = log(max.lambda *
                                                           lambda.min.ratio), length = nlambda.cou))
  predmat = array(NA, c(n, nlambda.sib, nlambda.cou))

  cvm <- matrix(NA, nlambda.sib, nlambda.cou)

  #Perform N-fold CV
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  print(paste0("Tuning lambdas: "))
  pb <- txtProgressBar(style = 3, width = floor(getOption("width")/2))
  for (i in seq(nfolds)) {
    setTxtProgressBar(pb, i/nfolds)
    which = (foldid == i)
    fitobj <- glmcmenet(xme = xme[!which, , drop = F], xcme = xcme[!which,
                                                                   , drop = F], y = y[!which], family=family,
                        lambda.sib = lambda.sib,lambda.cou = lambda.cou, gamma = parms2.min[1],
                        tau = parms2.min[2], act.vec = act.vec, penalty.factor=penalty.factor, group.penalty=group.penalty,
                        max.lambda = max.lambda, it.max = it.max.cv, screen_ind=screen_ind)
    xtest <- xmat[which, , drop = F]
    yhat <- predictcme(fitobj, xtest, type="response")
    predmat[which, , ] <- loss(fitobj,y[which],yhat,n=n,family=family,type.measure=type.measure)
  }
  cat("\n")
  cvm.lambda <- apply(predmat, c(2, 3), mean)
  whichmin = argmin(cvm.lambda) #select lambdas
  ind1.min = whichmin
  parms1.min = c(lambda.sib[whichmin[1]], lambda.cou[whichmin[2]])
  parms.min <- c(parms1.min, parms2.min)

  ##Summarize into list
  obj = list(y = y, family=family, lambda.sib = lambda.sib, lambda.cou = lambda.cou,
             gamma = gamma_vec, tau = tau_vec, penalty.factor=penalty.factor,
             cvm.gt = cvm.gt, cvm.lambda = cvm.lambda)
  obj$params = parms.min
  names(obj$params) <- c("lambda.sib", "lambda.cou", "gamma",
                         "tau")

  #Perform selection on full data with final parameters
  print(paste0("Fitting full data ..."))
  fitall <- glmcmenet(xme = xme, xcme = xcme, y,  family=family, lambda.sib = lambda.sib,
                      lambda.cou = lambda.cou, gamma = obj$params[3],
                      tau = obj$params[4], act.vec = act.vec, penalty.factor=penalty.factor, group.penalty=group.penalty,
                      max.lambda = max.lambda, it.max = it.max, screen_ind=screen_ind)
  obj$cme.fit <- fitall
  obj$select.idx <- which(fitall$coefficients[, which(lambda.sib ==
                                                        obj$params[1]), which(lambda.cou == obj$params[2])] !=
                            0)

  obj$select.names <- var.names[obj$select.idx]

  class(obj) = "cv.glmcme"
  return(obj)
}
