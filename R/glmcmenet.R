glmcmenet <- function (xme, xcme, y, family=c("gaussian","binomial", "poisson"),
                       lambda.sib = exp(seq(from = log(max.lambda),to = log(max.lambda * 1e-06), length = 10)),
                       lambda.cou = exp(seq(from = log(max.lambda),to = log(max.lambda * 1e-06), length = 10)),
                       max.lambda = lambda0.cme(cbind(xme,xcme), y), gamma = 1/(0.125 - tau) + 0.001, tau = 0.01,
                       act.vec = rep(1,ncol(xme) + ncol(xcme)), beta0 = rep(0, ncol(xme) + ncol(xcme)),
                       penalty.factor=rep(1,ncol(xme) + ncol(xcme)), group.penalty=rep(1,2*ncol(xme)),
                       it.max = 250, screen_ind=F)
{
  family <- match.arg(family)
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)
  if (length(gamma)>1){
    lambda.flg <- F
  } else {
    lambda.flg <- T
  }
  idx.constme <- which(apply(xme, 2, function(xx) {
    all(xx == mean(xx))
  }))
  idx.constcme <- which(apply(xcme, 2, function(xx) {
    all(xx == mean(xx))
  }))
  xme.sc <- scale(xme, center = T, scale = T)
  xme.sl <- attr(xme.sc, "scaled:scale")
  xme.ce <- attr(xme.sc, "scaled:center")
  xcme.sc <- scale(xcme, center = T, scale = T)
  xcme.sl <- attr(xcme.sc, "scaled:scale")
  xcme.ce <- attr(xcme.sc, "scaled:center")
  xme[, idx.constme] <- 0
  xcme[, idx.constcme] <- 0
  xme.sl[idx.constme] <- 1
  xcme.sl[idx.constcme] <- 1

  if (family=="gaussian") {
    ret <- cme_gaussian(xme.sc, xcme.sc, y, lambda.sib, lambda.cou, gamma,
                        tau, xme.sl, xcme.sl, beta0, act.vec, penalty.factor, group.penalty, max.lambda, it.max,
                        it_warm=3, reset=1, screen_ind)

    #Compute intercept
    if (lambda.flg){
      inter <- matrix(NA,nrow=length(lambda.sib),ncol=length(lambda.cou))
    }else{
      # inter <- matrix(NA,nrow=length(gamma),ncol=length(tau)) #ch
      inter <- matrix(NA,nrow=length(tau),ncol=length(gamma))
    }
    xmat <- cbind(xme,xcme)
    for (a in 1:nrow(inter)){
      for (b in 1:ncol(inter)){
        inter[a,b] <- mean(y - xmat%*%ret$coefficients[,a,b])
      }
    }

    #Save in object
    ret$intercept <- inter

  }else {
    ret <- cme(xme.sc, xcme.sc, y, family, lambda.sib, lambda.cou, gamma,
               tau, xme.sl, xcme.sl, beta0, act.vec, penalty.factor, group.penalty, max.lambda, it.max,
               it_warm=3, reset=1, screen_ind)
  }

  # if (lambda.flg) {
  #   inter <- matrix(NA, nrow = length(lambda.sib), ncol = length(lambda.cou))
  # }
  # else {
  #   inter <- matrix(NA, nrow = length(tau), ncol = length(gamma))
  # }
  # xmat <- cbind(xme, xcme)
  # for (a in 1:nrow(inter)) {
  #   for (b in 1:ncol(inter)) {
  #     inter[a, b] <- log(mean(y)/(1-mean(y)))-xmat %*% ret$coefficients[, a, b]
  #       ## mean(y - xmat %*% ret$coefficients[, a, b])
  #   }
  # }
  # ret$inter <- inter
  ret$xme <- xme
  ret$xcme <- xcme
  ret$y <- y
  ret$family <- family

  class(ret) = "glmcme"

  return(ret)
}
