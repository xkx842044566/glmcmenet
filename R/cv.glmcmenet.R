cv.glmcmenet <- function (xme, xcme, y, nfolds = 10, var.names = NULL, nlambda.sib = 10,
          nlambda.cou = 10, lambda.min.ratio = 1e-06, ngamma = 10,
          max.gamma = 150, ntau = 10, max.tau = 0.01, tau.min.ratio = 0.01,
          it.max = 250, it.max.cv = 25, warm.str = "lasso")
{
  pme <- ncol(xme)
  pcme <- ncol(xcme)
  n <- nrow(xme)
  xmat <- cbind(xme, xcme)
  min.tau <- max.tau * tau.min.ratio
  min.gamma <- max(max(apply(xmat,2,function(x){8*n/sum(x^2)})),1/(0.125 - max.tau) + 0.001)
  act.vec <- rep(-1, ncol(xme) + ncol(xcme))
  if (warm.str == "lasso") {
    # if (family == "gaussian"){
    #   cvlas <- cv.glmnet(cbind(xme, xcme), y,family = "gaussian")
    #   lasfit <- glmnet(cbind(xme, xcme), y,family = "gaussian")
    #   lasind <- which(lasfit$beta[, which(cvlas$lambda ==
    #                                         cvlas$lambda.1se)] != 0)
    #   act.vec <- rep(-1, ncol(xme) + ncol(xcme))
    #   act.vec[lasind] <- 1
    # }
    # if (family == "binomial"){
    cvlas <- cv.glmnet(cbind(xme, xcme), y,family = "binomial")
    lasfit <- glmnet(cbind(xme, xcme), y,family = "binomial")
    lasind <- which(lasfit$beta[, which(cvlas$lambda ==
                                          cvlas$lambda.1se)] != 0)
    act.vec <- rep(-1, ncol(xme) + ncol(xcme))
    act.vec[lasind] <- 1
    # }
    # else if (family == "poisson"){
    #   cvlas <- cv.glmnet(cbind(xme, xcme), y,family = "poisson")
    #   lasfit <- glmnet(cbind(xme, xcme), y,family = "poisson")
    #   lasind <- which(lasfit$beta[, which(cvlas$lambda ==
    #                                         cvlas$lambda.1se)] != 0)
    #   act.vec <- rep(-1, ncol(xme) + ncol(xcme))
    #   act.vec[lasind] <- 1
    # }
  }
  else if (warm.str == "hierNet") {
    ## need to change later
    act.vec <- rep(-1, ncol(xme) + ncol(xcme))
    warm.hn <- hierNet.path(xme, y)
    warm.cv <- hierNet.cv(warm.hn, xme, y)
    l.opt <- which(warm.hn$lamlist == warm.cv$lamhat.1se)
    me.sel <- (warm.hn$bp - warm.hn$bn)[, l.opt]
    me.idx <- which(me.sel != 0)
    int.sel <- warm.hn$th[, , l.opt]
    int.pidx <- which(int.sel > 0, arr.ind = T)
    int.pidx <- t(apply(int.pidx, 1, function(xx) {
      sort(xx)
    }))
    int.pidx <- unique(int.pidx)
    if (nrow(int.pidx) > 0) {
      for (ii in 1:nrow(int.pidx)) {
        inta <- int.pidx[ii, 1]
        intb <- int.pidx[ii, 2]
        if (inta %in% me.idx) {
          cmeind = (inta - 1) * (2 * (pme - 1)) + (intb -
                                                     2) * 2 + 1
          act.vec[pme + cmeind] = 1
          act.vec[pme + cmeind + 1] = 1
        }
        if (inta %in% me.idx) {
          cmeind = (intb - 1) * (2 * (pme - 1)) + (inta -
                                                     1) * 2 + 1
          act.vec[pme + cmeind] = 1
          act.vec[pme + cmeind + 1] = 1
        }
      }
    }
    int.nidx <- which(int.sel < 0, arr.ind = T)
    int.nidx <- t(apply(int.nidx, 1, function(xx) {
      sort(xx)
    }))
    int.nidx <- unique(int.nidx)
    if (nrow(int.nidx) > 0) {
      for (ii in 1:nrow(int.nidx)) {
        inta <- int.nidx[ii, 1]
        intb <- int.nidx[ii, 2]
        if (inta %in% me.idx) {
          cmeind = (inta - 1) * (2 * (pme - 1)) + (intb -
                                                     2) * 2 + 1
          act.vec[pme + cmeind] = 1
          act.vec[pme + cmeind + 1] = 1
        }
        if (inta %in% me.idx) {
          cmeind = (intb - 1) * (2 * (pme - 1)) + (inta -
                                                     1) * 2 + 1
          act.vec[pme + cmeind] = 1
          act.vec[pme + cmeind + 1] = 1
        }
      }
    }
  }
  max.lambda <- lambda0.cme(cbind(xme, xcme), y)
  lambda.sib <- exp(seq(from = log(max.lambda), to = log(max.lambda *
                                                           lambda.min.ratio), length = nlambda.sib))
  lambda.cou <- exp(seq(from = log(max.lambda), to = log(max.lambda *
                                                           lambda.min.ratio), length = nlambda.cou))
  gamma_vec = exp(seq(from = log(max.gamma), to = log(min.gamma),
                      length = ngamma - 1))
  gamma_vec = c(200, gamma_vec)
  tau_vec = rev(exp(seq(from = log(max.tau), to = log(min.tau),
                        length = ntau)))

  # For each replicate ...
  parms.min.bst <- c()
  min.err <- 1e+30
  cvm.gt.bst <- c()
  cvm.lambda.bst <- c()
  parms1.min <- c(median(lambda.sib), median(lambda.cou))

  ## Resample folds
  foldid = sample(rep(seq(nfolds), length = n))
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  beta0 <- rep(0, ncol(xme) + ncol(xcme))

  ## (1) Tuning gammas and tau:
  predmat = array(NA, c(n, length(tau_vec), length(gamma_vec)))

  #Perform N-fold CV
  print(paste0("Tuning gamma & tau: "))
  pb <- txtProgressBar(style = 3, width = floor(getOption("width")/2))
  for (i in seq(nfolds)) {
    setTxtProgressBar(pb, i/nfolds)
    which = (foldid == i)
    fitobj <- glmcmenet(xme = xme[!which, , drop = F], xcme = xcme[!which,
                                                                , drop = F], y = y[!which], lambda.sib = parms1.min[1],
                     lambda.cou = parms1.min[2], lambda.flg = F, gamma = gamma_vec,
                     tau = tau_vec, act.vec = act.vec, max.lambda = max.lambda,
                     it.max = it.max.cv)
    xtest <- xmat[which, , drop = F]
    predmat[which, , ] <- ifelse(predictcme(fitobj, xtest)$mu>0.5,1,0)!=y[which]

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
                                                                , drop = F], y = y[!which], lambda.sib = lambda.sib,
                     lambda.cou = lambda.cou, lambda.flg = T, gamma = parms2.min[1],
                     tau = parms2.min[2], act.vec = act.vec, max.lambda = max.lambda,
                     it.max = it.max.cv)
    xtest <- xmat[which, , drop = F]
    predmat[which, , ] <- ifelse(predictcme(fitobj, xtest)$mu>0.5,1,0)!=y[which]
      #ifelse(predictcme(fitobj, xtest)>0.5,1,0)!=y[which]
      #-2*(y[which]*log(predictcme(fitobj, xtest)$mu)+(1-y[which])*log(1-predictcme(fitobj, xtest)$mu))
  }
  cat("\n")
  cvm.lambda <- apply(predmat, c(2, 3), mean)
  whichmin = argmin(cvm.lambda) #select lambdas
  ind1.min = whichmin
  parms1.min = c(lambda.sib[whichmin[1]], lambda.cou[whichmin[2]])
  parms.min <- c(parms1.min, parms2.min)

  ##Summarize into list
  obj = list(y = y, lambda.sib = lambda.sib, lambda.cou = lambda.cou,
             gamma = gamma_vec, tau = tau_vec, cvm.gt = cvm.gt, cvm.lambda = cvm.lambda)
  obj$params = parms.min
  names(obj$params) <- c("lambda.sib", "lambda.cou", "gamma",
                         "tau")

  #Perform selection on full data with final parameters
  print(paste0("Fitting full data ..."))
  fitall <- glmcmenet(xme = xme, xcme = xcme, y, lambda.sib = lambda.sib,
                   lambda.cou = lambda.cou, lambda.flg = T, gamma = obj$params[3],
                   tau = obj$params[4], act.vec = act.vec, max.lambda = max.lambda,
                   it.max = it.max)
  obj$cme.fit <- fitall
  obj$select.idx <- which(fitall$coefficients[, which(lambda.sib ==
                                                        obj$params[1]), which(lambda.cou == obj$params[2])] !=
                            0)
  obj$select.names <- var.names[obj$select.idx]
  class(obj) = "cv.glmcme"
  return(obj)
}
