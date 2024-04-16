cv.glmcmenet <- function (xme, xcme, y, family = c("binomial", "poisson"), nfolds = 10, var.names = NULL, nlambda.sib = 20,
          nlambda.cou = 20, lambda.min.ratio = 1e-06, ngamma = 20,
          max.gamma = 150, ntau = 20, max.tau = 0.2, tau.min.ratio = 0.001,
          it.max = 250, it.max.cv = 25, type.measure=c("deviance","class"),
          warm.str = c("lasso","adaptive_lasso","elastic","ncvreg","NULL"),penalty.factor=rep(1,ncol(xme) + ncol(xcme)),
          screen_ind=F,str=F)
{
  pme <- ncol(xme)
  pcme <- ncol(xcme)
  n <- nrow(xme)
  xmat <- cbind(xme, xcme)
  min.tau <- max.tau * tau.min.ratio
  min.gamma <- max(max(apply(xmat,2,function(x){8*nrow(xmat)/sum(x^2)}))+ 0.001,1/(0.125 - min.tau) + 0.001)
  act.vec <- rep(-1, ncol(xme) + ncol(xcme))
  if (warm.str == "lasso") {
      cvlas <- cv.glmnet(cbind(xme, xcme), y,family = family,alpha=1,type.measure = type.measure)
      lasfit <- cvlas$glmnet.fit
      lasind <- which(lasfit$beta[, which(cvlas$lambda ==cvlas$lambda.min)] != 0)
      act.vec[lasind] <- 1
  } else if (warm.str == "adaptive_lasso") {
    cv.ridge <- cv.glmnet(cbind(xme,xcme),y, family=family, alpha=0)
    w3 <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)
                       [, 1][2:(ncol(cbind(xme,xcme))+1)] ))^1 ## Using gamma = 1
    w3[w3[,1] == Inf] <- 999999999
    cvaplas <- cv.glmnet(cbind(xme,xcme),y,family=family, alpha=1, type.measure=type.measure, penalty.factor=w3)
    aplasfit <- cvaplas$glmnet.fit
    aplasind <- which(aplasfit$beta[, which(cvaplas$lambda ==cvaplas$lambda.min)] != 0)[-1]
    act.vec[aplasind] <- 1
  }else if (warm.str == "elastic") {
    # control <- trainControl(method = "repeatedcv",
    #                         number = 5,
    #                         repeats = 5,
    #                         search = "random",
    #                         verboseIter = FALSE)
    # # Training ELastic Net Regression model
    # suppressWarnings(cvela <- train(y ~ .,
    #                     data = cbind(cbind(xme,xcme), y),
    #                     method = "glmnet",
    #                     family=family,
    #                     preProcess = c("center", "scale"),
    #                     tuneLength = 25,
    #                     trControl = control))
    #fitela <- glmnet(cbind(xme,xcme),y, family=family, alpha=cvela$bestTune$alpha,lambda = cvela$bestTune$lambda)#cv.elastic$finalModel
    cv.ela <- cv.glmnet(cbind(xme,xcme),y, family=family, alpha=0.25,type.measure=type.measure)#cv.elastic$finalModel
    fitela <- cv.ela$glmnet.fit
    elaind <- which(fitela$beta[,which(cv.ela$lambda==cv.ela$lambda.min)]!=0) #fitela$beta@i+1
    act.vec[elaind] <- 1
  } else if (warm.str == "ncvreg") {
    cvncv <- cv.ncvreg(cbind(xme,xcme),y,family = family,penalty="MCP")
    ncvfit <-cvncv$fit
    ncvind <- which(ncvfit$beta[,which(cvncv$lambda==cvncv$lambda.min)]!=0)[-1]-1
    act.vec[ncvind] <- 1
  } else if (warm.str == "NULL") {
    act.vec <- rep(1, ncol(xme) + ncol(xcme))
  }
  start_val <- get_start(cbind(xme, xcme), y,family)
  max.lambda <- 0.5*start_val$lambda_max
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

  # parms.min.bst <- c()
  # min.err <- 1e+30
  # cvm.gt.bst <- c()
  # cvm.lambda.bst <- c()
  parms1.min <- c(median(lambda.sib), median(lambda.cou))

  #if (!identical(sort(unique(y)), 0:1)) y <- as.double(y==max(y))
  ## Resample folds
  ## repeat{

  #foldid = sample(rep(seq(nfolds), length = n))
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
  }else if (family=="poisson") {
    foldid <- sample(1:n %% nfolds)
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
                     tau = tau_vec, act.vec = act.vec, penalty.factor=penalty.factor, max.lambda = max.lambda,
                     it.max = it.max.cv, screen_ind=F,str=str)
    xtest <- xmat[which, , drop = F]
    yhat <- predictcme(fitobj, xtest, type="response")
    predmat[which, , ] <- loss(fitobj,y[which],yhat,family=family,type.measure=type.measure)

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
                     tau = parms2.min[2], act.vec = act.vec, penalty.factor=penalty.factor, max.lambda = max.lambda,
                     it.max = it.max.cv, screen_ind=screen_ind,str=str)
    xtest <- xmat[which, , drop = F]
    yhat <- predictcme(fitobj, xtest, type="response")
    predmat[which, , ] <- loss(fitobj,y[which],yhat,family=family,type.measure=type.measure)
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
                   tau = obj$params[4], act.vec = act.vec, penalty.factor=penalty.factor, max.lambda = max.lambda,
                   it.max = it.max, screen_ind=screen_ind,str=str)
  obj$cme.fit <- fitall
  obj$select.idx <- which(fitall$coefficients[, which(lambda.sib ==
                                                        obj$params[1]), which(lambda.cou == obj$params[2])] !=
                            0)
  #if(length(obj$select.idx)>0){
  #  break
  #}
#}

  obj$select.names <- var.names[obj$select.idx]


  #refit model based on selected variables
  #temp<-cbind(cbind(xme,xcme)[,obj$select.idx],y)
  #colnames(temp)<-c(obj$select.names,"y")
  #refit <- glm(y~.,temp,family="binomial")
  #obj$cme.refit <- refit
  #obj$class.err.rate<-mean(ifelse(refit$fitted.values>0.5,1,0)!=y)

  class(obj) = "cv.glmcme"
  return(obj)
}
