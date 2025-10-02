test_that("glmcmenet runs (gaussian, tiny grid)", {
  library(MASS)
  n <- 50 #number of observations
  p <- 12 #number of main effects

  ## Simulate model matrix for MEs and CMEs
  set.seed(1)
  rho <- 0 #correlation
  ones <- matrix(1,p,p)
  covmtx <- rho*ones+(1-rho)*diag(p)
  latmtx <- mvrnorm(n,p,mu=rep(0,p),Sigma=covmtx) #equicorrelated cov. matrix
  memtx <- (latmtx>=0)-(latmtx<0) #simulate model matrix for MEs
  model.mtx <- full.model.mtx(memtx)$model.mtx #generate model matrix for MEs and CMEs
  glist <- grouplist(model.mtx)
  ## Set true model and generate response
  num.act <- 2 # two siblings active
  num.grp <- 4 # ... within four active groups
  ind <- c()
  for (ii in 1:num.grp){
   eff <- sample(seq(2*(p-1)),num.act)
   ind <- c(ind, p + eff + (ii-1)*(2*(p-1)))
  }
  colnames(model.mtx)[ind] # active CMEs

  des.mtx <- model.mtx[,ind]
  inter <- 0 #intercept
  betatrue <- rep(1, length(ind))
  xb <- inter + des.mtx %*% betatrue
  y  <- rbinom(nrow(des.mtx), 1, 1 / (1 + exp(-xb)))

  xme <- model.mtx[,1:p]
  xcme <- model.mtx[,(p+1):ncol(model.mtx)]

  # weights from ridge fit (recompute with current y)
  cv.ridge <- glmnet::cv.glmnet(cbind(xme, xcme), y, family = "binomial", alpha = 0, standardize = FALSE)
  coefs <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
  w  <- 1 / (abs(coefs) + 1 / n)       # element-wise
  w[!is.finite(w)] <- 9.999e8
  mg <- sapply(glist, function(idx) 1 / (sum(abs(coefs[idx])) + 1 / n))
  mg[!is.finite(mg)] <- 9.999e8

  cv.glmcme <- cv.glmcmenet(
    xme, xcme, y, family = "binomial", nfolds = 10,
    var.names=colnames(model.mtx), type.measure = "deviance",
    nlambda.sib = 10, nlambda.cou = 10, ngamma = 10, ntau = 10,
    warm.str = "elastic", elastic_alpha=0.25,
    group.penalty = mg, penalty.factor = w
  )

  fit.cme <- cv.glmcme$cme.fit
  sel.cme <- cv.glmcme$select.idx

  expect_s3_class(fit.cme, "glmcme")
  expect_equal(dim(fit.cme$coefficients), c(ncol(model.mtx), 10, 10))

  # Prediction
  set.seed(1000)
  ntst <- 20
  latmtx <- mvrnorm(ntst,p,mu=rep(0,p),Sigma=covmtx)
  memtx <- (latmtx>=0)-(latmtx<0)
  tst.mtx <- full.model.mtx(memtx)$model.mtx
  xbtst <- inter + tst.mtx[, ind] %*% betatrue
  ytst  <- rbinom(length(xbtst), 1, 1 / (1 + exp(-xbtst)))
  pred.glmcme <- predictcme(fit.cme, newx = tst.mtx, type = "response")[,
                                                  which(cv.glmcme$lambda.sib == cv.glmcme$params[1]),
                                                  which(cv.glmcme$lambda.cou == cv.glmcme$params[2])]
  expect_equal(length(pred.glmcme),  20)

  #'
})
