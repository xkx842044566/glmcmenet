predictcme <- function (fit.cme, newx)
{
  dm <- dim(fit.cme$coefficients)
  eta <- array(NA, c(nrow(newx), dm[2], dm[3]))
  mu <- array(NA, c(nrow(newx), dm[2], dm[3]))
  for (i in 1:dm[2]) {
    for (j in 1:dm[3]) {
      eta[, i, j] <- fit.cme$intercept[i, j] + newx %*% matrix(fit.cme$coefficients[,
                                                                                   i, j], ncol = 1)
      mu[, i, j] <- sapply(eta[, i, j],pbinomial) # probability, mean of y, only for losgistic
    }
  }

  return(list(mu = mu, eta = eta))
}
