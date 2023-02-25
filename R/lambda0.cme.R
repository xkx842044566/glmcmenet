lambda0.cme <- function (x, y)
{
  # if(family == "gaussian"){
  # max(abs(t(x) %*% (y - mean(y))))/length(y)
  # }else
  # if(family="binomial"){
  max(4*abs(t(x) %*% (y - 0.5)))/length(y)
  # }else if(family="poisson"){
  # max(abs(t(x) %*% (y - 1)))/length(y)
  # }
}
