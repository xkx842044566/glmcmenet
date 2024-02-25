loss <- function(y,yhat,family,type.measure=c("deviance","class")){
  if(type.measure=="deviance"){
      val <- array(NA, dim = dim(yhat))
    if (family=="binomial") {
      val <- -2*(y*log(yhat)+(1-y)*log(1-yhat))
     }  else if (family=="poisson") {
      yly <- y*log(y)
      yly[y==0] <- 0
      val <- 2*(yly - y + yhat - y*log(yhat))
     }
  }else if(type.measure=="class"){
    val <- array(NA, dim = dim(yhat))
    if (family=="binomial") {val <- (yhat < 0.5) == y} else NULL
  }
}
