loss <- function(y,yhat,family,type.measure=c("deviance","class")){
  if(type.measure=="deviance"){
      val <- array(NA, dim = dim(yhat))
    if (family=="binomial") {
      if (sum(y==1)) val[y==1,,] <- -2*log(yhat[y==1, , ,drop=FALSE])
      if (sum(y==0)) val[y==0,,] <- -2*log(1-yhat[y==0, , ,drop=FALSE])
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
