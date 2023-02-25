pbinomial <- function (eta){
  if (eta > 16) {
    return(0.999999);
  } else if (eta < -16) {
    return(1e-06);
  } else {
    return(exp(eta)/(1+exp(eta)));
  }
}
