pbinomial <- function (eta){
  if (eta > 16) {
    return(0.9999);
  } else if (eta < -16) {
    return(0.0001);
  } else {
    return(exp(eta)/(1+exp(eta)));
  }
}
