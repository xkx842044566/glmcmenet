pbinomial <- function (eta){
  #if(is.na(eta)) {0.5
  #} else if (eta > 16) {
  if (eta > 16) {
    return(0.9999);
  } else if (eta < -16) {
    return(0.0001);
  } else {
    return(exp(eta)/(1+exp(eta)));
  }
}
