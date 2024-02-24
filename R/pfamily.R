pfamily <- function (eta,family){
  #if(is.na(eta)) {0.5
  #} else if (eta > 16) {
  if (eta > 16) {
    return(0.9999);
  } else if (eta < -16) {
    return(0.0001);
  } else {
    if(family=="binomial"){
      return(exp(eta)/(1+exp(eta)));
    }else if(family=="poisson"){
      return(exp(eta));
    }
  }
}

