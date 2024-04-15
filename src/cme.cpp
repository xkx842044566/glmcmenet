// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

//Compares two doubles for equality
bool dbleq(double a, double b)
{
  return (fabs(a - b) < DBL_EPSILON);
}

//Provides ranked indices
template <typename T>  //a template function
vector<int> sort_idx(const vector<T> &v) { //that takes a single argument, v, which is a vector of elements of type T.

  // initialize original index locations
  vector<int> idx(v.size());  // return size of vector
  iota(idx.begin(), idx.end(), 0); //iota: Fills the range [first, last) with sequentially increasing values, starting with value and repetitively evaluating ++value.

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];}); //lambda expression, which is an anonymous function that can capture variables from its surrounding scope.

  return idx; //returns a vector of integers, which represents the ranked indices of the elements in v.
}

//Computes the indices for inverse pairs
vector<int> inv_ind(int pme){

  vector<int> ret(2*pme*(pme-1));

  arma::mat tmpmat(pme,pme); //temporary matrix for obtaining indices
  int count = 0;
  for (int i=0;i<pme;i++){
    for (int j=(i+1);j<pme;j++){
      tmpmat(i,j) = count;
      count++;
    }
  }
  for (int i=0;i<pme;i++){//symmetrize
    for (int j=0;j<i;j++){
      tmpmat(i,j) = tmpmat(j,i);
    }
  }

  //Retrieve indices
  int indct = 0;
  for (int i=0;i<pme;i++){
    for (int j=0;j<pme;j++){
      if (i != j){
        ret[indct] = tmpmat(i,j);
        ret[indct+1] = tmpmat(i,j);
        indct = indct + 2;
      }
    }
  }

  return (ret);
}

// Computate max absolute value
double fmax2 (double x, double y){
  double ret = 0.0;
  if (fabs(x)>fabs(y)){
    ret=x;
  }else{
    ret=y;
  }
  return (ret);
}

// Pr(y=1) for binomial
double pbinomial(double eta) {
  if (eta > 16) {
    return(0.9999);
  } else if (eta < -16) {
    return(0.0001);
  } else {
    return(exp(eta)/(1+exp(eta)));
  }
}

// E(y|X) for poisson
double ppoisson(double eta) {
  if (eta > 16) {
    return(0.9999);
  } else if (eta < -16) {
    return(0.0001);
  } else {
    return(exp(eta));
  }
}


// Cross product of y with jth column of X
double crossprod(vector<double>& X, vector<double>& yy, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*yy[i];
  return(val);
}

// Weighted cross product of y with jth column of x
double wcrossprod(vector<double>& X, vector<double>& yy, vector<double>& w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*yy[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(vector<double>&X, vector<double>&w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

// Sum of squares of jth column of X
double sqsum(vector<double>&X, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += pow(X[nn+i], 2);
  return(val);
}

double sum(vector<double>&x, int n) {
  double val=0;
  for (int i=0;i<n;i++) val += x[i];
  return(val);
}

int sum_act(vector<double>&x, vector<int>& indices) {
  int val=0;
  for(int index : indices) {
    if(index >= 0 && index < x.size()) {
      val += x[index];
    } else {
      cerr << "Warning: Index " << index << " out of bounds. Ignored." << endl;
    }
  }
  return(val);
}

vector<string> splitString(const string& str) {
  vector<string> parts;
  size_t pipe_pos = str.find("|");
  if (pipe_pos == string::npos) {
    // If no pipe is found, return the original string in parts
    parts.push_back(str);
    return parts;
  }

  // Extract parent effect (before "|")
  parts.push_back(str.substr(0, pipe_pos));

  // Extract child effect (between "|" and "+" or "-"), assuming one of these always exists at the end
  size_t sign_pos = str.length() - 1; // Position of '+' or '-' at the end
  parts.push_back(str.substr(pipe_pos + 1, sign_pos - pipe_pos - 1));

  // Optionally, extract the sign if needed
  parts.push_back(str.substr(sign_pos)); // "+" or "-"

  return parts;
}

int findind(const unordered_map<string, int>& effectIndexMap, const string& effect) {
  auto it = effectIndexMap.find(effect);
  if (it != effectIndexMap.end()) {
    return it->second; // Return the index of the effect
  } else {
    return -1; // Effect not found
  }
}


//Threshold function for ME (varying lambda)
double s_me(double inprod, double v, NumericVector& lambda, double gamma, NumericVector& delta){

  // inprod - inner product to threshold
  // lambda - penalties for reg, sib, cou and inv (ignore inv for now)
  // gamma - assumed fixed
  // delta - linearized penalties for reg, sib, cou and inv (ignore inv for now)
  // nn - the number of observations, n

  //Rank penalties and their ratios
  vector<double> lambda_r(2);
  vector<double> delta_r(2);
  vector<double> ratio(2);
  if (lambda[0] <= lambda[1]){
    lambda_r[0] = lambda[0];
    lambda_r[1] = lambda[1];
    delta_r[0] = delta[0];
    delta_r[1] = delta[1];
    ratio[0] = delta[0]/lambda[0];
    ratio[1] = delta[1]/lambda[1];
  }
  else{
    lambda_r[0] = lambda[1];
    lambda_r[1] = lambda[0];
    delta_r[0] = delta[1];
    delta_r[1] = delta[0];
    ratio[0] = delta[1]/lambda[1];
    ratio[1] = delta[0]/lambda[0];
  }

  //Compute thresholds
  double ret = 0.0;
  double sgn = 0.0;
  if (inprod < 0.0){
    sgn = -1.0;
  }
  else{
    sgn = 1.0;
  }

  if (abs(inprod) < (v*lambda_r[0]*gamma + delta_r[1]*(1-lambda_r[0]/lambda_r[1]))  ){
    if (abs(inprod) >= (delta_r[0]+delta_r[1]) ){
      ret = ( abs(inprod)-(delta_r[0]+delta_r[1]) ) / (1 - 1.0/gamma*(ratio[0]+ratio[1]));
    }
    else{
      ret = 0.0;
    }
  }
  else if (abs(inprod) < v*lambda_r[1]*gamma){
    if (abs(inprod) > (delta_r[1])){
      ret = ( abs(inprod)-(delta_r[1]) ) / (1 - 1.0/gamma*(ratio[1]) );
    }
    else{
      ret = 0.0;
    }
  }
  else{
    ret = abs(inprod);
  }
  ret=ret/v;
  return (sgn*ret);

}

// //Threshold function for CMEs (varying lambda)
// double s_cme(double inprod, NumericVector& lambda, double gamma, std::vector<double>& delta, int nn){
//
//   // inprod - inner product to threshold
//   // lambda - penalties for reg, sib, cou and inv
//   // gamma - assumed fixed
//   // delta - linearized penalties for reg, sib, cou and inv
//   // nn - the number of observations, n
//
//   double n = (double) nn;
//
//   //Rank penalties and their ratios
//   std::vector<double> lambda_r(4);
//   std::vector<double> ratio_r(4);
//   for (int i=0; i<3; i++){
//     lambda_r[i] = lambda[i];
//   }
//   std::vector<int> idx = sort_idx(lambda_r);
//   for (int i=0; i<3; i++){
//     ratio_r[i] = delta[idx[i]]/lambda_r[idx[i]];
//   }
//   std::sort(lambda_r.begin(),lambda_r.end());
//
//   //Compute thresholds
//   double ret;
//   double sgn = 0.0;
//   if (inprod < 0.0){
//     sgn = -1.0;
//   }
//   else{
//     sgn = 1.0;
//   }
//
//   if (abs(inprod) < lambda_r[0]*gamma){
//     if (abs(inprod) > (lambda_r[0]+lambda_r[1]+lambda_r[2]+lambda_r[3]) ){
//       ret = ( abs(inprod)-(lambda_r[0]+lambda_r[1]+lambda_r[2]+lambda_r[3]) ) / (1.0 - 1.0/gamma * (ratio_r[0]+ratio_r[1]+ratio_r[2]+ratio_r[3]) );
//     }
//   }
//   else if (abs(inprod) < lambda_r[1]*gamma){
//     if (abs(inprod) > (lambda_r[1]+lambda_r[2]+lambda_r[3])){
//       ret = ( abs(inprod)-(lambda_r[1]+lambda_r[2]+lambda_r[3]) ) / (1.0 - 1.0/gamma * (ratio_r[1]+ratio_r[2]+ratio_r[3]) );
//     }
//   }
//   else if (abs(inprod) < lambda_r[2]*gamma){
//     if (abs(inprod) > (lambda_r[2]+lambda_r[3])){
//       ret = ( abs(inprod)-(lambda_r[2]+lambda_r[3]) ) / (1.0 - 1.0/gamma * (ratio_r[2]+ratio_r[3]) );
//     }
//   }
//   else if (abs(inprod) < lambda_r[3]*gamma){
//     if (abs(inprod) > (lambda_r[3])){
//       ret = ( abs(inprod)-(lambda_r[3]) ) / (1.0 - 1.0/gamma * (ratio_r[3]) );
//     }
//   }
//   else{//inprod >= lambda_r[2]*gamma
//     ret = abs(inprod);
//   }
//
//   return (sgn*ret);
//
// }
//
// //Threshold function for MCP
// double s_mcp(double inprod, double lambda, double gamma){
//
//   double ret = 0.0;
//   double sgn = 0.0;
//   if (inprod < 0.0){
//     sgn = -1.0;
//   }
//   else{
//     sgn = 1.0;
//   }
//
//   if (abs(inprod) <= lambda){
//     ret = 0.0;
//   }
//   else if(abs(inprod) <= (lambda*gamma)){
//     ret = sgn * (abs(inprod)-lambda) / (1.0 - (1.0/gamma));
//   }
//   else{
//     ret = inprod;
//   }
//
//   return (ret);
// }

//MCP penalty
// [[Rcpp::export]]
double mcp(double beta, double lambda, double gamma){
  double ret = 0.0;
  if (abs(beta) <= (lambda*gamma) ){
    ret = abs(beta) - pow(beta,2.0)/(2.0*lambda*gamma);
  }
  else{
    ret = lambda*gamma/(2.0);
  }
  return(ret);
}

//KKT condition
bool kkt(double inprod,  NumericVector cur_delta){
  //Checks KKT condition for \beta=0.0
  bool ret;
  double lb = -inprod - cur_delta[0] - cur_delta[1];
  double ub = -inprod + cur_delta[0] + cur_delta[1];
  if ((0.0 >= lb)&&(0.0 <= ub)){
    ret = true; //kkt satisfied
  }
  else{
    // cout << "lb: " << lb << ", ub: " << ub << endl;
    ret = false;
  }
  return(ret);
}


//One run of coordinate descent with string detection
bool coord_des_onerun_str(int pme, int pcme, int nn, NumericVector& lambda, NumericVector& cur_delta,
                          bool dummy, double tau, double gamma,
                          vector<double>& X_me, vector<double>& X_cme, NumericVector& yy,
                          CharacterVector& names_me, CharacterVector& names_cme, unordered_map<string, int>& effectIndexMap,
                          CharacterVector& family,
                          vector<double>& delta_sib, vector<double>& delta_cou,
                          vector<bool>& act_me, vector<bool>& act_cme, double& inter,
                          vector<double>& beta_me, vector<double>& beta_cme, vector<double>& m_me, vector<double>& m_cme,
                          vector<double>& eta, vector<double>& resid, vector<double>& W, double dev){
  bool chng_flag = false;
  double cur_beta = 0.0;
  double cur_inter = 0.0;
  double xwr = 0.0;
  double xwx = 0.0;
  double inprod = 0.0;
  double v = 0.25;
  double mu = 0.0;
  //double dev = 0.0;

  // Extract the family type from CharacterVector
  string familyType = Rcpp::as<string>(family[0]);

  cur_inter= inter;

  //Update intercept
  xwr = crossprod(W, resid, nn, 0);
  xwx = sum(W,nn);
  inter = xwr/xwx + cur_inter;
  for (int i=0; i<nn; i++) {
    resid[i] -= inter - cur_inter;
    eta[i] += inter - cur_inter;
    if (familyType == "binomial"){
      mu = pbinomial(eta[i]) ;
      W[i] = fmax2(mu*(1-mu),0.0001);
    }else if(familyType == "poisson"){
      mu = ppoisson(eta[i]) ;
      W[i] = fmax2(mu,0.0001);
    }
  }

  //CD for main effects
  for (int j=0;j<pme;j++){
    //Only update if active
    if (act_me[j]){
      cur_beta = beta_me[j];

      // Updata covariates
      xwr = wcrossprod(X_me, resid, W, nn, j);
      xwx = wsqsum(X_me, W, nn, j);
      v = xwx/((double)nn);
      // inprod = inprod/((double)nn)+beta_me[j]; i.e, zj in proof
      inprod = xwr/((double)nn)+v*beta_me[j]; //checked to pod from update eqn (mod from above eqn)

      string me = Rcpp::as<string>(names_me[j]);
      int delta_ind = findind(effectIndexMap,me);

      //Update cur_delta
      cur_delta[0] = delta_sib[delta_ind] * m_me[j];
      cur_delta[1] = delta_cou[delta_ind] * m_me[j];

      //Perform ME thresholding
      beta_me[j] = s_me(inprod,v,lambda,gamma,cur_delta);

      // Update eta, mu, weight and delta
      if (!dbleq(beta_me[j],cur_beta)){ // if beta changed...

        //Update resid eta, mu, weight
        for (int k=0;k<nn;k++){
          resid[k] -= X_me[j*nn+k]*(beta_me[j]-cur_beta);
          eta[k] += X_me[j*nn+k]*(beta_me[j]-cur_beta);
          if (familyType == "binomial"){
            mu = pbinomial(eta[k]) ;
            W[k] = fmax2(mu*(1-mu),0.0001);
          }else if(familyType == "poisson"){
            mu = ppoisson(eta[k]) ;
            W[k] = fmax2(mu,0.0001);
          }
          //v += (X_me[j*nn+k]*W[k]*X_me[j*nn+k])/((double)nn);
        }
        // xwx = wsqsum(X_me, W, nn, j);
        // v = xwx/((double)nn);

        //Update deltas
        double offset_sib = mcp(beta_me[j],lambda[0],gamma)-mcp(cur_beta,lambda[0],gamma); // new - old
        double offset_cou = mcp(beta_me[j],lambda[1],gamma)-mcp(cur_beta,lambda[1],gamma);
        delta_sib[delta_ind] = delta_sib[delta_ind] * (exp(-(tau/lambda[0]) * m_me[j] * offset_sib ));
        delta_cou[delta_ind] = delta_cou[delta_ind] * (exp(-(tau/lambda[1]) * m_me[j] * offset_cou ));

        //Update flag
        chng_flag = true;

      }
    }
  }


  //Update intercept
  xwr = crossprod(W, resid, nn, 0);
  xwx = sum(W,nn);
  inter = xwr/xwx + cur_inter;
  for (int i=0; i<nn; i++) {
    resid[i] -= inter - cur_inter;
    eta[i] += inter - cur_inter;
    if (familyType == "binomial"){
      mu = pbinomial(eta[i]) ;
      W[i] = fmax2(mu*(1-mu),0.0001);
    }else if(familyType == "poisson"){
      mu = ppoisson(eta[i]) ;
      W[i] = fmax2(mu,0.0001);
    }
  }
  cur_beta = 0.0;
  inprod = 0.0;

  //CD for CME effects
  //for (int j=0;j<pme;j++){ //parent effect
  //  for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
  for (int j=0;j<pcme;j++){

    //int cmeind = j*(2*(pme-1))+k; //index for CME
    //int condind = 0; //index for condition

    //if (act_cme[cmeind]){
    //  int condind = floor((double)k/2.0);
    //  if (condind >= j){
    //    condind ++;
    //  }
    if (act_cme[j]){
      string cme = Rcpp::as<string>(names_cme[j]);
      auto parts = splitString(cme);
      int sib_ind = findind(effectIndexMap,parts[0]);
      int cou_ind = findind(effectIndexMap,parts[1]);

      cur_beta = beta_cme[j]; //beta_cme ordered by parent effects, then condition
      //cur_inter= inter;

      //Update cur_delta
      cur_delta[0] = delta_sib[sib_ind] * m_cme[j];
      cur_delta[1] = delta_cou[cou_ind] * m_cme[j];

      // Updata covariates
      xwr = wcrossprod(X_cme, resid, W, nn, j);
      xwx = wsqsum(X_cme, W, nn, j);
      v = xwx/((double)nn);
      // inprod = inprod/((double)nn)+beta_me[j]; i.e, zj in proof
      inprod = xwr/((double)nn)+v*beta_cme[j]; //checked to pod from update eqn (mod from above eqn)

      //Perform CME thresholding
      beta_cme[j] = s_me(inprod,v,lambda,gamma,cur_delta);

      //Update eta, mu, weight and delta
      if (!dbleq(beta_cme[j],cur_beta)){ // if beta changed...

        //Update resid eta, mu, weight
        for (int ll=0;ll<nn;ll++){
          resid[ll] -= X_cme[j*nn+ll]*(beta_cme[j]-cur_beta);
          eta[ll] += X_cme[j*nn+ll]*(beta_cme[j]-cur_beta);
          if (familyType == "binomial"){
            mu = pbinomial(eta[ll]) ;
            W[ll] = fmax2(mu*(1-mu),0.0001);
          }else if(familyType == "poisson"){
            mu = ppoisson(eta[ll]) ;
            W[ll] = fmax2(mu,0.0001);
          }
          //v += (X_me[j*nn+k]*W[k]*X_me[j*nn+k])/((double)nn);
        }

        //Update deltas
        double offset_sib = mcp(beta_cme[j],lambda[0],gamma)-mcp(cur_beta,lambda[0],gamma); // new - old
        double offset_cou = mcp(beta_cme[j],lambda[1],gamma)-mcp(cur_beta,lambda[1],gamma); // new - old
        delta_sib[sib_ind] = delta_sib[sib_ind] * (exp(-(tau/lambda[0]) * m_cme[j] * offset_sib )); // update delta for siblings
        delta_cou[cou_ind] = delta_cou[cou_ind] * (exp(-(tau/lambda[1]) * m_cme[j] * offset_cou )); // update delta for cousins

        //Update flag
        chng_flag = true;

      }



      // //Reduce A|B+ and A|B- to A
      //  if (abs(beta_cme[j]) > 0.0){ //if current CME is active
      //      string oppo_eff = parts[0] + "|" + parts[1] + (parts[2] == "+" ? "-" : "+");
      //      int j_oppo = -1; // Start with an invalid index
      //      for (int ii = 0; ii < pcme; ii++) {
      //        if (names_cme[ii] == oppo_eff) {
      //          j_oppo = ii;
      //          break; // Stop the loop once the effect is found
      //        }
      //      }
      //      int me_ind = -1; // Start with an invalid index
      //      for (int ii = 0; ii < pcme; ii++) {
      //        if (names_me[ii] == parts[0]) {
      //          me_ind = ii;
      //          break; // Stop the loop once the effect is found
      //        }
      //      }
      //      if (abs(beta_cme[j_oppo]) > 0.0){ //if cme .|.- is also in model...
      //
      //        double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;
      //
      //        if ( abs(beta_cme[j]) > abs(beta_cme[j_oppo]) ){// if abs(.|.+) > abs(.|.-)
      //          chg = beta_cme[j_oppo]; // change
      //          cur_beta_cme1 = beta_cme[j]; // current beta cme 1
      //          cur_beta_cme2 = beta_cme[j_oppo]; // current beta cme 2
      //          beta_cme[j] -= chg; // update larger CME
      //          beta_cme[j_oppo] = 0.0; // remove smaller CME
      //
      //          cur_beta_me = beta_me[me_ind]; // current beta me
      //          beta_me[me_ind] += chg; // update ME with smaller CME
      //        }else{// if abs(.|.+) < abs(.|.-)
      //          chg = beta_cme[j]; // change
      //          cur_beta_cme1 = beta_cme[j]; // current beta cme 1
      //          cur_beta_cme2 = beta_cme[j_oppo]; // current beta cme 2
      //          beta_cme[j_oppo] -= chg; // update larger CME
      //          beta_cme[j] = 0.0; // remove smaller CME
      //
      //          cur_beta_me = beta_me[me_ind]; // current beta me
      //          beta_me[me_ind] += chg; // update ME with smaller CME
      //        }
      //
      //        //Update deltas and flag
      //        double offset_sib = mcp(beta_me[me_ind],lambda[0],gamma)-mcp(cur_beta_me,lambda[0],gamma); // new - old (for me)
      //        double offset_cou = mcp(beta_me[me_ind],lambda[1],gamma)-mcp(cur_beta_me,lambda[1],gamma);
      //        delta_sib[sib_ind] = delta_sib[sib_ind] * (exp(-(tau/lambda[0]) * offset_sib ));
      //        delta_cou[sib_ind] = delta_cou[sib_ind] * (exp(-(tau/lambda[1]) * offset_cou ));
      //
      //        offset_sib = mcp(beta_cme[j],lambda[0],gamma)-mcp(cur_beta_cme1,lambda[0],gamma); // new - old (for .|.+)
      //        offset_cou = mcp(beta_cme[j],lambda[1],gamma)-mcp(cur_beta_cme1,lambda[1],gamma);
      //        delta_sib[sib_ind] = delta_sib[sib_ind] * (exp(-(tau/lambda[0]) * offset_sib ));
      //        delta_cou[cou_ind] = delta_cou[cou_ind] * (exp(-(tau/lambda[1]) * offset_cou ));
      //
      //        offset_sib = mcp(beta_cme[j_oppo],lambda[0],gamma)-mcp(cur_beta_cme2,lambda[0],gamma); // new - old (for .|.-)
      //        offset_cou = mcp(beta_cme[j_oppo],lambda[1],gamma)-mcp(cur_beta_cme2,lambda[1],gamma);
      //        delta_sib[sib_ind] = delta_sib[sib_ind] * (exp(-(tau/lambda[0]) * offset_sib ));
      //        delta_cou[cou_ind] = delta_cou[cou_ind] * (exp(-(tau/lambda[1]) * offset_cou ));
      //
      //        //residuals shouldn't change
      //
      //      }
      //   } // else{ //cme is .|.-
      //         if (abs(beta_cme[cmeind-1]) > 0.0){ //if cme .|.+ is also in model...

      //           double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

      //           if ( abs(beta_cme[cmeind]) > abs(beta_cme[cmeind-1]) ){// if abs(.|.+) < abs(.|.-)
      //             chg = beta_cme[cmeind-1]; // change
      //             cur_beta_me = beta_me[j]; // current beta me
      //             cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
      //             cur_beta_cme2 = beta_cme[cmeind-1]; // current beta cme 2
      //             beta_me[j] += chg; // update ME with smaller CME
      //             beta_cme[cmeind] -= chg; // update larger CME
      //             beta_cme[cmeind-1] = 0.0; // remove smaller CME
      //           }else{// if abs(.|.+) > abs(.|.-)
      //             chg = beta_cme[cmeind]; // change
      //             cur_beta_me = beta_me[j]; // current beta me
      //             cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
      //             cur_beta_cme2 = beta_cme[cmeind-1]; // current beta cme 2
      //             beta_me[j] += chg; // update ME with smaller CME
      //             beta_cme[cmeind-1] -= chg; // update larger CME
      //             beta_cme[cmeind] = 0.0; // remove smaller CME
      //           }

      //           //Update deltas and flag
      //           double offset_sib = mcp(beta_me[j],lambda[0],gamma)-mcp(cur_beta_me,lambda[0],gamma); // new - old (for me)
      //           double offset_cou = mcp(beta_me[j],lambda[1],gamma)-mcp(cur_beta_me,lambda[1],gamma);
      //           delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda[0]) * offset_sib ));
      //           delta_cou[j] = delta_cou[j] * (exp(-(tau/lambda[1]) * offset_cou ));

      //           offset_sib = mcp(beta_cme[cmeind],lambda[0],gamma)-mcp(cur_beta_cme1,lambda[0],gamma); // new - old (for .|.+)
      //           offset_cou = mcp(beta_cme[cmeind],lambda[1],gamma)-mcp(cur_beta_cme1,lambda[1],gamma);
      //           delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda[0]) * offset_sib ));
      //           delta_cou[condind] = delta_cou[condind] * (exp(-(tau/lambda[1]) * offset_cou ));

      //           offset_sib = mcp(beta_cme[cmeind-1],lambda[0],gamma)-mcp(cur_beta_cme2,lambda[0],gamma); // new - old (for .|.-)
      //           offset_cou = mcp(beta_cme[cmeind-1],lambda[1],gamma)-mcp(cur_beta_cme2,lambda[1],gamma);
      //           delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda[0]) * offset_sib ));
      //           delta_cou[condind] = delta_cou[condind] * (exp(-(tau/lambda[1]) * offset_cou ));

      //           //residuals shouldn't change

      //         }
      //       }
      //     }
    }

  }


  return (chng_flag);

}


//One run of coordinate descent under default structure
bool coord_des_onerun(int pme, int nn, NumericVector& K1,
                      NumericVector& lambda, NumericVector& cur_delta,
                      bool dummy, double tau, double gamma,
                      vector<double>& X, NumericVector& yy,
                      CharacterVector& family,
                      vector<double>& delta_sib, vector<double>& delta_cou,
                      vector<bool>& act, double& inter,
                      vector<double>& beta, vector<double>& mg,
                      vector<double>& eta, vector<double>& resid, vector<double>& W, double dev){
  bool chng_flag = false;
  double cur_beta = 0.0;
  double cur_inter = 0.0;
  double xwr = 0.0;
  double xwx = 0.0;
  double inprod = 0.0;
  double v = 0.25;
  double mu = 0.0;
  int J = K1.size() - 1;

  //double dev = 0.0;

  // Extract the family type from CharacterVector
  string familyType = Rcpp::as<string>(family[0]);

  cur_inter= inter;

  //Update intercept
  xwr = crossprod(W, resid, nn, 0);
  xwx = sum(W,nn);
  inter = xwr/xwx + cur_inter;
  for (int i=0; i<nn; i++) {
    resid[i] -= inter - cur_inter;
    eta[i] += inter - cur_inter;
    if (familyType == "binomial"){
      mu = pbinomial(eta[i]) ;
      W[i] = fmax2(mu*(1-mu),0.0001);
    }else if(familyType == "poisson"){
      mu = ppoisson(eta[i]) ;
      W[i] = fmax2(mu,0.0001);
    }
  }

  for (int g=0; g<J; g++) {

    //int K = K1[g+1] - K1[g];

    //CD for covariates in group g
    for (int j=K1[g];j<K1[g+1];j++){
      //Only update if active
      if (act[j]){
        cur_beta = beta[j];

        // Updata covariates
        xwr = wcrossprod(X, resid, W, nn, j);
        xwx = wsqsum(X, W, nn, j);
        v = xwx/((double)nn);
        // inprod = inprod/((double)nn)+beta_me[j]; i.e, zj in proof
        inprod = xwr/((double)nn)+v*beta[j]; //checked to pod from update eqn (mod from above eqn)

        //string me = Rcpp::as<string>(names_me[j]);
        //int delta_ind = findind(effectIndexMap,me);

        //Update cur_delta
        int gind;
        if (g % 2 == 0) { //if this is sibling group
          gind = floor((double)g/2.0);
          cur_delta[0] = delta_sib[gind] * mg[g];
          cur_delta[1] = delta_cou[gind] * mg[g+1];
        } else { //if this is cousin group
          gind = floor((double)g/2.0); //index for cousin group
          cur_delta[0] = delta_sib[gind] * mg[g-1];
          cur_delta[1] = delta_cou[gind] * mg[g];
        }

        //Perform ME thresholding
        beta[j] = s_me(inprod,v,lambda,gamma,cur_delta);

        // Update eta, mu, weight and delta
        if (!dbleq(beta[j],cur_beta)){ // if beta changed...

          //Update resid eta, mu, weight
          for (int k=0;k<nn;k++){
            resid[k] -= X[j*nn+k]*(beta[j]-cur_beta);
            eta[k] += X[j*nn+k]*(beta[j]-cur_beta);
            if (familyType == "binomial"){
              mu = pbinomial(eta[k]) ;
              W[k] = fmax2(mu*(1-mu),0.0001);
            }else if(familyType == "poisson"){
              mu = ppoisson(eta[k]) ;
              W[k] = fmax2(mu,0.0001);
            }
            //v += (X_me[j*nn+k]*W[k]*X_me[j*nn+k])/((double)nn);
          }
          // xwx = wsqsum(X_me, W, nn, j);
          // v = xwx/((double)nn);

          //Update deltas
          double offset_sib = mcp(beta[j],lambda[0],gamma)-mcp(cur_beta,lambda[0],gamma); // new - old
          double offset_cou = mcp(beta[j],lambda[1],gamma)-mcp(cur_beta,lambda[1],gamma);

          if (g % 2 == 0) {
            delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
            delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou ))
          } else {
            delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib ))
            delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
          }

          //Update flag
          chng_flag = true;

        }

        //Reduce A|B+ and A|B- to A
        //
        if (abs(beta[j]) > 0.0 && j != K1[g]){ //if current CME is active
          int id = j - K1[g];
          if (id % 2 == 1){ //cme is .|.+
            if (abs(beta[j+1]) > 0.0){ //if cme .|.- is also in model...

              double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

              if ( abs(beta[j]) > abs(beta[j+1]) ){// if abs(.|.+) > abs(.|.-)
                chg = beta[j+1]; // change
                cur_beta_me = beta[K1[g]]; // current beta me
                cur_beta_cme1 = beta[j]; // current beta cme 1
                cur_beta_cme2 = beta[j+1]; // current beta cme 2
                beta[K1[g]] += chg; // update ME with smaller CME
                beta[j] -= chg; // update larger CME
                beta[j+1] = 0.0; // remove smaller CME
              }else{// if abs(.|.+) < abs(.|.-)
                chg = beta[j]; // change
                cur_beta_me = beta[K1[g]]; // current beta me
                cur_beta_cme1 = beta[j]; // current beta cme 1
                cur_beta_cme2 = beta[j+1]; // current beta cme 2
                beta[K1[g]] += chg; // update ME with smaller CME
                beta[j+1] -= chg; // update larger CME
                beta[j] = 0.0; // remove smaller CME
              }

              //Update deltas and flag
              double offset_sib = mcp(beta[K1[g]],lambda[0],gamma)-mcp(cur_beta_me,lambda[0],gamma); // new - old (for me)
              double offset_cou = mcp(beta[K1[g]],lambda[1],gamma)-mcp(cur_beta_me,lambda[1],gamma);
              if (g % 2 == 0) {
                delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
                delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou ))
              } else {
                delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib ))
                delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
              }

              offset_sib = mcp(beta[j],lambda[0],gamma)-mcp(cur_beta_cme1,lambda[0],gamma); // new - old (for .|.+)
              offset_cou = mcp(beta[j],lambda[1],gamma)-mcp(cur_beta_cme1,lambda[1],gamma);
              if (g % 2 == 0) {
                delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
                delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou ))
              } else {
                delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib ))
                delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
              }

              offset_sib = mcp(beta[j+1],lambda[0],gamma)-mcp(cur_beta_cme2,lambda[0],gamma); // new - old (for .|.-)
              offset_cou = mcp(beta[j+1],lambda[1],gamma)-mcp(cur_beta_cme2,lambda[1],gamma);
              if (g % 2 == 0) {
                delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
                delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou ))
              } else {
                delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib ))
                delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
              }
              //residuals shouldn't change

            }
          }else{ //cme is .|.-
            if (abs(beta[j-1]) > 0.0){ //if cme .|.+ is also in model...

              double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

              if ( abs(beta[j]) > abs(beta[j-1]) ){// if abs(.|.+) < abs(.|.-)
                chg = beta[j-1]; // change
                cur_beta_me = beta[K1[g]]; // current beta me
                cur_beta_cme1 = beta[j]; // current beta cme 1
                cur_beta_cme2 = beta[j-1]; // current beta cme 2
                beta[K1[g]] += chg; // update ME with smaller CME
                beta[j] -= chg; // update larger CME
                beta[j-1] = 0.0; // remove smaller CME
              }else{// if abs(.|.+) > abs(.|.-)
                chg = beta[j]; // change
                cur_beta_me = beta[K1[g]]; // current beta me
                cur_beta_cme1 = beta[j]; // current beta cme 1
                cur_beta_cme2 = beta[j-1]; // current beta cme 2
                beta[K1[g]] += chg; // update ME with smaller CME
                beta[j-1] -= chg; // update larger CME
                beta[j] = 0.0; // remove smaller CME
              }

              //Update deltas and flag
              double offset_sib = mcp(beta[K1[g]],lambda[0],gamma)-mcp(cur_beta_me,lambda[0],gamma); // new - old (for me)
              double offset_cou = mcp(beta[K1[g]],lambda[1],gamma)-mcp(cur_beta_me,lambda[1],gamma);
              if (g % 2 == 0) {
                delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
                delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou )) *
              } else {
                delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib )) *
                delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
              }

              offset_sib = mcp(beta[j],lambda[0],gamma)-mcp(cur_beta_cme1,lambda[0],gamma); // new - old (for .|.+)
              offset_cou = mcp(beta[j],lambda[1],gamma)-mcp(cur_beta_cme1,lambda[1],gamma);
              if (g % 2 == 0) {
                delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
                delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou )) *
              } else {
                delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib )) *
                delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
              }

              offset_sib = mcp(beta[j-1],lambda[0],gamma)-mcp(cur_beta_cme2,lambda[0],gamma); // new - old (for .|.-)
              offset_cou = mcp(beta[j-1],lambda[1],gamma)-mcp(cur_beta_cme2,lambda[1],gamma);
              if (g % 2 == 0) {
                delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * offset_sib )) * mg[g];
                delta_cou[gind] = delta_cou[gind] * mg[g+1]; //(exp(-(tau/lambda[1]) * offset_cou )) *
              } else {
                delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * offset_sib )) *
                delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * offset_cou )) * mg[g];
              }

              //residuals shouldn't change

            }
          }
        }
      }
    }
  }


  return (chng_flag);

}


// [[Rcpp::export]]
List cme_str(NumericMatrix& XX_me, NumericMatrix& XX_cme, NumericVector& yy, CharacterVector& family,
             NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec,
             NumericVector& gamma_vec, NumericVector& tau_vec,
             NumericVector& XX_me_sl, NumericVector& XX_cme_sl, NumericVector& beta_vec, NumericVector& act_vec,NumericVector& multiplier,
             double lambda_max, int it_max, int it_warm, int reset, bool screen_ind) {
  // // [[Rcpp::plugins(openmp)]]
  //------------------------------------------------------------
  // XX - Full model matrix including both ME and CME effects (assume normalized)
  // yy - Response vector of length nn
  // family- family of GLM
  // lambda_sib_vec - Vector of sibling penalties (decr. sequence)
  // lambda_cou_vec - Vector of cousin penalties (decr. sequence)
  // tau - Exponential penalty parameter
  // gamma - MC+ non-convex penalty parameter
  // beta_vec - Initial beta value
  // it_max - Maximum iterations for coordinate descent
  //------------------------------------------------------------

  //Variable initialization
  int pme = XX_me.ncol(); //# of MEs
  int pcme = XX_cme.ncol(); //# of CMEs
  int nn = XX_me.nrow(); //# of observations
  int nlambdasib = lambda_sib_vec.size();
  int nlambdacou = lambda_cou_vec.size();
  int it_inner = 0;
  int it_max_reset = it_max / reset;
  bool cont = true;
  bool chng_flag = false;
  double v = 0.25;
  double mu = 0;

  // Extract the family type from CharacterVector
  string familyType = Rcpp::as<string>(family[0]);

  //Vectorize model matrices
  vector<double> X_me(nn*pme); //for ME
  vector<double> X_cme(nn*pcme); //for CME
  for (int i=0;i<pme;i++){
    for (int j=0;j<nn;j++){
      X_me[i*nn+j] = XX_me(j,i);
    }
  }
  for (int i=0;i<pcme;i++){
    for (int j=0;j<nn;j++){
      X_cme[i*nn+j] = XX_cme(j,i);
    }
  }
  // store colnames name of ME and CME
  CharacterVector names_me = colnames(XX_me);
  CharacterVector names_cme = colnames(XX_cme);

  //Check whether lambda is to be iterated or not
  bool lambda_it;
  int niter_1; //Number to iterate first
  int niter_2; //Number to iterate next
  if (gamma_vec.size()>1){ //Iterate on gamma and tau
    lambda_it = false;
    // niter_1 = gamma_vec.size(); //ch
    // niter_2 = tau_vec.size();
    niter_1 = tau_vec.size();
    niter_2 = gamma_vec.size();
  }
  else{
    lambda_it = true;
    niter_1 = nlambdasib;
    niter_2 = nlambdacou;
  }

  //Containers for beta and active set (alpha)
  arma::cube beta_cube(pme+pcme,niter_1,niter_2); //betas to return
  arma::cube delta_sib_cube(pme,niter_1,niter_2); //deltas to return
  arma::cube delta_cou_cube(pme,niter_1,niter_2); //deltas to return
  arma::mat nz(niter_1,niter_2);
  arma::mat inter_mat(niter_1,niter_2); //intercept to return
  arma::mat dev_mat(niter_1,niter_2); //deviation to return
  arma::mat beta_mat(pme+pcme,niter_1);
  arma::mat delta_sib_mat(pme,niter_1);
  arma::mat delta_cou_mat(pme,niter_1);

  vector<double> beta_me(pme,0.0); //for MEs
  for (int i=0;i<pme;i++){
    beta_me[i] = beta_vec[i];
  }
  vector<double> beta_cme(pcme,0.0); //for CMEs
  for (int i=0;i<pcme;i++){
    beta_cme[i] = beta_vec[pme+i];
  }
  //double cur_beta = 0.0; //running beta

  //Containers for siblings or cousion family
  vector<int> sibind;
  vector<int> couind;

  //Set all factors as active to begin
  vector<bool> act_me(pme,true); //Current active set
  vector<bool> act_cme(pcme,true);
  vector<bool> scr_me(pme,true); //Screened active set
  vector<bool> scr_cme(pcme,true);
  //bool kkt_bool;

  set<string> uniqueEffects;

  // Include main effects
  for (int i = 0; i < pme; i++) {
    uniqueEffects.insert(Rcpp::as<string>(names_me[i]));
  }

  // Include parent and child effects from conditional main effects
  for (int i = 0; i < pcme; i++) {
    auto parts = splitString(Rcpp::as<string>(names_cme[i]));
    for (const auto& part : parts) {
      uniqueEffects.insert(part);
    }
  }

  unordered_map<string, int> effectIndexMap;
  int index = 0;
  for (const auto& effect : uniqueEffects) {
    effectIndexMap[effect] = index++;
  }

  // Containers for linearized slopes Delta
  vector<double> delta_sib(uniqueEffects.size(), 0.0); // Linearized penalty for siblings (sib(A), sib(B), ...)
  vector<double> delta_cou(uniqueEffects.size(), 0.0); // Linearized penalty for cousins (cou(A), cou(B), ...)

  vector<double> m_me(pme,1);
  vector<double> m_cme(pcme,1);
  for (int i=0;i<pme;i++){
    m_me[i] = multiplier[i];
  }
  for (int i=0;i<pcme;i++){
    m_cme[i] = multiplier[pme+i];
  }

  //Containers for linearized slopes Delta
  //vector<double> delta_sib(pme); //Linearized penalty for siblings (sib(A), sib(B), ...)
  //vector<double> delta_cou(pme); //Linearized penalty for cousins (cou(A), cou(B), ...)
  NumericVector lambda(2); //Current penalties
  NumericVector cur_delta(2); //Current delta vector
  double gamma;
  double tau;
  lambda[0] = lambda_sib_vec[0];
  lambda[1] = lambda_cou_vec[0];
  gamma = gamma_vec[0];
  tau = tau_vec[0];

  //Update resid (eta)
  vector<double> eta(nn);
  vector<double> W(nn);
  vector<double> resid(nn); //Residual vector
  arma::cube resid_cube(nn,niter_1,niter_2);
  arma::mat resid_mat(nn,niter_1);
  arma::cube scr_cube(pme+pcme,niter_1,niter_2); //screening vector
  arma::mat scr_mat(pme+pcme,niter_1);

  double nullDev = 0.0;
  double dev = 0.0;
  double inter= 0.0; //for intercept
  double inprod = 0.0; //inner product
  double cj = 0.0;
  double vj = 0.0;
  double thresh = 0.0; //threshold for screening
  //int size = 0;
  int num_act = 0;
  int num_scr = 0;

  double ymean = 0.0;
  for (int i=0;i<nn;i++){
    ymean += (1.0/(double)nn)*yy[i];
  }
  if (familyType == "binomial") {
    inter = log(ymean/(1-ymean));
    for (int i=0; i<nn; i++) {
      eta[i] = log(ymean/(1-ymean)); //
      mu = pbinomial(eta[i]) ;
      W[i] = fmax2(mu*(1-mu),0.0001);
      resid[i] = (yy[i]-mu)/W[i];
      nullDev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
      dev = nullDev;
    }
  } else if (familyType == "poisson") {
    inter = log(ymean);
    for (int i=0;i<nn;i++) {
      eta[i] = log(ymean);
      mu = ppoisson(eta[i]) ;
      W[i] = fmax2(mu,0.0001);
      resid[i] = (yy[i]-mu)/W[i];
      if (yy[i]!=0) nullDev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
      else nullDev += 2*ymean;
      dev = nullDev;
    }
  }

  //vector<bool> kkt_v_me(pme,true);
  //vector<bool> kkt_v_cme(pcme,true); //KKT checks

  // Optimize for each penalty combination
  // #pragma omp parallel for
  for (int b=0; b<niter_2; b++){ //iterate over cousins...

    for (int a=0; a<niter_1; a++){ //iterate over siblings...

      // cout << "Tuning ... a = " << a << ", b = " << b << endl;

      //Update iteration variables
      if (lambda_it){
        lambda[0] = lambda_sib_vec[a];
        lambda[1] = lambda_cou_vec[b];
      }
      else{
        // gamma = gamma_vec[a]; //ch
        // tau = tau_vec[b];
        tau = tau_vec[a];
        gamma = gamma_vec[b];
      }

      //Update screen set
      if (lambda_it && screen_ind && a>0 && b>0){

        int num_scr = 0;
        for (int j=0;j<pme;j++){

          if(act_me[j] && scr_me[j]){

            string me = Rcpp::as<string>(names_me[j]);

            sibind.clear();
            couind.clear();

            // Iterate over each column name in XX_cme to find siblings for the current main effect
            for (int i = 0; i < pcme; i++) {
              string colName = Rcpp::as<string>(names_cme[i]);
              auto parts = splitString(colName); // Use splitString function here


              // Proceed only if "|" is found and there are two parts (before and after "|")
              // if (parts.size() == 2) {
              string beforePipe = parts[0];
              string afterPipe = parts[1];

              // If the main effect matches the part before "|", it's a sibling
              if (beforePipe == me) {
                sibind.push_back(i); // 0-based indexing
              }

              // If the main effect matches the part after "|", it's a cousin
              if (afterPipe == me) {
                couind.push_back(i); // 0-based indexing
              }
              //}
            }


            cj = wcrossprod(X_me, resid, W, nn, j);
            vj = wsqsum(X_me, W, nn, j)/((double)nn);

            if (sum_act(beta_cme,sibind)==0) {
              if (sum_act(beta_cme,couind)==0) {
                thresh = max(lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[0]-lambda_sib_vec[a-1]),
                             lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[1]-lambda_cou_vec[b-1]));
                if (abs(cj) > thresh) {
                  scr_me[j] = true;
                  num_scr ++;
                }else{
                  scr_me[j] = false;
                }
              } else if (sum_act(beta_cme,couind)>0) {
                thresh = lambda[0]+cur_delta[1]+vj*gamma/(vj*gamma-cur_delta[1]/lambda[1]-1)*(lambda[0]-lambda_sib_vec[a-1]);
                if (abs(cj) > thresh) {
                  scr_me[j] = true;
                  num_scr ++;
                }else{
                  scr_me[j] = false;
                }
              }
            }
            else {
              if (sum_act(beta_cme,couind)==0) {
                thresh = cur_delta[0]+lambda[1]+vj*gamma/(vj*gamma-cur_delta[0]/lambda[0]-1)*(lambda[1]-lambda_cou_vec[b-1]);
                if (abs(cj) > thresh) {
                  scr_me[j] = true;
                  num_scr ++;
                }else{
                  scr_me[j] = false;
                }
              }
            }
          }
        }

        //for (int j=0;j<pme;j++){ //parent effect
        //for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
        for (int j=0; j<pcme; j++){

          sibind.clear();
          couind.clear();

          if(act_cme[j] && scr_cme[j]){

            string cme = Rcpp::as<string>(names_cme[j]);
            auto parts = splitString(cme);
            // if (parts.size() != 2) continue; // Skip if not a conditional effect

            string parent = parts[0], child = parts[1];

            // siblings[colName].push_back(parent); // Add parent as a sibling
            for (int k = 0; k < pcme; k++) {
              if (j == k) continue; // Skip self
              string othercme = Rcpp::as<string>(names_cme[k]);
              auto otherParts = splitString(othercme);
              // if (otherParts.size() != 2) continue;
              if (otherParts[0] == parent) {
                sibind.push_back(k); // Siblings: same parent.
              }
              if (otherParts[1] == child) {
                couind.push_back(k); // Cousins: same child.
              }

            }


            cj = wcrossprod(X_cme, resid, W, nn, j);
            vj = wsqsum(X_cme, W, nn, j)/((double)nn);

            if (sum_act(beta_cme,sibind)==0) {
              if (sum_act(beta_cme,couind)==0) {
                thresh = max(lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[0]-lambda_sib_vec[a-1]),
                             lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[1]-lambda_cou_vec[b-1]));
                if (abs(cj) > thresh) {
                  scr_cme[j] = true;
                  num_scr ++;
                }else{
                  scr_cme[j] = false;
                }
              } else if (sum_act(beta_cme,couind)>0) {
                thresh = lambda[0]+cur_delta[1]+vj*gamma/(vj*gamma-cur_delta[1]/lambda[1]-1)*(lambda[0]-lambda_sib_vec[a-1]);
                if (abs(cj) > thresh) {
                  scr_cme[j] = true;
                  num_scr ++;
                }else{
                  scr_cme[j] = false;
                }
              }
            }
            else {
              if (sum_act(beta_cme,couind)==0) {
                thresh = cur_delta[0]+lambda[1]+vj*gamma/(vj*gamma-cur_delta[0]/lambda[0]-1)*(lambda[1]-lambda_cou_vec[b-1]);
                if (abs(cj) > thresh) {
                  scr_cme[j] = true;
                  num_scr ++;
                }else{
                  scr_cme[j] = false;
                }
              }
            }
          }
        }
      }



      //Return trivial solution of \beta=0 when \lambda_s + \lambda_c >= \lambda_max
      // if ( (lambda[0]+lambda[1]) >= lambda_max){
      if ( (a==0) || ( (lambda[0]+lambda[1]) >= lambda_max) ){
        for (int i=0;i<pme;i++){//reset beta
          beta_me[i] = 0.0;
        }
        for (int i=0;i<pcme;i++){
          beta_cme[i] = 0.0;
        }
        if (familyType == "binomial") {
          inter = log(ymean/(1-ymean));
          for (int i=0; i<nn; i++) {
            eta[i] = log(ymean/(1-ymean)); //
            mu = pbinomial(eta[i]) ;
            W[i] = fmax2(mu*(1-mu),0.0001);
            resid[i] = (yy[i]-mu)/W[i];
            dev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
          }
        } else if (familyType == "poisson") {
          inter = log(ymean);
          for (int i=0;i<nn;i++) {
            eta[i] = log(ymean);
            mu = ppoisson(eta[i]) ;
            W[i] = fmax2(mu,0.0001);
            resid[i] = (yy[i]-mu)/W[i];
            if (yy[i]!=0) dev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
            else dev += 2*ymean;
          }
        }
        num_act = 0;
        num_scr = 0;
        for (int i=0;i<pme;i++){//reset active flag
          act_me[i] = true;
          scr_me[i] = true;
          num_act ++;
          num_scr ++;
        }
        for (int i=0;i<pcme;i++){
          act_cme[i] = true;
          scr_cme[i] = true;
          num_act ++;
          num_scr ++;
        }
        // cout << "num_act: " << num_act << endl;
        if ( (lambda[0]+lambda[1]) >= lambda_max){
          goto cycend;
        }
      }

      // // RESET AFTER EACH RUN
      for (int i=0;i<pme;i++){//reset beta
        beta_me[i] = 0.0;
      }
      for (int i=0;i<pcme;i++){
        beta_cme[i] = 0.0;
      }
      if (familyType == "binomial") {
        inter = log(ymean/(1-ymean));
        for (int i=0; i<nn; i++) {
          eta[i] = log(ymean/(1-ymean)); //
          mu = pbinomial(eta[i]) ;
          W[i] = fmax2(mu*(1-mu),0.0001);
          resid[i] = (yy[i]-mu)/W[i];
          dev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
        }
      } else if (familyType == "poisson") {
        inter = log(ymean);
        for (int i=0;i<nn;i++) {
          eta[i] = log(ymean);
          mu = ppoisson(eta[i]) ;
          W[i] = fmax2(mu,0.0001);
          resid[i] = (yy[i]-mu)/W[i];
          if (yy[i]!=0) dev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
          else dev += 2*ymean;
        }
      }

      //Recompute deltas
      fill(delta_sib.begin(),delta_sib.end(),lambda[0]);
      fill(delta_cou.begin(),delta_cou.end(),lambda[1]);
      for (int j=0; j<pme; j++){
        // v = wsqsum(X_me, W, nn, j)/((double)nn);
        string me = Rcpp::as<string>(names_me[j]);
        int delta_ind = findind(effectIndexMap,me);
        delta_sib[delta_ind] = delta_sib[delta_ind] * ( exp( -(tau/lambda[0]) * m_me[j] * mcp(beta_me[j],lambda[0],gamma) ) );
        delta_cou[delta_ind] = delta_cou[delta_ind] * ( exp( -(tau/lambda[1]) * m_me[j] * mcp(beta_me[j],lambda[1],gamma) ) );
      }
      //for (int j=0;j<pme;j++){ //parent effect
      //  for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
      for (int j=0; j>pcme; j++){
        //int cmeind = j*(2*(pme-1))+k;
        //int condind = floor((double)k/2.0);
        //if (condind >= j){
        //  condind ++;
        //}
        string cme = Rcpp::as<string>(names_cme[j]);
        auto parts = splitString(cme);
        int sib_ind = findind(effectIndexMap,parts[0]);
        int cou_ind = findind(effectIndexMap,parts[1]);

        // v = wsqsum(X_cme, W, nn, cmeind)/((double)nn);
        delta_sib[sib_ind] = delta_sib[sib_ind] * (exp(-(tau/lambda[0]) * m_cme[j] * mcp(beta_cme[j],lambda[0],gamma) ));
        delta_cou[cou_ind] = delta_cou[cou_ind] * (exp(-(tau/lambda[1]) * m_cme[j] * mcp(beta_cme[j],lambda[1],gamma) ));
      }

      //Coordinate descent with warm active set resets
      for (int q=0; q<reset; q++){

        //Active set reset for it_warm iterations
        for (int m=0; m<it_warm; m++){
          if (lambda_it && screen_ind && a>0 && b>0){
            chng_flag = coord_des_onerun_str(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy, names_me, names_cme, effectIndexMap,
                                             family, delta_sib, delta_cou, scr_me, scr_cme, inter, beta_me, beta_cme, m_me, m_cme, eta, resid, W, dev);
          } else{
            chng_flag = coord_des_onerun_str(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy, names_me, names_cme, effectIndexMap,
                                             family, delta_sib, delta_cou, act_me, act_cme, inter, beta_me, beta_cme, m_me, m_cme, eta, resid, W, dev);
          }

        }



        //Update active set
        int num_act = 0;
        int num_scr = 0;
        for (int j=0;j<pme;j++){
          if ((abs(beta_me[j])>0.0||(act_vec[j]>0.0))){ //
            act_me[j] = true;
            num_act ++;
            scr_me[j] = true;
            num_scr ++;
          }
          else{
            scr_me[j] = false;
            act_me[j] = false;
          }
        }
        for (int j=0;j<pcme;j++){
          if ((abs(beta_cme[j])>0.0)||(act_vec[j+pme]>0.0)){ //
            act_cme[j] = true;
            num_act ++;
            scr_cme[j] = true;
            num_scr ++;
          }
          else{
            act_cme[j] = false;
            scr_cme[j] = false;
          }
        }

        //Cycle on active set
        it_inner = 0; //inner iteration count
        cont = true; //continue flag
        chng_flag = false; //change flag

        while (cont){
          // cout << it_inner << endl;

          //Increment count and update flags
          it_inner ++;
          if (lambda_it && screen_ind && a>0 && b>0){
            chng_flag = coord_des_onerun_str(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy,  names_me, names_cme, effectIndexMap,
                                             family, delta_sib, delta_cou, scr_me, scr_cme, inter, beta_me, beta_cme, m_me, m_cme, eta, resid, W, dev);
          } else {
            chng_flag = coord_des_onerun_str(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy,names_me, names_cme, effectIndexMap,
                                             family, delta_sib, delta_cou, act_me, act_cme, inter, beta_me, beta_cme, m_me, m_cme, eta, resid, W, dev);
          }

          //Update cont flag for termination
          if ( (it_inner >= it_max_reset)||(!chng_flag) ){
            cont = false;
          }
        }//end while


        // Rcout << accumulate(act_me.begin(),act_me.end(),0) << endl;
        // Rcout << accumulate(act_cme.begin(),act_cme.end(),0) << endl;
        //Update active set
        num_act = 0;
        num_scr = 0;
        for (int j=0;j<pme;j++){
          if ((abs(beta_me[j])>0.0)){ //||(act_vec[j]>0.0)
            act_me[j] = true;
            num_act ++;
            scr_me[j] = true;
            num_scr ++;
          }
          else{
            scr_me[j] = false;
            act_me[j] = false;
          }
        }
        for (int j=0;j<pcme;j++){
          if ((abs(beta_cme[j])>0.0)){ //||(act_vec[j+pme]>0.0)
            act_cme[j] = true;
            num_act ++;
            scr_cme[j] = true;
            num_scr ++;
          }
          else{
            act_cme[j] = false;
            scr_cme[j] = false;
          }
        }


      }

      cycend:


        //Copy into beta_mat, and warm-start next cycle
        int betacount = 0;
      int betanz = 0;
      for (int k=0;k<pme;k++){
        if (abs(beta_me[k])>0.0){
          betanz++;
        }
        beta_mat(betacount,a) = beta_me[k];
        betacount++;
      }
      for (int k=0;k<pcme;k++){
        if (abs(beta_cme[k])>0.0){
          betanz++;
        }
        beta_mat(betacount,a) = beta_cme[k];
        betacount++;
      }
      nz(a,b) = betanz;
      inter_mat(a,b)= inter;
      dev_mat(a,b) = dev;

      //Copy deltas
      for (int k=0;k<pme;k++){
        delta_sib_mat(k,a) = delta_sib[k];
        delta_cou_mat(k,a) = delta_cou[k];
      }

      //Copy residuals
      for (int k=0;k<nn;k++){
        resid_mat(k,a) = resid[k];
      }

      //Copy screening data
      for (int k=0;k<pme;k++){
        if (lambda_it && screen_ind && a>0 && b>0){
          scr_mat(k,a) = scr_me[k];
        } else {
          scr_mat(k,a) = act_me[k];
        }
      }
      for (int k=0;k<pcme;k++){
        if (lambda_it && screen_ind && a>0 && b>0){
          scr_mat(pme+k,a) = scr_cme[k];
        } else {
          scr_mat(pme+k,a) = act_cme[k];
        }
      }

    }//end nlambda.sib (a)

    //Copy into beta_cube
    beta_cube.slice(b) = beta_mat;
    delta_sib_cube.slice(b) = delta_sib_mat;
    delta_cou_cube.slice(b) = delta_cou_mat;
    resid_cube.slice(b) = resid_mat;
    scr_cube.slice(b) = scr_mat;

  }//end nlambda.cou (b)

  //Rescale betas and compute intercepts
  // arma::mat inter_mat(niter_1,niter_2);
  for (int b=0; b<niter_2; b++){ //iterate over cousins...
    for (int a=0; a<niter_1; a++){ //iterate over siblings...
      //Rescale betas to original scale
      inter_mat(a,b) = inter_mat(a,b);
      dev_mat(a,b) = dev_mat(a,b);
      for (int k=0;k<pme;k++){
        // beta_cube(k,a,b) = beta_cube(k,a,b);
        beta_cube(k,a,b) = beta_cube(k,a,b)/XX_me_sl(k);
      }
      for (int k=0;k<pcme;k++){
        // beta_cube(pme+k,a,b) = beta_cube(pme+k,a,b);
        beta_cube(pme+k,a,b) = beta_cube(pme+k,a,b)/XX_cme_sl(k);
      }
    }
  }

  return (List::create(Named("coefficients") = beta_cube,
                       Named("intercept") = inter_mat,
                       Named("residuals") = resid_cube,
                       Named("deviation") = dev_mat,
                       Named("nzero") = nz,
                       Named("lambda_sib") = lambda_sib_vec,
                       Named("lambda_cou") = lambda_cou_vec,
                       Named("act") = scr_cube,
                       Named("gamma") = gamma,
                       Named("tau") = tau
  ));

}


// [[Rcpp::export]]
List cme(NumericMatrix& XX, NumericVector& yy, CharacterVector& family,
         NumericVector& K1,
         NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec,
         NumericVector& gamma_vec, NumericVector& tau_vec,
         NumericVector& XX_sl, NumericVector& beta_vec, NumericVector& act_vec, NumericVector& multiplier,
         double lambda_max, int it_max, int it_warm, int reset, bool screen_ind) {
  // // [[Rcpp::plugins(openmp)]]
  //------------------------------------------------------------
  // XX - Full model matrix including both ME and CME effects (assume normalized)
  // yy - Response vector of length nn
  // family- family of GLM
  // lambda_sib_vec - Vector of sibling penalties (decr. sequence)
  // lambda_cou_vec - Vector of cousin penalties (decr. sequence)
  // tau - Exponential penalty parameter
  // gamma - MC+ non-convex penalty parameter
  // beta_vec - Initial beta value
  // it_max - Maximum iterations for coordinate descent
  //------------------------------------------------------------

  //Variable initialization
  int pp = XX.ncol(); //# of CMEs
  int nn = XX.nrow(); //# of observations
  int nlambdasib = lambda_sib_vec.size();
  int nlambdacou = lambda_cou_vec.size();
  int it_inner = 0;
  int it_max_reset = it_max / reset;
  bool cont = true;
  bool chng_flag = false;
  double v = 0.25;
  double mu = 0;
  int J = K1.size() - 1;
  int pme = J/2; //# of MEs

  // Extract the family type from CharacterVector
  string familyType = Rcpp::as<string>(family[0]);

  //Vectorize model matrices
  vector<double> X(nn*pp); //for ME
  //vector<double> X_cme(nn*pcme); //for CME
  for (int i=0;i<pp;i++){
    for (int j=0;j<nn;j++){
      X[i*nn+j] = XX(j,i);
    }
  }
  // for (int i=0;i<pcme;i++){
  //   for (int j=0;j<nn;j++){
  //     X_cme[i*nn+j] = XX_cme(j,i);
  //   }
  // }
  // store colnames name of ME and CME
  // CharacterVector names_me = colnames(XX_me);
  // CharacterVector names_cme = colnames(XX_cme);

  //Check whether lambda is to be iterated or not
  bool lambda_it;
  int niter_1; //Number to iterate first
  int niter_2; //Number to iterate next
  if (gamma_vec.size()>1){ //Iterate on gamma and tau
    lambda_it = false;
    // niter_1 = gamma_vec.size(); //ch
    // niter_2 = tau_vec.size();
    niter_1 = tau_vec.size();
    niter_2 = gamma_vec.size();
  }
  else{
    lambda_it = true;
    niter_1 = nlambdasib;
    niter_2 = nlambdacou;
  }

  //Containers for beta and active set (alpha)
  arma::cube beta_cube(pp,niter_1,niter_2); //betas to return
  arma::cube delta_sib_cube(J,niter_1,niter_2); //deltas to return
  arma::cube delta_cou_cube(J,niter_1,niter_2); //deltas to return
  arma::mat nz(niter_1,niter_2);
  arma::mat inter_mat(niter_1,niter_2); //intercept to return
  arma::mat dev_mat(niter_1,niter_2); //deviation to return
  arma::mat beta_mat(pp,niter_1);
  arma::mat delta_sib_mat(J,niter_1);
  arma::mat delta_cou_mat(J,niter_1);

  vector<double> beta(pp,0.0); //for MEs
  for (int i=0;i<pp;i++){
    beta[i] = beta_vec[i];
  }
  // vector<double> beta_cme(pcme,0.0); //for CMEs
  // for (int i=0;i<pcme;i++){
  //   beta_cme[i] = beta_vec[pme+i];
  // }
  //double cur_beta = 0.0; //running beta

  //Containers for siblings or cousion family
  //vector<int> sibind(2*(pme-1),0);
  //vector<int> couind(2*(pme-1),0);

  //Set all factors as active to begin
  vector<bool> act(pp,true); //Current active set
  //vector<bool> act_cme(pcme,true);
  vector<bool> scr(pp,true); //Screened active set
  //vector<bool> scr_cme(pcme,true);
  bool kkt_bool;

  //set<string> uniqueEffects;

  // Include main effects
  //for (int i = 0; i < pme; i++) {
  //  uniqueEffects.insert(Rcpp::as<string>(names_me[i]));
  //}

  //// Include parent and child effects from conditional main effects
  //for (int i = 0; i < pcme; i++) {
  //  auto parts = splitString(Rcpp::as<string>(names_cme[i]));
  //  for (const auto& part : parts) {
  //    uniqueEffects.insert(part);
  //  }
  //}

  //unordered_map<string, int> effectIndexMap;
  //int index = 0;
  //for (const auto& effect : uniqueEffects) {
  //  effectIndexMap[effect] = index++;
  //}

  //// Containers for linearized slopes Delta
  //vector<double> delta_sib(uniqueEffects.size(), 0.0); // Linearized penalty for siblings (sib(A), sib(B), ...)
  //vector<double> delta_cou(uniqueEffects.size(), 0.0); // Linearized penalty for cousins (cou(A), cou(B), ...)
  vector<double> mg(J,1);
  //vector<double> m_cme(pcme,1);
  for (int i=0;i<J;i++){
    mg[i] = multiplier[i];
  }
  // for (int i=0;i<pcme;i++){
  //   m_cme[i] = multiplier[pme+i];
  // }

  //Containers for linearized slopes Delta
  vector<double> delta_sib(pme); //Linearized penalty for siblings (sib(A), sib(B), ...)
  vector<double> delta_cou(pme); //Linearized penalty for cousins (cou(A), cou(B), ...)
  NumericVector lambda(2); //Current penalties
  NumericVector cur_delta(2); //Current delta vector
  double gamma;
  double tau;
  lambda[0] = lambda_sib_vec[0];
  lambda[1] = lambda_cou_vec[0];
  gamma = gamma_vec[0];
  tau = tau_vec[0];

  //Update resid (eta)
  vector<double> eta(nn);
  vector<double> W(nn);
  vector<double> resid(nn); //Residual vector
  arma::cube resid_cube(nn,niter_1,niter_2);
  arma::mat resid_mat(nn,niter_1);
  arma::cube scr_cube(pp,niter_1,niter_2); //screening vector
  arma::mat scr_mat(pp,niter_1);

  double nullDev = 0.0;
  double dev = 0.0;
  double inter= 0.0; //for intercept
  double inprod = 0.0; //inner product
  double cj = 0.0;
  double vj = 0.0;
  double thresh = 0.0; //threshold for screening
  int size = 0;
  int num_act = 0;
  int num_scr = 0;

  double ymean = 0.0;
  for (int i=0;i<nn;i++){
    ymean += (1.0/(double)nn)*yy[i];
  }
  if (familyType == "binomial") {
    inter = log(ymean/(1-ymean));
    for (int i=0; i<nn; i++) {
      eta[i] = log(ymean/(1-ymean)); //
      mu = pbinomial(eta[i]) ;
      W[i] = fmax2(mu*(1-mu),0.0001);
      resid[i] = (yy[i]-mu)/W[i];
      nullDev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
      dev = nullDev;
    }
  } else if (familyType == "poisson") {
    inter = log(ymean);
    for (int i=0;i<nn;i++) {
      eta[i] = log(ymean);
      mu = ppoisson(eta[i]) ;
      W[i] = fmax2(mu,0.0001);
      resid[i] = (yy[i]-mu)/W[i];
      if (yy[i]!=0) nullDev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
      else nullDev += 2*ymean;
      dev = nullDev;
    }
  }

  //vector<bool> kkt_v_me(pme,true);
  //vector<bool> kkt_v_cme(pcme,true); //KKT checks

  // Optimize for each penalty combination
  // #pragma omp parallel for
  for (int b=0; b<niter_2; b++){ //iterate over cousins...

    for (int a=0; a<niter_1; a++){ //iterate over siblings...

      // cout << "Tuning ... a = " << a << ", b = " << b << endl;

      //Update iteration variables
      if (lambda_it){
        lambda[0] = lambda_sib_vec[a];
        lambda[1] = lambda_cou_vec[b];
      }
      else{
        // gamma = gamma_vec[a]; //ch
        // tau = tau_vec[b];
        tau = tau_vec[a];
        gamma = gamma_vec[b];
      }

      // //Update screen set
      // if (lambda_it && screen_ind && a>0 && b>0){
      //
      //   int num_scr = 0;
      //   for (int j=0;j<pme;j++){
      //
      //     if (act_me[j] && scr_me[j]){
      //     int size = 2 * (pme - 1);
      //     vector<int> eff(size);
      //     vector<int> sibind(size);
      //     vector<int> couind(size); // Pre-allocate space for couind
      //
      //     // eff index
      //     for(int i = 0; i < size; ++i) {
      //       eff[i] = i + 1 ; // eff index
      //       sibind[i] = i + j*size;  // sibling index
      //     }
      //
      //     // cousin index
      //     for(int jj = 0; jj < size; ++jj) {
      //       int halfCeil = std::ceil(static_cast<double>(eff[jj]) / 2);
      //       if((j+1) > halfCeil) {
      //         couind[jj] = (halfCeil - 1) * size + (eff[jj] % 2 == 0 ? 1 : 0) + (j+1 - 2) * 2;
      //       } else {
      //         couind[jj] = halfCeil * size + (eff[jj] % 2 == 0 ? 1 : 0) + (j+1 - 1) * 2;
      //       }
      //     }
      //
      //
      //     cj = wcrossprod(X_me, resid, W, nn, j);
      //     vj = wsqsum(X_me, W, nn, j)/((double)nn);
      //
      //     if (sum_act(beta_cme,sibind)==0) {
      //       if (sum_act(beta_cme,couind)==0) {
      //         thresh = max(lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[0]-lambda_sib_vec[a-1]),
      //                      lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[1]-lambda_cou_vec[b-1]));
      //         if (abs(cj) > thresh) {
      //           scr_me[j] = true;
      //           num_scr ++;
      //         }else{
      //           scr_me[j] = false;
      //         }
      //       } else if (sum_act(beta_cme,couind)>0) {
      //         thresh = lambda[0]+cur_delta[1]+vj*gamma/(vj*gamma-cur_delta[1]/lambda[1]-1)*(lambda[0]-lambda_sib_vec[a-1]);
      //         if (abs(cj) > thresh) {
      //           scr_me[j] = true;
      //           num_scr ++;
      //         }else{
      //           scr_me[j] = false;
      //         }
      //       }
      //     }
      //     else {
      //       if (sum_act(beta_cme,couind)==0) {
      //         thresh = cur_delta[0]+lambda[1]+vj*gamma/(vj*gamma-cur_delta[0]/lambda[0]-1)*(lambda[1]-lambda_cou_vec[b-1]);
      //         if (abs(cj) > thresh) {
      //           scr_me[j] = true;
      //           num_scr ++;
      //         }else{
      //           scr_me[j] = false;
      //         }
      //       }
      //     }
      //   }
      //   }
      //
      //   for (int j=0;j<pme;j++){ //parent effect
      //    for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
      //   //for (int j=0; j<pcme; j++){
      //
      //     int cmeind = j*(2*(pme-1))+k; //index for CME
      //
      //     if (act_cme[cmeind] && scr_cme[cmeind]){
      //     int condind = floor((double)k/2.0);
      //     if (condind >= j){
      //       condind ++;
      //     }
      //
      //     int size = 2 * (pme - 1);
      //
      //     // cousin index
      //     for(int jj = 0; jj < size; jj++) {
      //
      //       sibind[jj] = jj + j*size;  // sibling index
      //
      //       int halfCeil = ceil(static_cast<double>(jj+1) / 2);
      //       if((condind+1) > halfCeil) {
      //         couind[jj] = (halfCeil - 1) * size + ((jj+1) % 2 == 0 ? 1 : 0) + (condind+1 - 2) * 2;
      //       } else {
      //         couind[jj] = halfCeil * size + ((jj+1) % 2 == 0 ? 1 : 0) + (condind+1 - 1) * 2;
      //       }
      //     }
      //
      //
      //     cj = wcrossprod(X_cme, resid, W, nn, cmeind);
      //     vj = wsqsum(X_cme, W, nn, cmeind)/((double)nn);
      //
      //     if (sum_act(beta_cme,sibind)==0) {
      //       if (sum_act(beta_cme,couind)==0) {
      //         thresh = max(lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[0]-lambda_sib_vec[a-1]),
      //                      lambda[0]+lambda[1]+vj*gamma/(vj*gamma-2)*(lambda[1]-lambda_cou_vec[b-1]));
      //         if (abs(cj) > thresh) {
      //           scr_cme[cmeind] = true;
      //           num_scr ++;
      //         }else{
      //           scr_cme[cmeind] = false;
      //         }
      //       } else if (sum_act(beta_cme,couind)>0) {
      //         thresh = lambda[0]+cur_delta[1]+vj*gamma/(vj*gamma-cur_delta[1]/lambda[1]-1)*(lambda[0]-lambda_sib_vec[a-1]);
      //         if (abs(cj) > thresh) {
      //           scr_cme[cmeind] = true;
      //           num_scr ++;
      //         }else{
      //           scr_cme[cmeind] = false;
      //         }
      //       }
      //     }
      //     else {
      //       if (sum_act(beta_cme,couind)==0) {
      //         thresh = cur_delta[0]+lambda[1]+vj*gamma/(vj*gamma-cur_delta[0]/lambda[0]-1)*(lambda[1]-lambda_cou_vec[b-1]);
      //         if (abs(cj) > thresh) {
      //           scr_cme[cmeind] = true;
      //           num_scr ++;
      //         }else{
      //           scr_cme[cmeind] = false;
      //         }
      //       }
      //     }
      //   }
      //  }
      //   }
      // }



      //Return trivial solution of \beta=0 when \lambda_s + \lambda_c >= \lambda_max
      // if ( (lambda[0]+lambda[1]) >= lambda_max){
      if ( (a==0) || ( (lambda[0]+lambda[1]) >= lambda_max) ){
        for (int i=0;i<pp;i++){//reset beta
          beta[i] = 0.0;
        }
        // for (int i=0;i<pcme;i++){
        //   beta_cme[i] = 0.0;
        // }
        if (familyType == "binomial") {
          inter = log(ymean/(1-ymean));
          for (int i=0; i<nn; i++) {
            eta[i] = log(ymean/(1-ymean)); //
            mu = pbinomial(eta[i]) ;
            W[i] = fmax2(mu*(1-mu),0.0001);
            resid[i] = (yy[i]-mu)/W[i];
            dev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
          }
        } else if (familyType == "poisson") {
          inter = log(ymean);
          for (int i=0;i<nn;i++) {
            eta[i] = log(ymean);
            mu = ppoisson(eta[i]) ;
            W[i] = fmax2(mu,0.0001);
            resid[i] = (yy[i]-mu)/W[i];
            if (yy[i]!=0) dev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
            else dev += 2*ymean;
          }
        }
        num_act = 0;
        num_scr = 0;
        for (int i=0;i<pp;i++){//reset active flag
          act[i] = true;
          scr[i] = true;
          num_act ++;
          num_scr ++;
        }
        // for (int i=0;i<pcme;i++){
        //   act_cme[i] = true;
        //   scr_cme[i] = true;
        //   num_act ++;
        //   num_scr ++;
        // }
        cout << "num_act: " << num_act << endl;
        if ( (lambda[0]+lambda[1]) >= lambda_max){
          goto cycend;
        }
      }

      // // RESET AFTER EACH RUN
      for (int i=0;i<pp;i++){//reset beta
        beta[i] = 0.0;
      }
      // for (int i=0;i<pcme;i++){
      //   beta_cme[i] = 0.0;
      // }
      if (familyType == "binomial") {
        inter = log(ymean/(1-ymean));
        for (int i=0; i<nn; i++) {
          eta[i] = log(ymean/(1-ymean)); //
          mu = pbinomial(eta[i]) ;
          W[i] = fmax2(mu*(1-mu),0.0001);
          resid[i] = (yy[i]-mu)/W[i];
          dev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
        }
      } else if (familyType == "poisson") {
        inter = log(ymean);
        for (int i=0;i<nn;i++) {
          eta[i] = log(ymean);
          mu = ppoisson(eta[i]) ;
          W[i] = fmax2(mu,0.0001);
          resid[i] = (yy[i]-mu)/W[i];
          if (yy[i]!=0) dev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
          else dev += 2*ymean;
        }
      }

      //Recompute deltas
      fill(delta_sib.begin(),delta_sib.end(),lambda[0]); //assigns each element the value lambda[0]
      fill(delta_cou.begin(),delta_cou.end(),lambda[1]); //assigns each element the value lambda[1]
      for (int g=0; g<J; g++) {
        for (int j=K1[g];j<K1[g+1];j++){
          // v = wsqsum(X_me, W, nn, j)/((double)nn);
          // string me = Rcpp::as<string>(names_me[j]);
          // int delta_ind = findind(effectIndexMap,me);
          //delta_sib[j] = delta_sib[j] * ( exp( -(tau/lambda[0]) * m_me[j] * mcp(beta_me[j],lambda[0],gamma) ) );
          //delta_cou[j] = delta_cou[j] * ( exp( -(tau/lambda[1]) * m_me[j] * mcp(beta_me[j],lambda[1],gamma) ) );
          int gind = floor((double)g/2.0); //index for sibling group
          if (g % 2 == 0) {
            delta_sib[gind] = delta_sib[gind] * (exp(-(tau/lambda[0]) * mcp(beta[j],lambda[0],gamma) )) * mg[g];
            delta_cou[gind] = delta_cou[gind] *  mg[g+1]; //(exp(-(tau/lambda[1]) * mcp(beta[j],lambda[1],gamma) )) *
          } else {
            delta_sib[gind] = delta_sib[gind] * mg[g-1]; //(exp(-(tau/lambda[0]) * mcp(beta[j],lambda[0],gamma) )) *
            delta_cou[gind] = delta_cou[gind] * (exp(-(tau/lambda[1]) * mcp(beta[j],lambda[1],gamma) )) * mg[g];
          }
        }
      }

      //Coordinate descent with warm active set resets
      for (int q=0; q<reset; q++){

        //Active set reset for it_warm iterations
        for (int m=0; m<it_warm; m++){
          if (lambda_it && screen_ind && a>0 && b>0){
            chng_flag = coord_des_onerun(pme, nn, K1, lambda, cur_delta, chng_flag, tau, gamma, X, yy,
                                         family, delta_sib, delta_cou, scr, inter, beta, mg, eta, resid, W, dev);
          } else{
            chng_flag = coord_des_onerun(pme, nn, K1, lambda, cur_delta, chng_flag, tau, gamma, X, yy,
                                         family, delta_sib, delta_cou, act, inter, beta, mg, eta, resid, W, dev);
          }

        }



        //Update active set
        int num_act = 0;
        int num_scr = 0;
        for (int j=0;j<pp;j++){
          if ((abs(beta[j])>0.0||(act_vec[j]>0.0))){ //
            act[j] = true;
            num_act ++;
            scr[j] = true;
            num_scr ++;
          }
          else{
            scr[j] = false;
            act[j] = false;
          }
        }

        //Cycle on active set
        it_inner = 0; //inner iteration count
        cont = true; //continue flag
        chng_flag = false; //change flag

        while (cont){
          cout << "it_inner: " << it_inner << endl;

          //Increment count and update flags
          it_inner ++;
          if (lambda_it && screen_ind && a>0 && b>0){
            chng_flag = coord_des_onerun(pme, nn, K1, lambda, cur_delta, chng_flag, tau, gamma, X, yy,
                                         family, delta_sib, delta_cou, scr, inter, beta, mg, eta, resid, W, dev);
          } else {
            chng_flag = coord_des_onerun(pme, nn, K1, lambda, cur_delta, chng_flag, tau, gamma, X, yy,
                                         family, delta_sib, delta_cou, act, inter, beta, mg, eta, resid, W, dev);
          }

          //Update cont flag for termination
          if ( (it_inner >= it_max_reset)||(!chng_flag) ){
            cont = false;
          }
        }//end while


        Rcout << accumulate(act.begin(),act.end(),0) << endl;
        //Update active set
        num_act = 0;
        num_scr = 0;
        for (int j=0;j<pp;j++){
          if ((abs(beta[j])>0.0)){ //||(act_vec[j]>0.0)
            act[j] = true;
            num_act ++;
            scr[j] = true;
            num_scr ++;
          }
          else{
            scr[j] = false;
            act[j] = false;
          }
        }

      }

      cycend:


        //Copy into beta_mat, and warm-start next cycle
        int betacount = 0;
      int betanz = 0;
      for (int k=0;k<pp;k++){
        if (abs(beta[k])>0.0){
          betanz++;
        }
        beta_mat(betacount,a) = beta[k];
        betacount++;
      }
      // for (int k=0;k<pcme;k++){
      //   if (abs(beta_cme[k])>0.0){
      //     betanz++;
      //   }
      //   beta_mat(betacount,a) = beta_cme[k];
      //   betacount++;
      // }
      nz(a,b) = betanz;
      inter_mat(a,b)= inter;
      dev_mat(a,b) = dev;

      //Copy deltas
      for (int k=0;k<J;k++){
        delta_sib_mat(k,a) = delta_sib[k];
        delta_cou_mat(k,a) = delta_cou[k];
      }

      //Copy residuals
      for (int k=0;k<nn;k++){
        resid_mat(k,a) = resid[k];
      }

      //Copy screening data
      for (int k=0;k<pp;k++){
        if (lambda_it && screen_ind && a>0 && b>0){
          scr_mat(k,a) = scr[k];
        } else {
          scr_mat(k,a) = act[k];
        }
      }
      // for (int k=0;k<pcme;k++){
      //   if (lambda_it && screen_ind && a>0 && b>0){
      //     scr_mat(pme+k,a) = scr_cme[k];
      //   } else {
      //     scr_mat(pme+k,a) = act_cme[k];
      //   }
      // }

    }//end nlambda.sib (a)

    //Copy into beta_cube
    beta_cube.slice(b) = beta_mat;
    delta_sib_cube.slice(b) = delta_sib_mat;
    delta_cou_cube.slice(b) = delta_cou_mat;
    resid_cube.slice(b) = resid_mat;
    scr_cube.slice(b) = scr_mat;

  }//end nlambda.cou (b)

  //Rescale betas and compute intercepts
  // arma::mat inter_mat(niter_1,niter_2);
  for (int b=0; b<niter_2; b++){ //iterate over cousins...
    for (int a=0; a<niter_1; a++){ //iterate over siblings...
      //Rescale betas to original scale
      inter_mat(a,b) = inter_mat(a,b);
      dev_mat(a,b) = dev_mat(a,b);
      for (int k=0;k<pp;k++){
        // beta_cube(k,a,b) = beta_cube(k,a,b);
        beta_cube(k,a,b) = beta_cube(k,a,b)/XX_sl(k);
      }
      // for (int k=0;k<pcme;k++){
      //   // beta_cube(pme+k,a,b) = beta_cube(pme+k,a,b);
      //   beta_cube(pme+k,a,b) = beta_cube(pme+k,a,b)/XX_cme_sl(k);
      // }
    }
  }

  return (List::create(Named("coefficients") = beta_cube,
                       Named("intercept") = inter_mat,
                       Named("residuals") = resid_cube,
                       Named("deviation") = dev_mat,
                       Named("nzero") = nz,
                       Named("lambda_sib") = lambda_sib_vec,
                       Named("lambda_cou") = lambda_cou_vec,
                       Named("act") = scr_cube,
                       Named("gamma") = gamma,
                       Named("tau") = tau
  ));

}
