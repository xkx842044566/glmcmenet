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

// Max of x
double max(vector<double>& x, int n) {
  double val = x[0];
  for (int i=1; i<n; i++) {
    if (x[i] > val) val = x[i];
  }
  return(val);
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
      ret = ( abs(inprod)-(delta_r[0]+delta_r[1]) ) / (v - 1.0/gamma*(ratio[0]+ratio[1]));
    }
    else{
      ret = 0.0;
    }
  }
  else if (abs(inprod) < v*lambda_r[1]*gamma){
    if (abs(inprod) > (delta_r[1])){
      ret = ( abs(inprod)-(delta_r[1]) ) / (v - 1.0/gamma*(ratio[1]) );
    }
    else{
      ret = 0.0;
    }
  }
  else{
    ret = abs(inprod)/v;
  }
  return (sgn*ret);

}

double s_me_gaussian(double inprod, NumericVector& lambda, double gamma, NumericVector& delta){

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

  if (abs(inprod) < (lambda_r[0]*gamma + delta_r[1]*(1-lambda_r[0]/lambda_r[1]))  ){
    if (abs(inprod) > (delta_r[0]+delta_r[1]) ){
      ret = ( abs(inprod)-(delta_r[0]+delta_r[1]) ) / (1.0 - 1.0/gamma*(ratio[0]+ratio[1]) );
    }
    else{
      ret = 0.0;
    }
  }
  else if (abs(inprod) < lambda_r[1]*gamma){
    if (abs(inprod) > (delta_r[1])){
      ret = ( abs(inprod)-(delta_r[1]) ) / (1.0 - 1.0/gamma*(ratio[1]) );
    }
    else{
      ret = 0.0;
    }
  }
  else{
    ret = abs(inprod);
  }

  return (sgn*ret);

}

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

//KKT condition check
bool kkt(double inprod,  NumericVector cur_delta){
  //Checks KKT condition for \beta=0.0
  bool ret;
  double lb = -inprod - cur_delta[0] - cur_delta[1];
  double ub = -inprod + cur_delta[0] + cur_delta[1];
  if ((0.0 >= lb)&&(0.0 <= ub)){
    ret = true; //kkt satisfied
  }
  else{
    ret = false; //kkt not satisfied
  }
  return(ret);
}

//screen rule condition check
bool screen_rule_me(double inprod, double j, double a, double b,
                    NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec,
                    double gamma_v, vector<double>& delta_sib, vector<double>& delta_cou,
                    vector<double>& m_me, vector<double>& m_sib, vector<double>& m_cou){

  double thresh=0;
  double thresh1=0;
  double thresh2=0;
  double thresh3=0;
  bool active=true;
  // if no active effect in sibling for setting (lambda_sib_l-1, lambda_cou_m)
  if(abs(lambda_sib_vec[a-1]-delta_sib[j])<1e-05){
    thresh1 = lambda_sib_vec[a] + delta_cou[j] + gamma_v*m_sib[j]/(gamma_v/m_me[j]-m_sib[j]-m_cou[j]*delta_cou[j]/lambda_cou_vec[b])*(lambda_sib_vec[a]-lambda_sib_vec[a-1]);
  }
  if(abs(lambda_cou_vec[b-1]-delta_cou[j])<1e-05){
    thresh2 = delta_sib[j] + lambda_cou_vec[b] + gamma_v*m_cou[j]/(gamma_v/m_me[j]-m_sib[j]*delta_sib[j]/lambda_sib_vec[a]-m_cou[j])*(lambda_cou_vec[b]-lambda_cou_vec[b-1]);
  }
  if(abs(lambda_sib_vec[a-1]-delta_sib[j])<1e-05 && abs(lambda_cou_vec[b-1]-delta_cou[j])<1e-05){
    thresh3 = max(lambda_sib_vec[a]+lambda_cou_vec[b]+gamma_v*m_sib[j]/(gamma_v/m_me[j]-m_sib[j]-m_cou[j])*(lambda_sib_vec[a]-lambda_sib_vec[a-1]),
                  lambda_sib_vec[a]+lambda_cou_vec[b]+gamma_v*m_cou[j]/(gamma_v/m_me[j]-m_sib[j]-m_cou[j])*(lambda_cou_vec[b]-lambda_cou_vec[b-1]));
  }
  thresh = max(max(thresh1,thresh2),thresh3);
  if (abs(inprod) < thresh) {
    active = false;
  } else {
    active = true;
  }

  return(active);
}


bool screen_rule_cme(double inprod, double j, double k, double cmeind, double condind, double a, double b,
                     NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec,
                     double gamma_v, vector<double>& delta_sib, vector<double>& delta_cou,
                     vector<double>& m_cme,vector<double>& m_sib, vector<double>& m_cou){

  double thresh=0;
  double thresh1=0;
  double thresh2=0;
  double thresh3=0;
  bool active=true;
  // if no active effect in sibling for setting (lambda_sib_l-1, lambda_cou_m)
  if(abs(lambda_sib_vec[a-1]-delta_sib[j])<1e-04){
    thresh1 = lambda_sib_vec[a] + delta_cou[condind] + gamma_v*m_sib[j]/(gamma_v/m_cme[cmeind]-m_sib[j]-m_cou[condind]*delta_cou[condind]/lambda_cou_vec[b])*(lambda_sib_vec[a]-lambda_sib_vec[a-1]);
  }
  //cout << "diff2:" << abs(lambda_cou_vec[b-1]-delta_cou[j]) << endl;
  if(abs(lambda_cou_vec[b-1]-delta_cou[j])<1e-04){
    thresh2 = delta_sib[j] + lambda_cou_vec[b] + gamma_v*m_cou[condind]/(gamma_v/m_cme[cmeind]-m_sib[j]*delta_sib[j]/lambda_sib_vec[a]-m_cou[condind])*(lambda_cou_vec[b]-lambda_cou_vec[b-1]);
  }
  if(abs(lambda_sib_vec[a-1]-delta_sib[j])<1e-04 && abs(lambda_cou_vec[b-1]-delta_cou[condind])<1e-04){
    thresh3 = max(lambda_sib_vec[a]+lambda_cou_vec[b]+gamma_v*m_sib[j]/(gamma_v/m_cme[cmeind]-m_sib[j]-m_cou[condind])*(lambda_sib_vec[a]-lambda_sib_vec[a-1]),
                  lambda_sib_vec[a]+lambda_cou_vec[b]+gamma_v*m_cou[condind]/(gamma_v/m_cme[cmeind]-m_sib[j]-m_cou[condind])*(lambda_cou_vec[b]-lambda_cou_vec[b-1]));
  }
  thresh = max(max(thresh1,thresh2),thresh3);
  if (abs(inprod) < thresh) {
    active = false;
  } else {
    active = true;
  }

  return(active);
}



bool coord_des_onerun(int pme, int pcme, int nn, NumericVector& lambda, NumericVector& cur_delta,
                      bool dummy, double tau, double gamma,
                      vector<double>& X_me, vector<double>& X_cme, NumericVector& yy,
                      CharacterVector& family,
                      vector<double>& delta_sib, vector<double>& delta_cou,
                      vector<bool>& act_me, vector<bool>& act_cme, double& inter,
                      vector<double>& beta_me, vector<double>& beta_cme,
                      vector<double>& m_me, vector<double>& m_cme,vector<double>& m_sib, vector<double>& m_cou,
                      vector<double>& eta, vector<double>& resid, vector<double>& W, double& dev){
  bool chng_flag = false;
  double cur_beta = 0.0;
  double cur_inter = 0.0;
  double xwr = 0.0;
  double xwx = 0.0;
  double inprod = 0.0;
  double v = 0.0;
  double mu = 0.0;
  NumericVector lambda_adp(2);
  //double dev = 0.0;

  // Extract the family type from CharacterVector
  string familyType = Rcpp::as<string>(family[0]);



  if (familyType == "binomial") {
    v = 0.25;
    for (int i=0; i<nn; i++) {
      mu = pbinomial(eta[i]);
      W[i] = fmax2(mu*(1-mu), 0.0001);
      resid[i] = (yy[i] - mu)/W[i];
      if (yy[i]==1) dev = dev - log(mu);
      if (yy[i]==0) dev = dev - log(1-mu);
    }
  } else if (familyType == "poisson") {
    v = exp(max(eta, nn));
    for (int i=0; i<nn; i++) {
      mu = exp(eta[i]);
      W[i] = mu;
      resid[i] = (yy[i] - mu)/W[i];
      if (yy[i]!=0) dev += yy[i]*log(yy[i]/mu);
    }
  }

  //Update intercept
  cur_inter= inter;

  xwr = crossprod(W, resid, nn, 0);
  xwx = sum(W,nn);
  inter = xwr/xwx + cur_inter;
  for (int i=0; i<nn; i++) {
    resid[i] -= inter - cur_inter;
    eta[i] += inter - cur_inter;
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
      inprod = xwr/((double)nn)+v*beta_me[j]; //checked to pod from update eqn (mod from above eqn)

      //Update cur_delta
      cur_delta[0] = delta_sib[j]*m_me[j];//*m_sib[j]*m_me[j];
      cur_delta[1] = delta_cou[j]*m_me[j];//*m_cou[j]*m_me[j];

      //adaptive lambda
      lambda_adp[0] = lambda[0]*m_sib[j];
      lambda_adp[1] = lambda[1]*m_cou[j];

      //adaptive Perform ME thresholding
      beta_me[j] = s_me(inprod,v,lambda_adp,gamma,cur_delta);

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
        }

        //Update adaptive deltas
        double offset_sib = mcp(beta_me[j],lambda_adp[0],gamma)-mcp(cur_beta,lambda_adp[0],gamma); // new - old
        double offset_cou = mcp(beta_me[j],lambda_adp[1],gamma)-mcp(cur_beta,lambda_adp[1],gamma);
        delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda_adp[0]) * m_me[j] * offset_sib ));
        delta_cou[j] = delta_cou[j] * (exp(-(tau/lambda_adp[1]) * m_me[j] * offset_cou ));

        //Update flag
        chng_flag = true;

      }
    }
  }

  cur_beta = 0.0;
  inprod = 0.0;

  //CD for CME effects
  for (int j=0;j<pme;j++){ //parent effect
    for (int k=0;k<(2*(pme-1));k++){ //conditioned effect

      int cmeind = j*(2*(pme-1))+k; //index for CME

      if (act_cme[cmeind]){
        int condind = floor((double)k/2.0);
        if (condind >= j){
          condind ++;
        }

        cur_beta = beta_cme[cmeind]; //beta_cme ordered by parent effects, then condition

        //Update cur_delta
        cur_delta[0] = delta_sib[j] * m_cme[cmeind];// * m_sib[j] * m_cme[cmeind];
        cur_delta[1] = delta_cou[condind] * m_cme[cmeind];// * m_cou[condind] * m_cme[cmeind];

        // Updata covariates
        xwr = wcrossprod(X_cme, resid, W, nn, cmeind);
        xwx = wsqsum(X_cme, W, nn, cmeind);
        v = xwx/((double)nn);
        inprod = xwr/((double)nn)+v*beta_cme[cmeind]; //checked to pod from update eqn (mod from above eqn)

        //adaptive lambda
        lambda_adp[0] = lambda[0]*m_sib[j];
        lambda_adp[1] = lambda[1]*m_cou[condind];

        //adaptive Perform CME thresholding
        beta_cme[cmeind] = s_me(inprod,v,lambda_adp,gamma,cur_delta);

        //Update eta, mu, weight and delta
        if (!dbleq(beta_cme[cmeind],cur_beta)){ // if beta changed...

          //Update resid eta, mu, weight
          for (int ll=0;ll<nn;ll++){
            resid[ll] -= X_cme[cmeind*nn+ll]*(beta_cme[cmeind]-cur_beta);
            eta[ll] += X_cme[cmeind*nn+ll]*(beta_cme[cmeind]-cur_beta);
            if (familyType == "binomial"){
              mu = pbinomial(eta[ll]) ;
              W[ll] = fmax2(mu*(1-mu),0.0001);
            }else if(familyType == "poisson"){
              mu = ppoisson(eta[ll]) ;
              W[ll] = fmax2(mu,0.0001);
            }
          }

           //Update deltas
          double offset_sib = mcp(beta_cme[cmeind],lambda_adp[0],gamma)-mcp(cur_beta,lambda_adp[0],gamma); // new - old
          double offset_cou = mcp(beta_cme[cmeind],lambda_adp[1],gamma)-mcp(cur_beta,lambda_adp[1],gamma); // new - old
          delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda_adp[0]) * m_cme[cmeind] * offset_sib )); // update delta for siblings
          delta_cou[condind] = delta_cou[condind] * (exp(-(tau/lambda_adp[1]) * m_cme[cmeind] * offset_cou )); // update delta for cousins

          //Update flag
          chng_flag = true;

        }



        //Reduce A|B+ and A|B- to A
        //
        if (abs(beta_cme[cmeind]) > 0.0){ //if current CME is active
          if (k % 2 == 0){ //cme is .|.+
            if (abs(beta_cme[cmeind+1]) > 0.0){ //if cme .|.- is also in model...

              double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

              if ( abs(beta_cme[cmeind]) > abs(beta_cme[cmeind+1]) ){// if abs(.|.+) > abs(.|.-)
                chg = beta_cme[cmeind+1]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind+1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind] -= chg; // update larger CME
                beta_cme[cmeind+1] = 0.0; // remove smaller CME
              }else{// if abs(.|.+) < abs(.|.-)
                chg = beta_cme[cmeind]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind+1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind+1] -= chg; // update larger CME
                beta_cme[cmeind] = 0.0; // remove smaller CME
              }

              //adaptive Update deltas and flag
              double offset_sib = mcp(beta_me[j],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_me,lambda[0]*m_sib[j],gamma); // new - old (for me)
              double offset_cou = mcp(beta_me[j],lambda[1]*m_cou[j],gamma)-mcp(cur_beta_me,lambda[1]*m_cou[j],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_me[j] * offset_sib ));
              delta_cou[j] = delta_cou[j] * (exp(-(tau/(lambda[1]*m_cou[j])) * m_me[j] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme1,lambda[0]*m_sib[j],gamma); // new - old (for .|.+)
              offset_cou = mcp(beta_cme[cmeind],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme1,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind+1],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme2,lambda[0]*m_sib[j],gamma); // new - old (for .|.-)
              offset_cou = mcp(beta_cme[cmeind+1],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme2,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind+1] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind+1] * offset_cou ));

              //residuals shouldn't change

            }
          }else{ //cme is .|.-
            if (abs(beta_cme[cmeind-1]) > 0.0){ //if cme .|.+ is also in model...

              double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

              if ( abs(beta_cme[cmeind]) > abs(beta_cme[cmeind-1]) ){// if abs(.|.+) < abs(.|.-)
                chg = beta_cme[cmeind-1]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind-1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind] -= chg; // update larger CME
                beta_cme[cmeind-1] = 0.0; // remove smaller CME
              }else{// if abs(.|.+) > abs(.|.-)
                chg = beta_cme[cmeind]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind-1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind-1] -= chg; // update larger CME
                beta_cme[cmeind] = 0.0; // remove smaller CME
              }

               //Update deltas and flag
              double offset_sib = mcp(beta_me[j],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_me,lambda[0]*m_sib[j],gamma); // new - old (for me)
              double offset_cou = mcp(beta_me[j],lambda[1]*m_cou[j],gamma)-mcp(cur_beta_me,lambda[1]*m_cou[j],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_me[j] * offset_sib ));
              delta_cou[j] = delta_cou[j] * (exp(-(tau/(lambda[1]*m_cou[j])) * m_me[j] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme1,lambda[0]*m_sib[j],gamma); // new - old (for .|.+)
              offset_cou = mcp(beta_cme[cmeind],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme1,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind-1],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme2,lambda[0]*m_sib[j],gamma); // new - old (for .|.-)
              offset_cou = mcp(beta_cme[cmeind-1],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme2,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind-1] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind-1] * offset_cou ));

              //residuals shouldn't change

            }
          }
        }
      }
    }
  }


  return (chng_flag);

}



bool coord_des_onerun_gaussian(int pme, int pcme, int nn, NumericVector& lambda, NumericVector& cur_delta,
                               bool dummy, double tau, double gamma,
                               vector<double>& X_me, vector<double>& X_cme, NumericVector& yy,
                               vector<double>& delta_sib, vector<double>& delta_cou,
                               vector<bool>& act_me, vector<bool>& act_cme, //double& inter,
                               vector<double>& beta_me, vector<double>& beta_cme,
                               vector<double>& m_me, vector<double>& m_cme,vector<double>& m_sib, vector<double>& m_cou,
                               vector<double>& resid){
  bool chng_flag = false;
  double cur_beta = 0.0;
  //double cur_inter = 0.0;
  double inprod = 0.0;
  NumericVector lambda_adp(2);


  //cur_inter= inter;

  for (int j=0;j<pme;j++){
    //Only update if active
    if (act_me[j]){

      cur_beta = beta_me[j];

      //Compute inner product
      inprod = 0.0;
      for (int k=0;k<nn;k++){
        inprod += (resid[k]*X_me[j*nn+k]);
      }
      // inprod = inprod/((double)nn)+beta_me[j];
      inprod = inprod/((double)nn)+(((double)nn)-1)/((double)nn)*beta_me[j]; //checked to pod from update eqn (mod from above eqn)

      //Update cur_delta
      cur_delta[0] = delta_sib[j]*m_me[j];//*m_sib[j]*m_me[j];
      cur_delta[1] = delta_cou[j]*m_me[j];//*m_cou[j]*m_me[j];

      //adaptive lambda
      lambda_adp[0] = lambda[0]*m_sib[j];
      lambda_adp[1] = lambda[1]*m_cou[j];

      //Perform ME thresholding
      beta_me[j] = s_me_gaussian(inprod,lambda_adp,gamma,cur_delta); //s_me_gaussian(inprod,lambda,gamma,cur_delta);

      // Update eta, mu, weight and delta
      if (!dbleq(beta_me[j],cur_beta)){ // if beta changed...

        //Update resid eta, mu, weight
        for (int k=0;k<nn;k++){
          resid[k] = resid[k] - X_me[j*nn+k]*(beta_me[j]-cur_beta);
        }

        //Update deltas
        double offset_sib = mcp(beta_me[j],lambda_adp[0],gamma)-mcp(cur_beta,lambda_adp[0],gamma); // new - old
        double offset_cou = mcp(beta_me[j],lambda_adp[1],gamma)-mcp(cur_beta,lambda_adp[1],gamma);
        delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda_adp[0]) * m_me[j] * offset_sib ));
        delta_cou[j] = delta_cou[j] * (exp(-(tau/lambda_adp[1]) * m_me[j] * offset_cou ));

        //Update flag
        chng_flag = true;
      }
    }
  }

  cur_beta = 0.0;
  inprod = 0.0;
  //CD for covariates in group g
  for (int j=0;j<pme;j++){ //parent effect
    for (int k=0;k<(2*(pme-1));k++){ //conditioned effect

      int cmeind = j*(2*(pme-1))+k; //index for CME

      //Only update if active
      if (act_cme[cmeind]){
        //int condind = 0; //index for condition
        int condind = floor((double)k/2.0);
        if (condind >= j){
          condind ++;
        }

        cur_beta = beta_cme[cmeind];

        //Update cur_delta
        cur_delta[0] = delta_sib[j] * m_cme[cmeind]; //* m_sib[j] * m_cme[cmeind];
        cur_delta[1] = delta_cou[condind] * m_cme[cmeind]; //* m_cou[condind] * m_cme[cmeind];

        //Compute inner product
        inprod = 0.0;
        for (int l=0;l<nn;l++){
          inprod += (resid[l]*X_cme[cmeind*nn+l]);
        }
        // inprod = inprod/((double)nn)+beta_cme[cmeind];
        inprod = inprod/((double)nn)+(((double)nn)-1)/((double)nn)*beta_cme[cmeind];

        //adaptive lambda
        lambda_adp[0] = lambda[0]*m_sib[j];
        lambda_adp[1] = lambda[1]*m_cou[condind];

        //Perform ME thresholding
        beta_cme[cmeind] = s_me_gaussian(inprod,lambda_adp,gamma,cur_delta); //s_me_gaussian(inprod,lambda,gamma,cur_delta);

        // Update eta, mu, weight and delta
        if (!dbleq(beta_cme[cmeind],cur_beta)){ // if beta changed...

          //Update residual vector
          for (int ll=0;ll<nn;ll++){
            resid[ll] = resid[ll] - X_cme[cmeind*nn+ll]*(beta_cme[cmeind]-cur_beta);
          }

          //Update deltas
          double offset_sib = mcp(beta_cme[cmeind],lambda_adp[0],gamma)-mcp(cur_beta,lambda_adp[0],gamma); // new - old
          double offset_cou = mcp(beta_cme[cmeind],lambda_adp[1],gamma)-mcp(cur_beta,lambda_adp[1],gamma); // new - old
          delta_sib[j] = delta_sib[j] * (exp(-(tau/lambda_adp[0]) * m_cme[cmeind] * offset_sib )); // update delta for siblings
          delta_cou[condind] = delta_cou[condind] * (exp(-(tau/lambda_adp[1]) * m_cme[cmeind] * offset_cou )); // update delta for cousins

          //Update flag
          chng_flag = true;

        }

        //Reduce A|B+ and A|B- to A
        //
        if (abs(beta_cme[cmeind]) > 0.0){ //if current CME is active
          if (k % 2 == 0){ //cme is .|.+
            if (abs(beta_cme[cmeind+1]) > 0.0){ //if cme .|.- is also in model...

              double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

              if ( abs(beta_cme[cmeind]) > abs(beta_cme[cmeind+1]) ){// if abs(.|.+) > abs(.|.-)
                chg = beta_cme[cmeind+1]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind+1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind] -= chg; // update larger CME
                beta_cme[cmeind+1] = 0.0; // remove smaller CME
              }else{// if abs(.|.+) < abs(.|.-)
                chg = beta_cme[cmeind]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind+1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind+1] -= chg; // update larger CME
                beta_cme[cmeind] = 0.0; // remove smaller CME
              }

              //adaptive Update deltas and flag
              double offset_sib = mcp(beta_me[j],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_me,lambda[0]*m_sib[j],gamma); // new - old (for me)
              double offset_cou = mcp(beta_me[j],lambda[1]*m_cou[j],gamma)-mcp(cur_beta_me,lambda[1]*m_cou[j],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_me[j] * offset_sib ));
              delta_cou[j] = delta_cou[j] * (exp(-(tau/(lambda[1]*m_cou[j])) * m_me[j] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme1,lambda[0]*m_sib[j],gamma); // new - old (for .|.+)
              offset_cou = mcp(beta_cme[cmeind],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme1,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind+1],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme2,lambda[0]*m_sib[j],gamma); // new - old (for .|.-)
              offset_cou = mcp(beta_cme[cmeind+1],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme2,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind+1] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind+1] * offset_cou ));

              //residuals shouldn't change

            }
          }else{ //cme is .|.-
            if (abs(beta_cme[cmeind-1]) > 0.0){ //if cme .|.+ is also in model...

              double chg, cur_beta_me, cur_beta_cme1, cur_beta_cme2;

              if ( abs(beta_cme[cmeind]) > abs(beta_cme[cmeind-1]) ){// if abs(.|.+) < abs(.|.-)
                chg = beta_cme[cmeind-1]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind-1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind] -= chg; // update larger CME
                beta_cme[cmeind-1] = 0.0; // remove smaller CME
              }else{// if abs(.|.+) > abs(.|.-)
                chg = beta_cme[cmeind]; // change
                cur_beta_me = beta_me[j]; // current beta me
                cur_beta_cme1 = beta_cme[cmeind]; // current beta cme 1
                cur_beta_cme2 = beta_cme[cmeind-1]; // current beta cme 2
                beta_me[j] += chg; // update ME with smaller CME
                beta_cme[cmeind-1] -= chg; // update larger CME
                beta_cme[cmeind] = 0.0; // remove smaller CME
              }


              //Update deltas and flag
              double offset_sib = mcp(beta_me[j],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_me,lambda[0]*m_sib[j],gamma); // new - old (for me)
              double offset_cou = mcp(beta_me[j],lambda[1]*m_cou[j],gamma)-mcp(cur_beta_me,lambda[1]*m_cou[j],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_me[j] * offset_sib ));
              delta_cou[j] = delta_cou[j] * (exp(-(tau/(lambda[1]*m_cou[j])) * m_me[j] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme1,lambda[0]*m_sib[j],gamma); // new - old (for .|.+)
              offset_cou = mcp(beta_cme[cmeind],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme1,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind] * offset_cou ));

              offset_sib = mcp(beta_cme[cmeind-1],lambda[0]*m_sib[j],gamma)-mcp(cur_beta_cme2,lambda[0]*m_sib[j],gamma); // new - old (for .|.-)
              offset_cou = mcp(beta_cme[cmeind-1],lambda[1]*m_cou[condind],gamma)-mcp(cur_beta_cme2,lambda[1]*m_cou[condind],gamma);
              delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind-1] * offset_sib ));
              delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind-1] * offset_cou ));

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
List cme(NumericMatrix& XX_me, NumericMatrix& XX_cme, NumericVector& yy, CharacterVector& family,
         NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec,
         NumericVector& gamma_vec, NumericVector& tau_vec,
         NumericVector& XX_me_sl, NumericVector& XX_cme_sl, NumericVector& beta_vec, NumericVector& act_vec,
         NumericVector& multiplier, NumericVector& multiplier_g,
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
  double inprod = 0.0;
  double v = 0.0; (void)v;
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

  //Check whether lambda is to be iterated or not
  bool lambda_it;
  int niter_1; //Number to iterate first
  int niter_2; //Number to iterate next
  if (gamma_vec.size()>1){ //Iterate on gamma and tau
    lambda_it = false;
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

  //Set all factors as active to begin
  vector<bool> act_me(pme,true); //Current active set
  vector<bool> act_cme(pcme,true);
  vector<bool> scr_me(pme,true); //Screened active set
  vector<bool> scr_cme(pcme,true);
  bool kkt_bool;

  vector<double> m_me(pme,1);
  vector<double> m_cme(pcme,1);
  for (int i=0;i<pme;i++){
    m_me[i] = multiplier[i];
  }
  for (int i=0;i<pcme;i++){
    m_cme[i] = multiplier[pme+i];
  }

  vector<double> m_sib(pme,1);
  vector<double> m_cou(pme,1);
  for (int i=0;i<pme;i++){
    m_sib[i] = multiplier_g[2*i];
  }
  for (int i=0;i<pme;i++){
    m_cou[i] = multiplier_g[2*i+1];
  }


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
  arma::cube scr_cube(pme+pcme,niter_1,niter_2); //screening vector
  arma::mat scr_mat(pme+pcme,niter_1);

  double nullDev = 0.0;
  double dev = 0.0;
  double inter= 0.0; //for intercept
  double cj = 0.0;
  double vj = 0.0;
  double thresh = 0.0; (void)thresh; //threshold for screening
  int num_act = 0; (void)num_act;
  int num_scr = 0; (void)num_scr;

  double ymean = 0.0;
  for (int i=0;i<nn;i++){
    ymean += (1.0/(double)nn)*yy[i];
  }
  if (familyType == "binomial") {
    inter = log(ymean/(1-ymean));
    for (int i=0; i<nn; i++) {
      eta[i] = log(ymean/(1-ymean)); //
      mu = ymean ;
      W[i] = fmax2(mu*(1-mu),0.0001);
      resid[i] = (yy[i]-mu)/W[i];
      nullDev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
      dev = nullDev;
    }
  } else if (familyType == "poisson") {
    inter = log(ymean);
    for (int i=0;i<nn;i++) {
      eta[i] = log(ymean);
      mu = ymean ;
      W[i] = mu;
      resid[i] = (yy[i]-mu)/W[i];
      if (yy[i]!=0) nullDev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
      else nullDev += 2*ymean;
      dev = nullDev;
    }
  }

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
        tau = tau_vec[a];
        gamma = gamma_vec[b];
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
            mu = ymean ;
            W[i] = fmax2(mu*(1-mu),0.0001);
            resid[i] = (yy[i]-mu)/W[i];
            dev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
          }
        } else if (familyType == "poisson") {
          inter = log(ymean);
          for (int i=0;i<nn;i++) {
            eta[i] = log(ymean);
            mu = ymean ;
            W[i] = mu;
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

      // //Update strong set
      if (lambda_it && screen_ind && a>0 && b>0){

        for (int j=0;j<pme;j++){

          cj = wcrossprod(X_me, resid, W, nn, j)/((double)nn);
          vj = wsqsum(X_me, W, nn, j)/((double)nn);

          scr_me[j] = screen_rule_me(cj, j, a, b, lambda_sib_vec, lambda_cou_vec, gamma*vj, delta_sib, delta_cou,
                                     m_me, m_sib, m_cou);
          if(!scr_me[j]){num_scr --;}
        }

        for (int j=0;j<pme;j++){ //parent effect
          for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
            //for (int j=0; j<pcme; j++){

            int cmeind = j*(2*(pme-1))+k; //index for CME

            int condind = floor((double)k/2.0);
            if (condind >= j){
              condind ++;
            }

            cj = wcrossprod(X_cme, resid, W, nn, cmeind)/((double)nn);
            vj = wsqsum(X_cme, W, nn, cmeind)/((double)nn);

            scr_cme[cmeind] = screen_rule_cme(cj, j, k, cmeind, condind, a, b, lambda_sib_vec, lambda_cou_vec, gamma*vj, delta_sib, delta_cou,
                                              m_cme, m_sib, m_cou);
            if(!scr_cme[cmeind]){num_scr --;}
          }
        }
      }

      //cout << "strong set:" << num_scr << endl;

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
          mu = ymean ;
          W[i] = fmax2(mu*(1-mu),0.0001);
          resid[i] = (yy[i]-mu)/W[i];
          dev -= 2*yy[i]*log(ymean) + 2*(1-yy[i])*log(1-ymean);
        }
      } else if (familyType == "poisson") {
        inter = log(ymean);
        for (int i=0;i<nn;i++) {
          eta[i] = log(ymean);
          mu = ymean ;
          W[i] = mu;
          resid[i] = (yy[i]-mu)/W[i];
          if (yy[i]!=0) dev += 2*(yy[i]*log(yy[i]/ymean) + ymean - yy[i]);
          else dev += 2*ymean;
        }
      }

      num_scr = 0;
      for (int i=0;i<pme;i++){//reset active flag
        scr_me[i] = true;
        num_scr ++;
      }
      for (int i=0;i<pcme;i++){
        scr_cme[i] = true;
        num_scr ++;
      }


      //adaptive Recompute deltas
      for (int j=0; j<pme; j++){
        delta_sib[j] = lambda[0] * m_sib[j];
        delta_cou[j] = lambda[1] * m_cou[j];
      }
      for (int j=0; j<pme; j++){
        delta_sib[j] = delta_sib[j] * ( exp( -(tau/(lambda[0]*m_sib[j])) * m_me[j] * mcp(beta_me[j],lambda[0]*m_sib[j],gamma) ) );
        delta_cou[j] = delta_cou[j] * ( exp( -(tau/(lambda[1]*m_cou[j])) * m_me[j] * mcp(beta_me[j],lambda[1]*m_cou[j],gamma) ) );
      }
      for (int j=0;j<pme;j++){ //parent effect
        for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
          int cmeind = j*(2*(pme-1))+k;
          int condind = floor((double)k/2.0);
          if (condind >= j){
            condind ++;
          }
          delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind] * mcp(beta_cme[cmeind],lambda[0]*m_sib[j],gamma) ));
          delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind] * mcp(beta_cme[cmeind],lambda[1]*m_cou[condind],gamma) ));
        }
      }

      //Coordinate descent with warm active set resets
      for (int q=0; q<reset; q++){

        //Active set reset for it_warm iterations
        for (int m=0; m<it_warm; m++){
           chng_flag = coord_des_onerun(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy,
                                       family, delta_sib, delta_cou, act_me, act_cme, inter, beta_me, beta_cme, m_me, m_cme, m_sib, m_cou, eta, resid, W, dev);

        }



        //Update active set
        int num_act = 0; (void)num_act;
        for (int j=0;j<pme;j++){
          if ((abs(beta_me[j])>0.0||(act_vec[j]>0.0))){ //
            act_me[j] = true;
            num_act ++;
          }
          else{
            act_me[j] = false;
          }
        }
        for (int j=0;j<pcme;j++){
          if ((abs(beta_cme[j])>0.0)||(act_vec[j+pme]>0.0)){ //
            act_cme[j] = true;
            num_act ++;
          }
          else{
            act_cme[j] = false;
          }
        }

        //cout << "warm act:" << num_act << endl;

        //Cycle on active set
        it_inner = 0; //inner iteration count

        while (it_inner < it_max_reset){

          //while (it_inner < it_max_reset){

          cont = true; //continue flag
          chng_flag = false; //change flag

          while (cont){
            // cout << it_inner << endl;

            //Increment count and update flags on active sets
            it_inner ++;
            chng_flag = coord_des_onerun(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy,
                                         family, delta_sib, delta_cou, act_me, act_cme, inter, beta_me, beta_cme, m_me, m_cme, m_sib, m_cou, eta, resid, W, dev);

            //Update cont flag for termination
            if ( (it_inner >= it_max_reset)||(!chng_flag)){
              cont = false;
            }
          }//end while

          //cout << "it_inner: " << it_inner << endl;

          //check violation in complementary active set\strong set, update active set
          int violations = 0;
          for (int j=0; j<pme; j++) {
            if (act_cme[j]==false && scr_me[j]==false) {
              inprod = wcrossprod(X_me, resid, W, nn, j)/((double)nn); //checked to pod from update eqn (mod from above eqn)
              kkt_bool = kkt(inprod, cur_delta);
              if (!kkt_bool) {
                act_me[j] = scr_me[j] = true;
                violations++;
              }
            }
          }
          for (int j=0; j<pcme; j++) {
            if (act_cme[j]==false && scr_cme[j]==false) {
              inprod = wcrossprod(X_cme, resid, W, nn, j)/((double)nn); //checked to pod from update eqn (mod from above eqn)
              kkt_bool = kkt(inprod, cur_delta);
              if (!kkt_bool) {
                act_cme[j] = scr_cme[j] = true;
                violations++;
              }
            }
          }
          //cout << "violations:" << violations << endl;
          //until no violation
          if (violations==0) break;
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

      //Copy active data
      for (int k=0;k<pme;k++){
        scr_mat(k,a) = act_me[k];
      }
      for (int k=0;k<pcme;k++){
        scr_mat(pme+k,a) = act_cme[k];
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
List cme_gaussian(NumericMatrix& XX_me, NumericMatrix& XX_cme, NumericVector& yy,
                  NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec,
                  NumericVector& gamma_vec, NumericVector& tau_vec,
                  NumericVector& XX_me_sl, NumericVector& XX_cme_sl, NumericVector& beta_vec, NumericVector& act_vec,
                  NumericVector& multiplier, NumericVector& multiplier_g,
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
  double mu = 0; (void)mu;

  // Extract the family type from CharacterVector

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

  //Check whether lambda is to be iterated or not
  bool lambda_it;
  int niter_1; //Number to iterate first
  int niter_2; //Number to iterate next
  if (gamma_vec.size()>1){ //Iterate on gamma and tau
    lambda_it = false;
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

  //Set all factors as active to begin
  vector<bool> act_me(pme,true); //Current active set
  vector<bool> act_cme(pcme,true);
  vector<bool> scr_me(pme,true); //Screened active set
  vector<bool> scr_cme(pcme,true);
  bool kkt_bool;


  //// Containers for linearized slopes Delta
  vector<double> m_me(pme,1);
  vector<double> m_cme(pcme,1);
  for (int i=0;i<pme;i++){
    m_me[i] = multiplier[i];
  }
  for (int i=0;i<pcme;i++){
    m_cme[i] = multiplier[pme+i];
  }

  vector<double> m_sib(pme,1);
  vector<double> m_cou(pme,1);
  for (int i=0;i<pme;i++){
    m_sib[i] = multiplier_g[2*i];
  }
  for (int i=0;i<pme;i++){
    m_cou[i] = multiplier_g[2*i+1];
  }


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
  //vector<double> eta(nn);
  vector<double> resid(nn); //Residual vector
  arma::cube resid_cube(nn,niter_1,niter_2);
  arma::mat resid_mat(nn,niter_1);
  arma::cube scr_cube(pme+pcme,niter_1,niter_2); //screening vector
  arma::mat scr_mat(pme+pcme,niter_1);


  //double inter= 0.0; //for intercept
  double inprod = 0.0; //inner product
  double cj = 0.0;
  double thresh = 0.0; (void)thresh; //threshold for screening
  int size = 0; (void)size;
  int num_act = 0; (void)num_act;
  int num_scr = 0; (void)num_scr;

  double ymean = 0.0;
  for (int i=0;i<nn;i++){
    ymean += (1.0/(double)nn)*yy[i];
  }
  for (int i=0; i<nn; i++) {
    resid[i] = yy(i) - ymean;
  }


  // Optimize for each penalty combination
  // #pragma omp parallel for
  for (int b=0; b<niter_2; b++){ //iterate over cousins...

    for (int a=0; a<niter_1; a++){ //iterate over siblings...

      //cout << "Tuning ... a = " << a << ", b = " << b << endl;//Update iteration variables
      if (lambda_it){
        lambda[0] = lambda_sib_vec[a];
        lambda[1] = lambda_cou_vec[b];
      }
      else{
       tau = tau_vec[a];
        gamma = gamma_vec[b];
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
        for (int i=0;i<nn;i++){//reset residuals
          resid[i] = yy(i) - ymean;
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
        //cout << "num_act: " << num_act << endl;
        if ( (lambda[0]+lambda[1]) >= lambda_max){
          goto cycend;
        }
      }

      // //Update screen/strong set
      num_scr = 0;
      for (int i=0;i<pme;i++){//reset screen flag
        scr_me[i] = true;
        num_scr ++;
      }
      for (int i=0;i<pcme;i++){
        scr_cme[i] = true;
        num_scr ++;
      }
      if (lambda_it && screen_ind && a>0 && b>0){

        for (int j=0;j<pme;j++){

          cj = crossprod(X_me, resid, nn, j)/(double(nn));
          scr_me[j] = screen_rule_me(cj, j, a, b, lambda_sib_vec, lambda_cou_vec, gamma, delta_sib, delta_cou,
                                     m_me, m_sib, m_cou);
          if(!scr_me[j]){
            num_scr --;
          }
        }

        for (int j=0;j<pme;j++){ //parent effect
          for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
            //for (int j=0; j<pcme; j++){

            int cmeind = j*(2*(pme-1))+k; //index for CME

            int condind = floor((double)k/2.0);
            if (condind >= j){
              condind ++;
            }

            cj = crossprod(X_cme, resid, nn, cmeind)/(double(nn));

            scr_cme[cmeind] = screen_rule_cme(cj, j, k, cmeind, condind, a, b, lambda_sib_vec, lambda_cou_vec, gamma, delta_sib, delta_cou,
                                              m_cme, m_sib, m_cou);
            if(!scr_cme[cmeind]){
              num_scr --;
            }
          }
        }
      }

      //cout << "strong set:" << num_scr << endl;


      // // RESET AFTER EACH RUN
      for (int i=0;i<pme;i++){//reset beta
        beta_me[i] = 0.0;
      }
      for (int i=0;i<pcme;i++){
        beta_cme[i] = 0.0;
      }
      for (int i=0;i<nn;i++){//reset residuals
        resid[i] = yy(i) - ymean;
      }

    //adaptive Recompute deltas
      for (int j=0; j<pme; j++){
        delta_sib[j] = lambda[0] * m_sib[j];
        delta_cou[j] = lambda[1] * m_cou[j];
      }
      for (int j=0; j<pme; j++){
        delta_sib[j] = delta_sib[j] * ( exp( -(tau/(lambda[0]*m_sib[j])) * m_me[j] * mcp(beta_me[j],lambda[0]*m_sib[j],gamma) ) );
        delta_cou[j] = delta_cou[j] * ( exp( -(tau/(lambda[1]*m_cou[j])) * m_me[j] * mcp(beta_me[j],lambda[1]*m_cou[j],gamma) ) );
      }
      for (int j=0;j<pme;j++){ //parent effect
        for (int k=0;k<(2*(pme-1));k++){ //conditioned effect
          int cmeind = j*(2*(pme-1))+k;
          int condind = floor((double)k/2.0);
          if (condind >= j){
            condind ++;
          }
          delta_sib[j] = delta_sib[j] * (exp(-(tau/(lambda[0]*m_sib[j])) * m_cme[cmeind] * mcp(beta_cme[cmeind],lambda[0]*m_sib[j],gamma) ));
          delta_cou[condind] = delta_cou[condind] * (exp(-(tau/(lambda[1]*m_cou[condind])) * m_cme[cmeind] * mcp(beta_cme[cmeind],lambda[1]*m_cou[condind],gamma) ));
        }
      }

      //Coordinate descent with warm active set resets
      for (int q=0; q<reset; q++){

        //Active set reset for it_warm iterations
        for (int m=0; m<it_warm; m++){
          chng_flag = coord_des_onerun_gaussian(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy,
                                                delta_sib, delta_cou, act_me, act_cme, beta_me, beta_cme, m_me, m_cme, m_sib, m_cou, resid);

        }

        //Update active set
        int num_act = 0; (void)num_act;
        for (int j=0;j<pme;j++){
          if ((abs(beta_me[j])>0.0||(act_vec[j]>0.0))){ //
            act_me[j] = true;
            num_act ++;
          }
          else{
            act_me[j] = false;
          }
        }
        for (int j=0;j<pcme;j++){
          if ((abs(beta_cme[j])>0.0)||(act_vec[j+pme]>0.0)){ //
            act_cme[j] = true;
            num_act ++;
          }
          else{
            act_cme[j] = false;
          }
        }

        //cout << "warm act:" << num_act << endl;

        //Cycle on active set
        it_inner = 0; //inner iteration count

        while (it_inner < it_max_reset){


          cont = true; //continue flag
          chng_flag = false; //change flag

          while (cont){

            //Increment count and update flags
            it_inner ++;
            chng_flag = coord_des_onerun_gaussian(pme, pcme, nn, lambda, cur_delta, chng_flag, tau, gamma, X_me, X_cme, yy,
                                                  delta_sib, delta_cou, act_me, act_cme, beta_me, beta_cme, m_me, m_cme, m_sib, m_cou, resid);

            //Update cont flag for termination
            if ( (it_inner >= it_max_reset)||(!chng_flag) ){
              cont = false;
            }
          }//end while

          //cout << "it_inner: " << it_inner << endl;

          //Rcout << accumulate(act.begin(),act.end(),0) << endl;

          //check violation in complementary active set\strong set, update active set
          int violations = 0;
          for (int j=0; j<pme; j++) {
            if (act_cme[j]==false && scr_me[j]==false) {
              inprod = crossprod(X_me, resid, nn, j)/(double(nn));
              kkt_bool = kkt(inprod, cur_delta);
              if (!kkt_bool) {
                act_me[j] = scr_me[j] = true;
                violations++;
              }
            }
          }
          for (int j=0; j<pcme; j++) {
            if (act_cme[j]==false && scr_cme[j]==false) {
              inprod = crossprod(X_cme, resid, nn, j)/(double(nn));
              kkt_bool = kkt(inprod, cur_delta);
              if (!kkt_bool) {
                act_cme[j] = scr_cme[j] = true;
                violations++;
              }
            }
          }
          //cout << "violations:" << violations << endl;
          //until no violation
          if (violations==0) break;
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
        scr_mat(k,a) = act_me[k];
      }
      for (int k=0;k<pcme;k++){
        scr_mat(pme+k,a) = act_cme[k];
      }

    }//end nlambda.sib (a)

    //Copy into beta_cube
    beta_cube.slice(b) = beta_mat;
    delta_sib_cube.slice(b) = delta_sib_mat;
    delta_cou_cube.slice(b) = delta_cou_mat;
    resid_cube.slice(b) = resid_mat;
    scr_cube.slice(b) = scr_mat;

  }//end nlambda.cou (b)


  for (int b=0; b<niter_2; b++){ //iterate over cousins...
    for (int a=0; a<niter_1; a++){ //iterate over siblings...
      //Rescale betas to original scale
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
                       Named("residuals") = resid_cube,
                       Named("nzero") = nz,
                       Named("lambda_sib") = lambda_sib_vec,
                       Named("lambda_cou") = lambda_cou_vec,
                       Named("act") = scr_cube,
                       Named("gamma") = gamma,
                       Named("tau") = tau
  ));

}
