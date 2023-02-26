// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mcp
double mcp(double beta, double lambda, double gamma);
RcppExport SEXP _glmcmenet_mcp(SEXP betaSEXP, SEXP lambdaSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(mcp(beta, lambda, gamma));
    return rcpp_result_gen;
END_RCPP
}
// cme
List cme(NumericMatrix& XX_me, NumericMatrix& XX_cme, NumericVector& yy, NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec, NumericVector& gamma_vec, NumericVector& tau_vec, NumericVector& XX_me_sl, NumericVector& XX_cme_sl, NumericVector& beta_vec, NumericVector& act_vec, double lambda_max, int it_max, int it_warm, int reset, bool screen_ind);
RcppExport SEXP _glmcmenet_cme(SEXP XX_meSEXP, SEXP XX_cmeSEXP, SEXP yySEXP, SEXP lambda_sib_vecSEXP, SEXP lambda_cou_vecSEXP, SEXP gamma_vecSEXP, SEXP tau_vecSEXP, SEXP XX_me_slSEXP, SEXP XX_cme_slSEXP, SEXP beta_vecSEXP, SEXP act_vecSEXP, SEXP lambda_maxSEXP, SEXP it_maxSEXP, SEXP it_warmSEXP, SEXP resetSEXP, SEXP screen_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX_me(XX_meSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX_cme(XX_cmeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_sib_vec(lambda_sib_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_cou_vec(lambda_cou_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_me_sl(XX_me_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_cme_sl(XX_cme_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta_vec(beta_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type act_vec(act_vecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_warm(it_warmSEXP);
    Rcpp::traits::input_parameter< int >::type reset(resetSEXP);
    Rcpp::traits::input_parameter< bool >::type screen_ind(screen_indSEXP);
    rcpp_result_gen = Rcpp::wrap(cme(XX_me, XX_cme, yy, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, lambda_max, it_max, it_warm, reset, screen_ind));
    return rcpp_result_gen;
END_RCPP
}
