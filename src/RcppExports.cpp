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
// cme_str
List cme_str(NumericMatrix& XX_me, NumericMatrix& XX_cme, NumericVector& yy, CharacterVector& family, NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec, NumericVector& gamma_vec, NumericVector& tau_vec, NumericVector& XX_me_sl, NumericVector& XX_cme_sl, NumericVector& beta_vec, NumericVector& act_vec, NumericVector& multiplier, double lambda_max, int it_max, int it_warm, int reset, bool screen_ind);
RcppExport SEXP _glmcmenet_cme_str(SEXP XX_meSEXP, SEXP XX_cmeSEXP, SEXP yySEXP, SEXP familySEXP, SEXP lambda_sib_vecSEXP, SEXP lambda_cou_vecSEXP, SEXP gamma_vecSEXP, SEXP tau_vecSEXP, SEXP XX_me_slSEXP, SEXP XX_cme_slSEXP, SEXP beta_vecSEXP, SEXP act_vecSEXP, SEXP multiplierSEXP, SEXP lambda_maxSEXP, SEXP it_maxSEXP, SEXP it_warmSEXP, SEXP resetSEXP, SEXP screen_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX_me(XX_meSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX_cme(XX_cmeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_sib_vec(lambda_sib_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_cou_vec(lambda_cou_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_me_sl(XX_me_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_cme_sl(XX_cme_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta_vec(beta_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type act_vec(act_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type multiplier(multiplierSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_warm(it_warmSEXP);
    Rcpp::traits::input_parameter< int >::type reset(resetSEXP);
    Rcpp::traits::input_parameter< bool >::type screen_ind(screen_indSEXP);
    rcpp_result_gen = Rcpp::wrap(cme_str(XX_me, XX_cme, yy, family, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, lambda_max, it_max, it_warm, reset, screen_ind));
    return rcpp_result_gen;
END_RCPP
}
// cme
List cme(NumericMatrix& XX_me, NumericMatrix& XX_cme, NumericVector& yy, CharacterVector& family, NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec, NumericVector& gamma_vec, NumericVector& tau_vec, NumericVector& XX_me_sl, NumericVector& XX_cme_sl, NumericVector& beta_vec, NumericVector& act_vec, NumericVector& multiplier, NumericVector& multiplier_g, double lambda_max, int it_max, int it_warm, int reset, bool screen_ind);
RcppExport SEXP _glmcmenet_cme(SEXP XX_meSEXP, SEXP XX_cmeSEXP, SEXP yySEXP, SEXP familySEXP, SEXP lambda_sib_vecSEXP, SEXP lambda_cou_vecSEXP, SEXP gamma_vecSEXP, SEXP tau_vecSEXP, SEXP XX_me_slSEXP, SEXP XX_cme_slSEXP, SEXP beta_vecSEXP, SEXP act_vecSEXP, SEXP multiplierSEXP, SEXP multiplier_gSEXP, SEXP lambda_maxSEXP, SEXP it_maxSEXP, SEXP it_warmSEXP, SEXP resetSEXP, SEXP screen_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX_me(XX_meSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX_cme(XX_cmeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_sib_vec(lambda_sib_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_cou_vec(lambda_cou_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_me_sl(XX_me_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_cme_sl(XX_cme_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta_vec(beta_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type act_vec(act_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type multiplier(multiplierSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type multiplier_g(multiplier_gSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_warm(it_warmSEXP);
    Rcpp::traits::input_parameter< int >::type reset(resetSEXP);
    Rcpp::traits::input_parameter< bool >::type screen_ind(screen_indSEXP);
    rcpp_result_gen = Rcpp::wrap(cme(XX_me, XX_cme, yy, family, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, multiplier_g, lambda_max, it_max, it_warm, reset, screen_ind));
    return rcpp_result_gen;
END_RCPP
}
// cme_wls
List cme_wls(NumericMatrix& XX, NumericVector& yy, CharacterVector& family, NumericVector& K1, NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec, NumericVector& gamma_vec, NumericVector& tau_vec, NumericVector& XX_sl, NumericVector& beta_vec, NumericVector& act_vec, NumericVector& multiplier, double lambda_max, int it_max, int it_warm, int reset, bool screen_ind);
RcppExport SEXP _glmcmenet_cme_wls(SEXP XXSEXP, SEXP yySEXP, SEXP familySEXP, SEXP K1SEXP, SEXP lambda_sib_vecSEXP, SEXP lambda_cou_vecSEXP, SEXP gamma_vecSEXP, SEXP tau_vecSEXP, SEXP XX_slSEXP, SEXP beta_vecSEXP, SEXP act_vecSEXP, SEXP multiplierSEXP, SEXP lambda_maxSEXP, SEXP it_maxSEXP, SEXP it_warmSEXP, SEXP resetSEXP, SEXP screen_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_sib_vec(lambda_sib_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_cou_vec(lambda_cou_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_sl(XX_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta_vec(beta_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type act_vec(act_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type multiplier(multiplierSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_warm(it_warmSEXP);
    Rcpp::traits::input_parameter< int >::type reset(resetSEXP);
    Rcpp::traits::input_parameter< bool >::type screen_ind(screen_indSEXP);
    rcpp_result_gen = Rcpp::wrap(cme_wls(XX, yy, family, K1, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_sl, beta_vec, act_vec, multiplier, lambda_max, it_max, it_warm, reset, screen_ind));
    return rcpp_result_gen;
END_RCPP
}
// cme_gaussian
List cme_gaussian(NumericMatrix& XX, NumericVector& yy, NumericVector& K1, NumericVector& lambda_sib_vec, NumericVector& lambda_cou_vec, NumericVector& gamma_vec, NumericVector& tau_vec, NumericVector& XX_sl, NumericVector& beta_vec, NumericVector& act_vec, NumericVector& multiplier, double lambda_max, int it_max, int it_warm, int reset, bool screen_ind);
RcppExport SEXP _glmcmenet_cme_gaussian(SEXP XXSEXP, SEXP yySEXP, SEXP K1SEXP, SEXP lambda_sib_vecSEXP, SEXP lambda_cou_vecSEXP, SEXP gamma_vecSEXP, SEXP tau_vecSEXP, SEXP XX_slSEXP, SEXP beta_vecSEXP, SEXP act_vecSEXP, SEXP multiplierSEXP, SEXP lambda_maxSEXP, SEXP it_maxSEXP, SEXP it_warmSEXP, SEXP resetSEXP, SEXP screen_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yy(yySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_sib_vec(lambda_sib_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda_cou_vec(lambda_cou_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type tau_vec(tau_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XX_sl(XX_slSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta_vec(beta_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type act_vec(act_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type multiplier(multiplierSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_warm(it_warmSEXP);
    Rcpp::traits::input_parameter< int >::type reset(resetSEXP);
    Rcpp::traits::input_parameter< bool >::type screen_ind(screen_indSEXP);
    rcpp_result_gen = Rcpp::wrap(cme_gaussian(XX, yy, K1, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_sl, beta_vec, act_vec, multiplier, lambda_max, it_max, it_warm, reset, screen_ind));
    return rcpp_result_gen;
END_RCPP
}
