# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

mcp <- function(beta, lambda, gamma) {
    .Call(`_glmcmenet_mcp`, beta, lambda, gamma)
}

cme <- function(XX_me, XX_cme, yy, family, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, multiplier_g, lambda_max, it_max, it_warm, reset, screen_ind) {
    .Call(`_glmcmenet_cme`, XX_me, XX_cme, yy, family, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, multiplier_g, lambda_max, it_max, it_warm, reset, screen_ind)
}

cme_gaussian <- function(XX_me, XX_cme, yy, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, multiplier_g, lambda_max, it_max, it_warm, reset, screen_ind) {
    .Call(`_glmcmenet_cme_gaussian`, XX_me, XX_cme, yy, lambda_sib_vec, lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, multiplier_g, lambda_max, it_max, it_warm, reset, screen_ind)
}

