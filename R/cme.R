cme <- function (XX_me, XX_cme, yy, family, lambda_sib_vec, lambda_cou_vec,
                 gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier, multiplier_g,
                 lambda_max, it_max, it_warm, reset, screen_ind)
{
  .Call(`_glmcmenet_cme`, XX_me, XX_cme, yy, family, lambda_sib_vec,
        lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl,
        beta_vec, act_vec, multiplier,  multiplier_g, lambda_max, it_max, it_warm, reset,
        screen_ind)
}


cme_str <- function (XX_me, XX_cme, yy, family, lambda_sib_vec, lambda_cou_vec,
                 gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec, multiplier,
                 lambda_max, it_max, it_warm, reset, screen_ind)
{
  .Call(`_glmcmenet_cme_str`, XX_me, XX_cme, yy, family, lambda_sib_vec,
        lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl,
        beta_vec, act_vec, multiplier, lambda_max, it_max, it_warm, reset,
        screen_ind)
}


cme_wls <- function (XX, yy, family, K1, lambda_sib_vec, lambda_cou_vec,
                     gamma_vec, tau_vec, XX_sl, beta_vec, act_vec, multiplier,
                     lambda_max, it_max, it_warm, reset, screen_ind)
{
  .Call(`_glmcmenet_cme_wls`, XX, yy, family, K1, lambda_sib_vec,
        lambda_cou_vec, gamma_vec, tau_vec, XX_sl,
        beta_vec, act_vec, multiplier, lambda_max, it_max, it_warm, reset,
        screen_ind)
}

cme_gaussian <- function (XX, yy, K1, lambda_sib_vec, lambda_cou_vec,
                          gamma_vec, tau_vec, XX_sl, beta_vec, act_vec, multiplier,
                          lambda_max, it_max, it_warm, reset, screen_ind)
{
  .Call(`_glmcmenet_cme_gaussian`, XX, yy, K1, lambda_sib_vec,
        lambda_cou_vec, gamma_vec, tau_vec, XX_sl,
        beta_vec, act_vec, multiplier, lambda_max, it_max, it_warm, reset,
        screen_ind)
}
