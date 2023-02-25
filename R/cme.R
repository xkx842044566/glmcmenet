cme <- function (XX_me, XX_cme, yy, lambda_sib_vec, lambda_cou_vec,
          gamma_vec, tau_vec, XX_me_sl, XX_cme_sl, beta_vec, act_vec,
          lambda_max, it_max, it_warm, reset, screen_ind)
{
  .Call(`_glmcmenet_cme`, XX_me, XX_cme, yy, lambda_sib_vec,
        lambda_cou_vec, gamma_vec, tau_vec, XX_me_sl, XX_cme_sl,
        beta_vec, act_vec, lambda_max, it_max, it_warm, reset,
        screen_ind)
}
