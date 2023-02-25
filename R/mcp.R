mcp <- function (beta, v, lambda, gamma)
{
  .Call(`_glmcmenet_mcp`, beta, v, lambda, gamma)
}
