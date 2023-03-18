mcp <- function (beta, lambda, gamma)
{
  .Call(`_glmcmenet_mcp`, beta,lambda, gamma)
}
