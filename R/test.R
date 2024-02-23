if (requireNamespace("tinytest", quietly=TRUE)) {
  if (length(unclass(packageVersion("glmcmenet"))[[1]]) == 4 | Sys.getenv('R_FORCE_TEST') == 'TRUE') {
    tinytest::test_package("glmcmenet", pattern="^[^_].*\\.[rR]$")
  }
}
