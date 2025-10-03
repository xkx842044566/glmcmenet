
# glmcmenet

<!-- badges: start -->
<!-- badges: end -->

`glmcmenet` fits generalized linear models with **Adaptive Bi-Level Variable
Selection of Conditional Main Effects (CMEs)**. It provides fast C++ (Rcpp /
Armadillo) implementations for:

- bi-level penalties that couple *sibling* and *cousin* CMEs,
- families: **gaussian**, **binomial**, **poisson**,
- cross-validation over `(gamma, tau)` and `(lambda.sib, lambda.cou)`,
- prediction for new ME/CME design matrices.

This repository contains the code for the method used in the paper:

> “Adaptive Bi-Level Variable Selection of Conditional Main Effects for
> Generalized Linear Models” by Kexin Xie and Xinwei Deng.


## Installation:
* the latest version: `devtools::install_github("xkx842044566/glmcmenet")`

## Reproducibility / Simulation Code

A companion repository with ready-to-run scripts that reproduce the paper’s simulation results is available here:

* [glmcmenet-simulation](https://github.com/xkx842044566/glmcmenet-simulations) *

## Report bugs：
* open an [issue](https://github.com/xkx842044566/glmcmenet/issues) or send an email to Kexin Xie at <kexinx@vt.edu>


