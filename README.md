## The Curious Problem of the Normal Inverse Mean: Robustness and Shrinkage

### Soham Ghosh, Uttaran Chatterjee, and Jyotishka Datta

[arXiv:2410.20641](https://arxiv.org/abs/2410.20641)

This repository contains the `R` and `Stan` codes for Ghosh et al. (2025).

### Contents

- `Rcodes/experiments/` — R and Stan scripts to reproduce all simulation 
  experiments (Setup A and Setup B) in Section 6.
- `Rcodes/experiments/` — Code for the real-data analysis using the Gaia 
  Data Release 1 (GDR1) dataset in Section 7.
- MCMC convergence diagnostics (ACF plots, trace plots) and posterior 
  predictive checks (Appendix A.3--A.4).
- Sensitivity analyses with respect to measurement error and prior scale 
  (Appendix A.5--A.6).

### Dependencies

R packages: `rstan`, `tidyverse`, `ggplot2`  
Stan version: 2.x
