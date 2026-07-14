## The Curious Problem of the Normal Inverse Mean: Robustness and Shrinkage

### Soham Ghosh, Uttaran Chatterjee, and Jyotishka Datta

[arXiv:2410.20641](https://arxiv.org/abs/2410.20641)

This repository contains the `R` and `Stan` codes for Ghosh et al. (2026).

### Contents

- `Rcodes/experiments/simulation_codes.r` - R and Stan scripts to reproduce all simulation 
  experiments (Setup A and Setup B) in Section 6.
- `Rcodes/experiments/GDR1study.R` - R Code for the real-data analysis using the Gaia 
  Data Release 1 (GDR1) dataset in Section 7.
- `Rcodes/experiments/MCMC_diagnostics.r` - MCMC convergence diagnostics (ACF plots, trace plots) and posterior 
  predictive checks (Appendix A.3--A.4).

### Dependencies

R packages: `rstan`, `tidyverse`, `ggplot2`  
Stan version: 2.x
