## The Curious Problem of the Normal Inverse Mean: Robustness and Shrinkage

### Soham Ghosh, Uttaran Chatterjee, and Jyotishka Datta

[arXiv:2410.20641](https://arxiv.org/abs/2410.20641)

This repository contains the `R` and `Stan` codes for Ghosh et al. (2026).

### Contents

- - `Rcodes/experiments/simulation_codes.R` — R and Stan scripts reproducing:
  - The likelihood plot (Figure 2)
  - Setup A: squared error comparison across priors (Figure 4)
  - Setup B: RMSE and mean squared relative error under varying noise (Figure 5a--b)
  - Fixed-width rejection loss experiment (Figure 5c)
  - Near-zero posterior predictive discrepancy curves, 
    Jensen--Shannon divergence and $L^1$ discrepancy (Figure 9)
  - Sensitivity to prior scale $L \in \{500, 1000, 5000\}$ (Figure 11)
- `Rcodes/experiments/GDR1study.R` - R Code for the real-data analysis using the Gaia 
  Data Release 1 (GDR1) dataset in Section 7.
- `Rcodes/experiments/MCMC_diagnostics.r` - MCMC convergence diagnostics (ACF plots, trace plots) and posterior 
  predictive checks (Appendix A.3--A.4).

### Dependencies

R packages: `rstan`, `tidyverse`, `ggplot2`  
Stan (https://github.com/stan-dev). 
