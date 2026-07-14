## The Curious Problem of the Normal Inverse Mean: Robustness and Shrinkage

### Soham Ghosh, Uttaran Chatterjee, and Jyotishka Datta

[arXiv:2410.20641](https://arxiv.org/abs/2410.20641)

### Contents

- `Rcodes/experiments/simulation_codes.R` — Reproduces all simulation experiments 
  and figures in the paper, including:
  - Likelihood plot (Figure 2)
  - Setup A: squared error comparison across priors (Figure 4)
  - Setup B: RMSE and mean squared relative error under varying noise (Figure 5a--b)
  - Fixed-width rejection loss experiment (Figure 5c)
  - Near-zero posterior predictive discrepancy curves using Jensen--Shannon 
    divergence and L1 discrepancy (Figure 9)
  - Sensitivity to prior scale L ∈ {500, 1000, 5000} (Figure 11)

- `Rcodes/experiments/MCMC_diagnostics.R` — Runs all six Stan models with 4 chains, 
  5,000 iterations (2,000 warmup), and produces:
  - Trace plots for all priors (Figure 8)
  - ACF plots for all priors (Figure 7)
  - R-hat and ESS summary table (Table 2, Appendix A.3)

- `Rcodes/experiments/GDR1study.R` — Real-data analysis using the Gaia Data Release 1 
  (GDR1) dataset (Section 7). Reads `gdr1set01.csv`, computes posterior mode distance 
  estimates for 100 stars under six candidate priors.

### Dependencies

R packages: `rstan`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `bayesplot`, 
`posterior`, `purrr`, `MCMCpack`, `pracma`, `ggdensity`, `reshape2`, 
`latex2exp`, `gridExtra`, `PolynomF`, `mvtnorm`, `plyr`, `parallel`, `here`
