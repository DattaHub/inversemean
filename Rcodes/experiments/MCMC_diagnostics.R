############################################################
## 4 chains
## 5,000 iterations per chain
## 2,000 warmup iterations per chain
## 3,000 retained draws per chain
############################################################

library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
rstan_options(auto_write = TRUE)

n_chains <- 4L
n_iter   <- 5000L
n_warmup <- 2000L
n_cores  <- min(n_chains, parallel::detectCores())

options(mc.cores = n_cores)

set.seed(123)

J <- 200L

omega_obs <- seq(
  from = 0.01,
  to = 8,
  length.out = J
)

sigma_omega <- 0.045

stan_data_common <- list(
  J = J,
  Y = omega_obs,
  sigma_omega = sigma_omega
)


# Gamma(3, 10)
gamma_shape <- 3
gamma_rate  <- 10

# Inverse-Gamma(4, 1)
ig_shape <- 4
ig_scale <- 1

# RG+(0, 10)
rg_sd <- sqrt(10)

# Weibull(shape = 0.5, scale = 1).
weibull_shape <- 0.5
weibull_scale <- 1

# Half-Cauchy scale.
hc_scale <- 1

# Product Half-Cauchy multiplicative scale.
phc_scale <- 1

############################
## Half-Cauchy
############################

hc_code <- "
data {
  int<lower=1> J;
  vector[J] Y;
  real<lower=0> sigma_omega;
  real<lower=0> hc_scale;
}

parameters {
  vector<lower=0>[J] u;
}

model {
  u ~ cauchy(0, hc_scale);

  for (j in 1:J) {
    Y[j] ~ normal(1.0 / u[j], sigma_omega);
  }
}
"

############################
## Reciprocal-Gaussian
############################

rg_code <- "
functions {
  real reciprocal_gaussian_lpdf(real x, real sigma) {
    return
      0.5 * log(2.0 / pi())
      - log(sigma)
      - 2.0 * log(x)
      - 1.0 / (2.0 * square(sigma) * square(x));
  }
}

data {
  int<lower=1> J;
  vector[J] Y;
  real<lower=0> sigma_omega;
  real<lower=0> rg_sd;
}

parameters {
  vector<lower=0>[J] u;
}

model {
  for (j in 1:J) {
    u[j] ~ reciprocal_gaussian(rg_sd);
    Y[j] ~ normal(1.0 / u[j], sigma_omega);
  }
}
"

############################
## Weibull
############################

weibull_code <- "
data {
  int<lower=1> J;
  vector[J] Y;
  real<lower=0> sigma_omega;
  real<lower=0> weibull_shape;
  real<lower=0> weibull_scale;
}

parameters {
  vector<lower=0>[J] u;
}

model {
  u ~ weibull(weibull_shape, weibull_scale);

  for (j in 1:J) {
    Y[j] ~ normal(1.0 / u[j], sigma_omega);
  }
}
"

############################
## Gamma
############################

gamma_code <- "
data {
  int<lower=1> J;
  vector[J] Y;
  real<lower=0> sigma_omega;
  real<lower=0> gamma_shape;
  real<lower=0> gamma_rate;
}

parameters {
  vector<lower=0>[J] u;
}

model {
  u ~ gamma(gamma_shape, gamma_rate);

  for (j in 1:J) {
    Y[j] ~ normal(1.0 / u[j], sigma_omega);
  }
}
"

############################
## Inverse-Gamma
############################

ig_code <- "
data {
  int<lower=1> J;
  vector[J] Y;
  real<lower=0> sigma_omega;
  real<lower=0> ig_shape;
  real<lower=0> ig_scale;
}

parameters {
  vector<lower=0>[J] u;
}

model {
  u ~ inv_gamma(ig_shape, ig_scale);

  for (j in 1:J) {
    Y[j] ~ normal(1.0 / u[j], sigma_omega);
  }
}
"

############################
## Product Half-Cauchy
############################

phc_code <- "
data {
  int<lower=1> J;
  vector[J] Y;
  real<lower=0> sigma_omega;
  real<lower=0> phc_scale;
}

parameters {
  vector<lower=0, upper=pi()/2>[J] theta_1;
  vector<lower=0, upper=pi()/2>[J] theta_2;
}

transformed parameters {
  vector<lower=0>[J] u;

  for (j in 1:J) {
    real log_tau;
    real log_lambda;

    log_tau =
      log(sin(theta_1[j])) -
      log(cos(theta_1[j]));

    log_lambda =
      log(sin(theta_2[j])) -
      log(cos(theta_2[j]));

    u[j] =
      phc_scale *
      exp(log_tau + log_lambda);
  }
}

model {
  theta_1 ~ uniform(0, pi()/2);
  theta_2 ~ uniform(0, pi()/2);

  for (j in 1:J) {
    Y[j] ~ normal(1.0 / u[j], sigma_omega);
  }
}
"

############################################################

hc_mod <- stan_model(
  model_code = hc_code,
  model_name = "half_cauchy_distance"
)

rg_mod <- stan_model(
  model_code = rg_code,
  model_name = "reciprocal_gaussian_distance"
)

weibull_mod <- stan_model(
  model_code = weibull_code,
  model_name = "weibull_distance"
)

gamma_mod <- stan_model(
  model_code = gamma_code,
  model_name = "gamma_distance"
)

ig_mod <- stan_model(
  model_code = ig_code,
  model_name = "inverse_gamma_distance"
)

phc_mod <- stan_model(
  model_code = phc_code,
  model_name = "product_half_cauchy_distance"
)


run_stan_fit <- function(
    model,
    data,
    seed,
    adapt_delta = 0.95,
    max_treedepth = 12L
) {
  
  sampling(
    object = model,
    data = data,
    chains = n_chains,
    iter = n_iter,
    warmup = n_warmup,
    cores = n_cores,
    seed = seed,
    refresh = 100,
    init = "random",
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
}


hc_fit <- run_stan_fit(
  model = hc_mod,
  data = c(
    stan_data_common,
    list(hc_scale = hc_scale)
  ),
  seed = 101,
  adapt_delta = 0.97,
  max_treedepth = 13L
)


rg_fit <- run_stan_fit(
  model = rg_mod,
  data = c(
    stan_data_common,
    list(rg_sd = rg_sd)
  ),
  seed = 102,
  adapt_delta = 0.97,
  max_treedepth = 13L
)


weibull_fit <- run_stan_fit(
  model = weibull_mod,
  data = c(
    stan_data_common,
    list(
      weibull_shape = weibull_shape,
      weibull_scale = weibull_scale
    )
  ),
  seed = 103,
  adapt_delta = 0.95,
  max_treedepth = 12L
)

gamma_fit <- run_stan_fit(
  model = gamma_mod,
  data = c(
    stan_data_common,
    list(
      gamma_shape = gamma_shape,
      gamma_rate = gamma_rate
    )
  ),
  seed = 104,
  adapt_delta = 0.95,
  max_treedepth = 12L
)


ig_fit <- run_stan_fit(
  model = ig_mod,
  data = c(
    stan_data_common,
    list(
      ig_shape = ig_shape,
      ig_scale = ig_scale
    )
  ),
  seed = 105,
  adapt_delta = 0.97,
  max_treedepth = 13L
)


phc_fit <- run_stan_fit(
  model = phc_mod,
  data = c(
    stan_data_common,
    list(phc_scale = phc_scale)
  ),
  seed = 106,
  adapt_delta = 0.995,
  max_treedepth = 15L
)


fits <- list(
  "Half-Cauchy" = hc_fit,
  "Reciprocal-Gaussian" = rg_fit,
  "Gamma" = gamma_fit,
  "Inverse-Gamma" = ig_fit,
  "Product Half-Cauchy" = phc_fit,
  "Weibull" = weibull_fit
)

treedepth_limits <- c(
  "Half-Cauchy" = 13L,
  "Reciprocal-Gaussian" = 13L,
  "Gamma" = 12L,
  "Inverse-Gamma" = 13L,
  "Product Half-Cauchy" = 15L,
  "Weibull" = 12L
)

saveRDS(
  fits,
  file = "stan_fits_4chains.rds"
)


selected_index <- 10L
selected_parameter <- paste0(
  "u[",
  selected_index,
  "]"
)

selected_true_parallax <- omega_obs[selected_index]
selected_true_distance <- 1 / selected_true_parallax

selected_fractional_error <-
  sigma_omega * selected_true_distance


extract_chain_draws <- function(
    fit,
    parameter
) {
  
  draw_array <- rstan::extract(
    object = fit,
    pars = parameter,
    permuted = FALSE,
    inc_warmup = FALSE
  )
  
  draw_matrix <- draw_array[, , 1, drop = TRUE]
  
  n_saved <- nrow(draw_matrix)
  n_chain <- ncol(draw_matrix)
  
  data.frame(
    Iteration = rep(
      seq_len(n_saved),
      times = n_chain
    ),
    Chain = factor(
      rep(
        seq_len(n_chain),
        each = n_saved
      ),
      labels = paste(
        "Chain",
        seq_len(n_chain)
      )
    ),
    Value = as.vector(draw_matrix)
  )
}


chain_colors <- c(
  "Chain 1" = "#0072B2",
  "Chain 2" = "#D55E00",
  "Chain 3" = "#009E73",
  "Chain 4" = "#CC79A7"
)


make_trace_panel <- function(
    fit,
    prior_name,
    parameter
) {
  
  trace_df <- extract_chain_draws(
    fit = fit,
    parameter = parameter
  )
  
  ggplot(
    trace_df,
    aes(
      x = Iteration,
      y = Value,
      color = Chain,
      group = Chain
    )
  ) +
    geom_line(
      size = 0.28,
      alpha = 0.72
    ) +
    scale_color_manual(
      values = chain_colors
    ) +
    labs(
      title = prior_name,
      x = "Post-warmup iteration",
      y = "Posterior distance",
      color = NULL
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "bottom",
      legend.text = element_text(
        size = 9
      ),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = 1,
          alpha = 1
        ),
        nrow = 1
      )
    )
}

trace_panels <- lapply(
  names(fits),
  function(prior_name) {
    
    make_trace_panel(
      fit = fits[[prior_name]],
      prior_name = prior_name,
      parameter = selected_parameter
    )
  }
)

trace_plot_all <- wrap_plots(
  trace_panels,
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(
    title = ""
  ) &
  theme(
    legend.position = "bottom"
  )

print(trace_plot_all)

# ggsave(
#   filename =
#     "trace_allpriors_4chains_superimposed.pdf",
#   plot = trace_plot_all,
#   width = 14,
#   height = 8.5,
#   units = "in"
# )

# ggsave(
#   filename =
#     "trace_allpriors_4chains_superimposed.png",
#   plot = trace_plot_all,
#   width = 14,
#   height = 8.5,
#   units = "in",
#   dpi = 300
# )


make_acf_data <- function(
    fit,
    parameter,
    lag_max = 50L
) {
  
  draw_df <- extract_chain_draws(
    fit = fit,
    parameter = parameter
  )
  
  acf_df <- draw_df %>%
    group_by(Chain) %>%
    group_modify(
      ~ {
        acf_object <- stats::acf(
          .x$Value,
          lag.max = lag_max,
          plot = FALSE,
          demean = TRUE
        )
        
        data.frame(
          Lag = as.numeric(
            acf_object$lag
          ),
          ACF = as.numeric(
            acf_object$acf
          )
        )
      }
    ) %>%
    ungroup()
  
  acf_df %>%
    filter(Lag > 0)
}


make_acf_panel <- function(
    fit,
    prior_name,
    parameter,
    lag_max = 50L
) {
  
  acf_df <- make_acf_data(
    fit = fit,
    parameter = parameter,
    lag_max = lag_max
  )
  
  ggplot(
    acf_df,
    aes(
      x = Lag,
      y = ACF,
      color = Chain,
      group = Chain
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      size = 0.35
    ) +
    geom_line(
      size = 0.65,
      alpha = 0.85
    ) +
    geom_point(
      size = 0.8,
      alpha = 0.8
    ) +
    scale_color_manual(
      values = chain_colors
    ) +
    labs(
      title = prior_name,
      x = "Lag",
      y = "Autocorrelation",
      color = NULL
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "bottom",
      legend.text = element_text(
        size = 9
      ),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = 1,
          alpha = 1
        ),
        nrow = 1
      )
    )
}

acf_panels <- lapply(
  names(fits),
  function(prior_name) {
    
    make_acf_panel(
      fit = fits[[prior_name]],
      prior_name = prior_name,
      parameter = selected_parameter,
      lag_max = 50L
    )
  }
)

acf_plot_all <- wrap_plots(
  acf_panels,
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(
    title = ""
  ) &
  theme(
    legend.position = "bottom"
  )

print(acf_plot_all)

# ggsave(
#   filename =
#     "acf_allpriors_4chains_superimposed.pdf",
#   plot = acf_plot_all,
#   width = 14,
#   height = 8.5,
#   units = "in"
# )

# ggsave(
#   filename =
#     "acf_allpriors_4chains_superimposed.png",
#   plot = acf_plot_all,
#   width = 14,
#   height = 8.5,
#   units = "in",
#   dpi = 300
# )

summarize_rhat_ess <- function(
    fit,
    prior_name
) {
  
  fit_summary <- rstan::summary(
    fit,
    pars = "u"
  )$summary
  
  data.frame(
    Prior = prior_name,
    
    Mean_Rhat = mean(
      fit_summary[, "Rhat"],
      na.rm = TRUE
    ),
    
    Max_Rhat = max(
      fit_summary[, "Rhat"],
      na.rm = TRUE
    ),
    
    Mean_ESS = mean(
      fit_summary[, "n_eff"],
      na.rm = TRUE
    ),
    
    Median_ESS = median(
      fit_summary[, "n_eff"],
      na.rm = TRUE
    ),
    
    Min_ESS = min(
      fit_summary[, "n_eff"],
      na.rm = TRUE
    )
  )
}

rhat_ess_table <- bind_rows(
  lapply(
    names(fits),
    function(prior_name) {
      
      summarize_rhat_ess(
        fit = fits[[prior_name]],
        prior_name = prior_name
      )
    }
  )
)

print(rhat_ess_table)

