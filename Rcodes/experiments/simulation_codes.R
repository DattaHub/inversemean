################ Code submission for "The Curious Problem of the Normal Inverse Mean: Robustness and Shrinkage" ##########################
#############################################################################################################################
#############################################################################################################################
library(plyr)
library(dplyr)
library(ggplot2)
library(rstan)
library(parallel)
library(reshape2)
library(here)
library(posterior)
library(purrr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
library(gridExtra)
library(PolynomF)
library(mvtnorm) 
library(latex2exp)
library(patchwork) 
library(MCMCpack)   
library(pracma)     
################# Code for reproducing the likelihood plot ######################
f <- function(r, omega) {
  exp(-0.5 * (1/r - omega)^2)
}
omegas <- c(0.01, 0.1, 0.5, 1, 5, 10, 1000)

# --- Panel (a): Small r values ---
r_small <- seq(1e-3, 5, length.out = 1000)
df_small <- expand.grid(r = r_small, omega = omegas) %>%
  mutate(y = map2_dbl(r, omega, f))

plot_a <-
  ggplot(df_small, aes(x = r, y = y, colour = as.factor(omega))) +
  geom_line(linewidth = 1) +
  labs(
    title = "",
    x = "r",
    y = TeX("$L(r ; \\omega)$"),
    colour = TeX("$\\omega$")
  ) +
  theme_minimal(base_size = 16)

# --- Panel (b): Large r values ---
r_large <- seq(1e-3, 100, length.out = 2000) 
df_large <- expand.grid(r = r_large, omega = omegas) %>%
  mutate(y = map2_dbl(r, omega, f))

plot_b <-
  ggplot(df_large, aes(x = r, y = y, colour = as.factor(omega))) +
  geom_line(linewidth = 1) +
  labs(
    title = "",
    x = "r",
    y = NULL, 
    colour = TeX("$\\omega$")
  ) +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_blank()) 

# combined_plot <- plot_a + plot_b
# ggsave("Enlarged-Likelihood-plot_two_panel.png", plot = combined_plot, width = 12, height = 5, dpi = 400)
print(combined_plot)

####### Codes for reproducing the simulation results (setup A) #########
########################################################################

#################### Compare Squared Error #######################
##################################################################
J = 200
y= seq(0, 8, length.out = J)
c.data = list('J'=J, 'Y' = y)
stan.iters = 5000

####### (1) Half-Cauchy prior ##########
cauchy_code <- "data {
int<lower=0> J;
vector[J] Y;
}
parameters {
real<lower=0> u[J];
}
model {
for (i in 1:J){
u[i] ~ cauchy(0, 1);
Y[i] ~ normal(1/u[i], 0.045);
}
}"

cauchy.mod <- stan_model(model_code = cauchy_code)
cauchy.fit <- sampling(cauchy.mod, data=c.data, iter=5000, chains=1, verbose=TRUE, warmup=2500, init=0, seed=123, refresh=10)
cauchy.array <- rstan::extract(cauchy.fit, permuted = TRUE)
cauchy.u.df <- as.data.frame(cauchy.array[["u"]])
cauchy.u.summary <- summarise_draws(cauchy.u.df) 


indivmean.dat.Cauchy <- data.frame(value=cauchy.u.summary$median, upper = cauchy.u.summary$q95, lower = cauchy.u.summary$q5,method = "Cauchy (0,1)")

indivmean.dat.Cauchy <- cbind(indivmean.dat.Cauchy, "Y" = y)
Cauchy.data <- indivmean.dat.Cauchy[-1,]
Cauchy.data$rtrue <- 1 / Cauchy.data$Y

se_cauchy <- (Cauchy.data$value - Cauchy.data$rtrue)^2
bias_cauchy <- Cauchy.data$value - Cauchy.data$rtrue

traceplot(cauchy.fit, pars = "u")
print(summary(cauchy.fit)$summary[,"Rhat"])

 
posterior_samples <- extract(cauchy.fit)

J <- c.data$J
Y_rep <- matrix(NA, nrow = nrow(posterior_samples$u), ncol = J)

for (i in 1:J) {
  u_samples <- posterior_samples$u[, i]
  Y_rep[, i] <- rnorm(nrow(Y_rep), mean = 1/u_samples, sd = 0.045)
}

#ppc_dens_overlay(c.data$Y, Y_rep)

ppc_ecdf_overlay(c.data$Y, Y_rep)

ppc_dens_overlay(c.data$Y, Y_rep[1:25, ])

######################################################
###### (2) Reciprocal Gaussian #########
RG_code <- "functions {
real RG_log (real x, real lambda){
real prob;
real lprob;
prob <- (sqrt(2/pi())*1/lambda*(1/x^2)*exp(-1/(2*lambda^2*x^2)));
lprob <- log(prob); 
return lprob;
}
}

data {
int<lower=0> J;
vector[J] Y;
}
parameters {
real<lower=0> u[J];
}
model {
for (i in 1:J){
u[i] ~ RG(1);
Y[i] ~ normal(1/u[i], 0.045);
}
}"
RG.mod <- stan_model(model_code = RG_code)
RG.fit <- sampling(RG.mod, data=c.data, iter=5000, chains=1, verbose=TRUE, warmup=2500, init=0, seed=123, refresh=10)
RG.array <- rstan::extract(RG.fit, permuted = TRUE)
RG.u.df <- as.data.frame(RG.array[["u"]])
RG.u.modes <- apply(RG.u.df,2,max)
RG.u.summary <- summarise_draws(RG.u.df) 
RG.u.summary <- cbind(RG.u.summary, RG.u.modes)
indivmean.dat.RG <- data.frame(value=RG.u.summary$median, upper = RG.u.summary$q95, lower = RG.u.summary$q5,method = "RG (0,1)")
indivmean.dat.RG <- cbind(indivmean.dat.RG, "Y" = y)


RG.data <- indivmean.dat.RG[-1,]
RG.data$rtrue <- 1 / RG.data$Y
se_RG <- (RG.data$value - RG.data$rtrue)^2
bias_RG <- RG.data$value - RG.data$rtrue

posterior_samples <- extract(RG.fit)

J <- c.data$J
Y_rep <- matrix(NA, nrow = nrow(posterior_samples$u), ncol = J)

for (i in 1:J) {
  u_samples <- posterior_samples$u[, i]
  Y_rep[, i] <- rnorm(nrow(Y_rep), mean = 1/u_samples, sd = 0.045)
}

ppc_dens_overlay(c.data$Y, Y_rep)

ppc_ecdf_overlay(c.data$Y, Y_rep)
######################################################
######## (3) Weibull #########
weibull_code <- "data {
int<lower=0> J;
vector[J] Y;
}
parameters {
real<lower=0> u[J];
}
model {
for (i in 1:J){
u[i] ~ weibull(0.5, 1);
Y[i] ~ normal(1/u[i], 0.045);
}
}"

w.mod <- stan_model(model_code = weibull_code)
w.fit <- sampling(w.mod, data=c.data, iter=5000, chains=1, verbose=TRUE, warmup=2500, init=0, seed=123, refresh=10)
w.array <- rstan::extract(w.fit, permuted = TRUE)
w.u.df <- as.data.frame(w.array[["u"]])
w.u.summary <- summarise_draws(w.u.df) 
indivmean.dat.Weibull <- data.frame(value=w.u.summary$mean, upper = w.u.summary$q95, lower = w.u.summary$q5,method = "Weibull(\u03b1=0.5, \u03C3=1)")
indivmean.dat.Weibull <- cbind(indivmean.dat.Weibull, "Y" = y)

Weibull.data <- indivmean.dat.Weibull[-1,]
Weibull.data$rtrue <- 1 / Weibull.data$Y
se_Weibull <- (Weibull.data$value - Weibull.data$rtrue)^2
bias_Weibull <- Weibull.data$value - Weibull.data$rtrue

posterior_samples <- extract(w.fit)

J <- c.data$J
Y_rep <- matrix(NA, nrow = nrow(posterior_samples$u), ncol = J)

for (i in 1:J) {
  u_samples <- posterior_samples$u[, i]
  Y_rep[, i] <- rnorm(nrow(Y_rep), mean = 1/u_samples, sd = 0.045)
}

ppc_dens_overlay(c.data$Y, Y_rep)

ppc_ecdf_overlay(c.data$Y, Y_rep)

ppc_pit_ecdf(c.data$Y, Y_rep, prob = 0.99, plot_diff = FALSE)

######################################################
####### (4) Gamma/ ED #######
gamma_code <- "data {
int<lower=0> J;
vector[J] Y;
}
parameters {
real<lower=0> u[J];
}
model {
for (i in 1:J){
u[i] ~ gamma(3, 10);
Y[i] ~ normal(1/u[i], 0.045);
}
}"
gamma.mod <- stan_model(model_code = gamma_code)
gamma.fit <- sampling(gamma.mod, data=c.data, iter=5000, chains=1, verbose=TRUE, warmup=2500, init=0, seed=123, refresh=10)
gamma.array <- rstan::extract(gamma.fit, permuted = TRUE)
gamma.u.df <- as.data.frame(gamma.array[["u"]])
gamma.u.summary <- summarise_draws(gamma.u.df) 
indivmean.dat.Gamma <- data.frame(value=gamma.u.summary$mean, upper = gamma.u.summary$q95, lower = gamma.u.summary$q5,method = "Gamma (3,10)")
indivmean.dat.Gamma <- cbind(indivmean.dat.Gamma, "Y" = y)

Gamma.data <- indivmean.dat.Gamma[-1,]
Gamma.data$rtrue <- 1 / Gamma.data$Y
se_Gamma <- (Gamma.data$value - Gamma.data$rtrue)^2
bias_Gamma <- Gamma.data$value - Gamma.data$rtrue

posterior_samples <- extract(gamma.fit)

J <- c.data$J
Y_rep <- matrix(NA, nrow = nrow(posterior_samples$u), ncol = J)

for (i in 1:J) {
  u_samples <- posterior_samples$u[, i]
  Y_rep[, i] <- rnorm(nrow(Y_rep), mean = 1/u_samples, sd = 0.045)
}
ppc_ecdf_overlay(c.data$Y, Y_rep)
ppc_dens_overlay(c.data$Y, Y_rep[1:25, ])


#######################################################
###### (5) Inverse Gamma ########
IG_code <- "data {
int<lower=0> J;
vector[J] Y;
}
parameters {
real<lower=0> u[J];
}
model {
for (i in 1:J){
u[i] ~ inv_gamma(1, 1/(Y[i]+0.01));
Y[i] ~ normal(1/u[i], 0.045);
}
}"


IG.mod <- stan_model(model_code = IG_code)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
IG.fit <- sampling(IG.mod, data=c.data, iter=5000, chains=1, verbose=TRUE, warmup=2500, init=0, seed=1, refresh=10)
IG.array <- rstan::extract(IG.fit, permuted = TRUE)
IG.u.df <- as.data.frame(IG.array[["u"]])
IG.u.summary <- summarise_draws(IG.u.df) 
indivmean.dat.IG <- data.frame(value=IG.u.summary$median, upper = IG.u.summary$q95, lower = IG.u.summary$q5,method = "IG (4,1)")
indivmean.dat.IG <- cbind(indivmean.dat.IG, "Y" = y)


IG.data <- indivmean.dat.IG[-1,]
IG.data$rtrue <- 1 / IG.data$Y
se_IG <- (IG.data$value - IG.data$rtrue)^2
bias_IG <- IG.data$value - IG.data$rtrue

posterior_samples <- extract(IG.fit)
J <- c.data$J
Y_rep <- matrix(NA, nrow = nrow(posterior_samples$u), ncol = J)

for (i in 1:J) {
  u_samples <- posterior_samples$u[, i]
  Y_rep[, i] <- rnorm(nrow(Y_rep), mean = 1/u_samples, sd = 0.045)
}

ppc_ecdf_overlay(c.data$Y, Y_rep)
ppc_dens_overlay(c.data$Y, Y_rep[1:25, ])

################################################################
########## (6) Product Half-Cauchy #############
HS_code <- "
data {
  int<lower=0> J;
  vector[J] Y;
}
parameters {
  real<lower=0> tau[J];
  real<lower=0> lambda[J];
}
transformed parameters {
  real<lower=0> u[J];
  for (i in 1:J) {
    u[i] = tau[i] * lambda[i];
  }
}
model {
  for (i in 1:J){
    // Horseshoe prior
    tau[i] ~ cauchy(0, 1);
    lambda[i] ~ cauchy(0, 1);    
    Y[i] ~ normal(1 / u[i], 0.045);
  }
}
"

HS.mod <- stan_model(model_code = HS_code)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
HS.fit <- sampling(HS.mod, data=c.data, iter=5000, chains=1, verbose=TRUE, warmup=2500, init=0, seed=123, refresh=10)
HS.array <- rstan::extract(HS.fit, permuted = TRUE)
HS.u.df <- as.data.frame(HS.array[["u"]])
HS.u.summary <- summarise_draws(HS.u.df) 
indivmean.dat.HS <- data.frame(value=HS.u.summary$median, upper = HS.u.summary$q95, lower = HS.u.summary$q5,method = "HS")
indivmean.dat.HS <- cbind(indivmean.dat.HS, "Y" = y)


HS.data <- indivmean.dat.HS[-1,]
HS.data$rtrue <- 1 / HS.data$Y

se_HS <- (HS.data$value - HS.data$rtrue)^2
bias_HS <- HS.data$value - HS.data$rtrue



Cauchy_mse <- data.frame(rtrue = 1 / Cauchy.data$Y, MSE = se_cauchy, Method = "Half-Cauchy")
RG_mse <- data.frame(rtrue = 1 / RG.data$Y, MSE = se_RG, Method = "Reciprocal-Gaussian")
Weibull_mse <- data.frame(rtrue = 1 / Weibull.data$Y, MSE = se_Weibull, Method = "Weibull")
Gamma_mse <- data.frame(rtrue = 1 / Gamma.data$Y, MSE = se_Gamma, Method = "Gamma")
IG_mse <- data.frame(rtrue = 1 / IG.data$Y, MSE = se_IG, Method = "Inverse-Gamma")
HS_mse <- data.frame(rtrue = 1 / HS.data$Y, MSE = se_HS, Method = "Product Half-Cauchy")

combined_mse_data <- rbind(Cauchy_mse,Gamma_mse,IG_mse,RG_mse,Weibull_mse,HS_mse)
combined_mse_data$True_Fractional_Parallax <- 0.045 * combined_mse_data$rtrue

mse_plot <- ggplot(combined_mse_data, aes(x = True_Fractional_Parallax, y = MSE, color = Method)) +
  geom_point(alpha = 0.6) +  
  geom_line() +              
  scale_color_manual(values = c("Gamma" = "firebrick", "Inverse-Gamma" = "deepskyblue", "Half-Cauchy" = "gold","Reciprocal-Gaussian" = "orange","Weibull"="purple","Product Half-Cauchy"="green")) +
  labs(title = "",
       x = "True Fractional Parallax",
       y = "Squared Error",
       color = "Priors") + 
  theme_minimal(base_size = 16) +
  theme(
    legend.position  = "right",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(
    override.aes = list(
      linetype = 1,
      size     = 1.8,
      alpha    = 1,
      shape    = 16
    )
  ))

print(mse_plot)
# ggsave("MSE_sim1.png", plot = mse_plot, units = "in", width = 6, height = 5, dpi=400)
# dev.off()

##################### Codes for simulating Setup B #######################
##########################################################################
set.seed(42)
n_sims   <- 500
r_lim    <- 1e4
f_true_grid <- seq(0.6, 2.0, by = 0.1)

r_grid <- seq(1, 2 * r_lim, length.out = 5000)
delta  <- diff(r_grid)[1]  

prior_gamma <- function(r) {
  dgamma(r, shape = 3, scale = 1e3)
}

prior_inv_gamma <- function(r) {
  dinvgamma(r, shape = 3, scale = 1e3)
}

prior_half_cauchy <- function(r) {
  2 * dcauchy(r, location = 0, scale = 1e3) * (r >= 0)
}

prior_rg <- function(r) {
  lam <- 1e3
  sqrt(2 / pi) * (1 / (lam * r^2)) * exp(-1 / (2 * (lam^2) * r^2))
}

prior_phc <- function(r) {
  scale <- 1e3
  v <- log(r / scale)
  v[v < 0] <- 0
  v / r^2
}

prior_weibull <- function(r) {
  dweibull(r, shape = 0.8, scale = 1e3)
}

priors <- list(
  Gamma     = prior_gamma,
  `Inverse-Gamma` = prior_inv_gamma,
  `Half-Cauchy` = prior_half_cauchy,
  `Reciprocal-Gaussian`         = prior_rg,
  `Product Half-Cauchy`        = prior_phc,
  `Weibull`    = prior_weibull
)

rmse_results <- matrix(NA, 
                       nrow = length(f_true_grid), 
                       ncol = length(priors),
                       dimnames = list(NULL, names(priors)))

for (f_idx in seq_along(f_true_grid)) {
  f_true <- f_true_grid[f_idx]
  
  u0      <- runif(n_sims)
  r_true  <- r_lim * u0^(1/3)
  pi_true <- 1 / r_true
  sigma   <- f_true * pi_true
  
  pi_obs  <- rnorm(n_sims, mean = pi_true, sd = sigma)
  
  for (p_idx in seq_along(priors)) {
    pfun <- priors[[p_idx]]
    r_est <- numeric(n_sims)
    
    for (i in seq_len(n_sims)) {
      # Unnormalized posterior
      lik  <- dnorm(pi_obs[i], mean = 1 / r_grid, sd = sigma[i])
      post_unnorm <- pfun(r_grid) * lik
      post <- post_unnorm / trapz(r_grid, post_unnorm)
      
      cdf  <- cumsum(post) * delta
      cdf  <- cdf / max(cdf)
      
      # Posterior median
      r_est[i] <- approx(cdf, r_grid, xout = 0.5)$y
    }
    
    # Compute RMSE
    #rmse_results[f_idx, p_idx] <- sqrt(mean((r_est - r_true)^2))
    rmse_results[f_idx, p_idx] <- mean((1-(r_est/r_true))^2)
  }
}


# Plot RMSE vs f_true
df <- data.frame(f_true = f_true_grid, rmse_results)
df_long <- pivot_longer(df, cols = -f_true,
                        names_to  = "Prior",
                        values_to = "RMSE")

ggplot(df_long, aes(x = f_true, y = RMSE, color = Prior)) +
  geom_point(alpha = 0.4) +               
  geom_smooth(method = "loess", se = FALSE, span = 0.5) +
  labs(
    x     = expression(f[true]),
    y     = "RMSE",
    title = "Comparing RMSE for different priors"
  ) +
  theme_minimal()
colnames(rmse_results) <- c("Gamma","Inverse-Gamma","Half-Cauchy","Reciprocal-Gaussian","Product Half-Cauchy","Weibull")

df <- data.frame(f_true = f_true_grid,
                 rmse_results,
                 check.names = FALSE)
df_long <- pivot_longer(df, cols = -f_true,
                        names_to  = "Prior",
                        values_to = "RMSE")

df_long$Prior <- factor(df_long$Prior,
                        levels = c("Gamma",
                                   "Half-Cauchy",
                                   "Inverse-Gamma",
                                   "Product Half-Cauchy",
                                   "Reciprocal-Gaussian",
                                   "Weibull")
)

my_cols <- c(
  "Gamma"               = "#B22222",
  "Half-Cauchy"         = "#FFD700",
  "Inverse-Gamma"       = "#1E90FF",
  "Product Half-Cauchy" = "#228B22",
  "Reciprocal-Gaussian" = "#FF8C00",
  "Weibull"             = "#8A2BE2"
)

p_re <- ggplot(df_long, aes(x = f_true, y = RMSE, color = Prior)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.5, linewidth = 1) +
  scale_color_manual(values = my_cols) +
  labs(
    x     = "True Fractional Parallax",
    y     = "Mean Squared Relative Error",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(
    override.aes = list(
      linetype = 1,
      size     = 1.5,
      alpha    = 1,
      shape    = 16
    )
  ))

ggsave("RMSE_RE.png", plot = p_re, units = "in", width = 6, height = 4, dpi=400)
dev.off()

plot(f_true_grid, rmse_results[,1], type = "n",
     xlab = expression(f[true]), ylab = "RMSE",
     main = "",
     ylim = range(rmse_results))

cols <- rainbow(ncol(rmse_results))
legend("topleft", legend = colnames(rmse_results), col = cols, lwd = 2)

# for each prior, fit a smoothing spline & add
for (j in seq_len(ncol(rmse_results))) {
  ss <- smooth.spline(f_true_grid, rmse_results[,j], spar = 0.3)
  lines(ss, col = cols[j], lwd = 2)
}

library(stats)

set.seed(0)
J        <- 2000
f_values <- seq(0.01, 1, length.out = 20)
r0       <- 2000

# Prior definitions
r_lim <- 1000
L     <- 1000
prior_list <- list(
  `Improper Uniform` = function(r) ifelse(r > 0, 1, 0),
  `Proper Uniform`   = function(r) ifelse(r > 0 & r <= r_lim, 1/r_lim, 0),
  `Constant-Volume`  = function(r) ifelse(r > 0 & r <= r_lim, r^2, 0),
  `ED (Gamma)` = function(r) r^2 * exp(-r / 2000)
)

posterior_mode <- function(omega, sigma, prior_pdf) {
  neg_logp <- function(r) {
    if (r <= 0) return(Inf)
    ll <- 0.5 * ((omega - 1/r)^2) / (sigma^2)
    lp <- -log(prior_pdf(r) + 1e-300)
    ll + lp
  }
  opt <- optimize(neg_logp, interval = c(1e-6, 5000))
  return(opt$minimum)
}

bias <- lapply(prior_list, function(.) numeric(length(f_values)))
sd   <- lapply(prior_list, function(.) numeric(length(f_values)))

for (i in seq_along(f_values)) {
  f     <- f_values[i]
  sigma <- f * (1/r0)
  omegas <- rnorm(J, mean = 1/r0, sd = sigma)
  
  for (nm in names(prior_list)) {
    prior_pdf <- prior_list[[nm]]
    modes <- vapply(omegas,
                    FUN = function(w) posterior_mode(w, sigma, prior_pdf),
                    FUN.VALUE = 0.0)
    bias[[nm]][i] <- mean(modes - r0)
    sd[[nm]][i]   <- sd(modes)
  }
}

op <- par(mfrow = c(1,2), mar = c(4.5,5,2,1))  

abs_bias_vals <- lapply(bias, abs)
ylim_bias <- range(unlist(abs_bias_vals))
plot(f_values, abs_bias_vals[[1]], type = "o", col = 1, ylim = ylim_bias,
     xlab = "Fractional error f", ylab = "Absolute bias of posterior mode",
     cex.lab = 1.4,    
     cex.axis = 1.2)    
for (k in 2:length(abs_bias_vals)) {
  lines(f_values, abs_bias_vals[[k]], type = "o", col = k)
}
legend("topleft", legend = names(bias), col = seq_along(bias),
       lty = 1, pch = 1,
       cex = 1.0)        

# SD
ylim_sd <- range(unlist(sd))
plot(f_values, sd[[1]], type = "o", col = 1, ylim = ylim_sd,
     xlab = "Fractional error f", ylab = "SD of posterior mode",
     cex.lab = 1.4,
     cex.axis = 1.2)
for (k in 2:length(sd)) {
  lines(f_values, sd[[k]], type = "o", col = k)
}
legend("topleft", legend = names(sd), col = seq_along(sd),
       lty = 1, pch = 1,
       cex = 1.0)

par(op)



df_sd <- data.frame(
  f = f_values,
  `Improper Uniform`  = sd[[ "Improper Uniform" ]],
  `Proper Uniform`    = sd[[ "Proper Uniform"   ]],
  `Constant-Volume`   = sd[[ "Constant-Volume" ]],
  `ED (Gamma)`                = sd[[ "ED (Gamma)" ]], check.names = FALSE
)

df_long <- df_sd %>%
  pivot_longer(
    cols = -f,
    names_to  = "Prior",
    values_to = "SD"
  )

init_plot <- ggplot(df_long, aes(x = f, y = SD, color = Prior)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x      = "Fractional error f",
    y      = "SD of posterior mode",
    color  = "Prior"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
# ggsave("init_priors.png", plot = init_plot, units = "in", width = 6, height = 5, dpi=400)
# dev.off()





