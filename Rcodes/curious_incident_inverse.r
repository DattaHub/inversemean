######################### vectorized
rm(list=ls())
library(plyr)
library(rstan)
library(parallel)
# library(rbenchmark)
library(reshape2)
library(here)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Use CMDSTAN 

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

cauchy.code <- file.path(here("cauchy_code.stan"))
cauchy.mod <- cmdstan_model(cauchy.code)

J = 200
y= seq(0, 30, length.out = J)


c.data = list('J'=J, 'Y' = y)
stan.iters = 5000

cauchy.fit <- cauchy.mod$sample(
  data = c.data, seed = 123,  chains = 1,
  parallel_chains = 1,  refresh = 10,
  iter_sampling = stan.iters, iter_warmup = stan.iters/2,
  init = 0)

library(posterior)

cauchy.array <- cauchy.fit$draws("u")
cauchy.u.df <- as_draws_df(cauchy.array) # as_draws_matrix() for matrix
cauchy.u.summary <- summarise_draws(cauchy.u.df) ## you can also use custome set of summary like ("mean", "median", Mode, "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")

# cauchy.u.summary$mean

plot(y, cauchy.u.summary$mean)

# cauchy.u.smpls = rstan::extract(cauchy.res, pars=c("u"), permuted=TRUE)[[1]]
# cauchy.u.mean = apply(cauchy.u.smpls,2,mean)


gamma.code <- file.path(here("gamma_code.stan"))
gamma.mod <- cmdstan_model(gamma.code)

gamma.fit <- gamma.mod$sample(
  data = c.data, seed = 123,  chains = 1,
  parallel_chains = 1,  refresh = 10,
  iter_sampling = stan.iters, iter_warmup = stan.iters/2,  init = 0)

library(posterior)

gamma.array <- gamma.fit$draws("u")
gamma.u.df <- as_draws_df(gamma.array) # as_draws_matrix() for matrix
gamma.u.summary <- summarise_draws(gamma.u.df) ## you can also use custome set of summary like ("mean", "median", Mode, "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")

# gamma.u.summary$mean

plot(y, gamma.u.summary$mean)

invmean.dat <- rbind(data.frame(value = gamma.u.summary$mean, upper = gamma.u.summary$q95, lower = gamma.u.summary$q5, method = "Gamma"),
                     data.frame(value = cauchy.u.summary$mean, upper = cauchy.u.summary$q95, lower = cauchy.u.summary$q5,method = "Cauchy"))


invmean.dat <- cbind(invmean.dat, "Y" = y)

library(ggplot2)

ggplot(invmean.dat, aes(x = Y, y = value, fill = method)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper, group = method),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick", size = 1) + facet_wrap(.~method, ncol = 1, scales = "free_y") + 
  theme_bw()

# ggplot(invmean.dat, aes(x = U, y = value, fill = method)) +
#   geom_ribbon(aes(ymin = lower,
#                   ymax = upper, group = method),    # shadowing cnf intervals
#               fill = "steelblue2") + 
#   geom_line(color = "firebrick", size = 1) + facet_wrap(.~method, scales = "free_y") + 
#   theme_bw()


# gamma.res = sampling(gamma.fit, 
#                       data = c.data, 
#                       iter = stan.iters,
#                       warmup = floor(stan.iters/2),
#                       thin = 2,
#                       pars = c('u'),
#                       init = 0,
#                       seed = seed.val, 
#                       chains = 1)
# 
# gamma.u.smpls = rstan::extract(gamma.res, pars=c("u"), permuted=TRUE)[[1]]
# gamma.u.mean = apply(gamma.u.smpls,2,mean)


