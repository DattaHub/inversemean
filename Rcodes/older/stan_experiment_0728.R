## x ~ N(u^{-1}, 1); u ~ C(0,1) vs u ~ G(3,10^3)
## Plot posterior mean E(u |x ) for a range of x in (0,30)

## This is the Stan code used for SS & Max simulation studies
library(rstan)
#set_cppo("fast")
library(ggplot2)
library(plyr)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
########## HS Code ############
cauchy.code = "
data {
real Y; 
}
parameters {
real<lower=0> u[J];
}
model {
u ~ cauchy(0, 1);
Y ~ normal(1/u, 1);
}  
"
cauchy.fit = stan_model(model_code=cauchy.code, model_name="Cauchy")
J = 100
y = seq(0,30,length.out = J)
u.mean = rep(0,J)
for (i in 1:J)
{
  c.data = list('Y' = y[i])
  stan.iters = 1000
  seed.val = 495
  
  cauchy.res = sampling(cauchy.fit, 
                        data = c.data, 
                        iter = stan.iters,
                        warmup = floor(stan.iters/2),
                        thin = 2,
                        pars = c('u'),
                        init = 0,
                        seed = seed.val, 
                        chains = 1)
  
  u.smpls = rstan::extract(cauchy.res, pars=c("u"), permuted=TRUE)[[1]]
  u.mean[i] = mean(u.smpls)
  
}

