//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
int<lower=0> J;
vector[J] Y;
}
parameters {
real<lower=0> u[J];
}
model {
for (i in 1:J){
u[i] ~ gamma(3, 1000);
Y[i] ~ normal(1/u[i], 1);
}
} 