
functions {
 real custom_exp_likelihood(int d, int x_, int t, real r, real b) {
   real ll;
   real lambda;
   lambda = exp(r + b * x_);
   ll = d * log(lambda) - lambda * t;
   return ll;
 }
}

data {
  int<lower=0> NR;
  int<lower=0> J;
  int x[NR];
  int status[NR];
  int interval[NR];
  int interval_fu_time[NR];
}

parameters {
  real interval_beta[J];
  real x_beta;
  real<lower=0> baseline_sigma;
}

model {
    x_beta ~ normal(0, 100);
    baseline_sigma ~ normal(0, 100);
  //  interval_beta ~ normal(0, baseline_sigma); 

  for (i in 1:NR)
    target += custom_exp_likelihood(status[i], x[i], interval_fu_time[i], interval_beta[interval[i]], x_beta);
}
