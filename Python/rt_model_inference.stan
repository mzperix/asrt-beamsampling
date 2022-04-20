data {
  int<lower=1> N; // number of observations
	
  real<upper=0> logP[N]; // subjective probabilities for the observations
  real<lower=0> rt[N]; // reaction time measurements

  // Gamma priors for the positive parameters
  real<lower=0> tau0_shape;
  real<lower=0> tau0_scale;
  real<lower=0> mu_shape;
  real<lower=0> mu_scale;
  real<lower=0> sigma_shape;
  real<lower=0> sigma_scale;
  real<lower=0> l_shape;
  real<lower=0> l_scale;

  // Beta prior for the lapse probability
  real p_lapse_alpha;
  real p_lapse_beta;

  real<lower=0.0> measurement_error;
}

transformed data {
  
  real tau0_rate;
  real mu_rate;
  real sigma_rate;
  real l_rate;

  tau0_rate = 1 / tau0_scale;
  mu_rate = 1 / mu_scale;
  sigma_rate = 1 / sigma_scale;
  l_rate = 1 / l_scale;

}

parameters {
  real<lower=0.0> tau0;
  real<lower=0.0> mu;
  real<lower=0.0> sigma;
  real<lower=0.0> l;
  real<lower=0.0,upper=1.0> p_lapse;
  real<lower=0.0> r[N];
  real<lower=0.0> rt_lapse[N];
}

transformed parameters{
  real<lower=0.0> rt_base[N];
  for (i in 1:N)
    rt_base[i] = 1/r[i];
}

model {
  //real rt_lapse;
  //int z_lapse;

  tau0 ~ gamma(tau0_shape, tau0_rate);
  mu ~ gamma(mu_shape, mu_rate);
  sigma ~ gamma(sigma_shape, sigma_rate);
  l ~ gamma(l_shape, l_rate);
  p_lapse ~ beta(p_lapse_alpha,p_lapse_beta);
  
  for (i in 1:N){
    r[i] ~ normal(mu, sigma) T[0,];
    //z_lapse ~ bernoulli(p_lapse);
    rt_lapse[i] ~ uniform(0,l);
    rt[i] ~ normal(rt_base[i]*(tau0-logP[i]),measurement_error) T[0,]; # chronos device error
  }
}