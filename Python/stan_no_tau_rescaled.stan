data{
    int<lower=1> N;
    real x[N];

    real tau0;
    
    real<lower=0.0> mu_shape;
    real<lower=0.0> mu_scale;
    real<lower=0.0> sigma_shape;
    real<lower=0.0> sigma_scale;
    
    real<lower=0.0> error_sigma;
    real<lower=0.0> x_max;
}

parameters{
    real<lower = 0.0> mu;
    real<lower = 0.0> sigma;
    real<lower = 1/x_max> r[N];
}

transformed parameters{
    real<lower = 0.0> mu_resc;
    real<lower = 0.0> sigma_resc;
    real<lower = 1/x_max/mu_scale> r_resc[N];

    mu_resc = mu / mu_scale;
    sigma_resc = sigma / sigma_scale;
    for (i in 1:N){
        r_resc[i] = r[i] / mu_scale;
    }
    
}

model{
    mu_resc ~ gamma(mu_shape, 1);
    sigma_resc ~ gamma(mu_shape, 1);
    for (i in 1:N){
        r_resc[i] ~ normal(mu_resc, sigma_resc) T[1/x_max/mu_scale,];
        x[i] ~ normal(tau0/r[i], error_sigma) T[0,x_max];
    }
}