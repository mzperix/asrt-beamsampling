data{
    int<lower=1> N;
    real rt[N];
    
    real logp[N];

    real<lower=0.0> tau0_shape;
    real<lower=0.0> tau0_scale;
    
    real<lower=0.0> mu_shape;
    real<lower=0.0> mu_scale;
    real<lower=0.0> sigma_shape;
    real<lower=0.0> sigma_scale;

    real<lower=0.0> error_sigma;
    real<lower=0.0> rt_max;
}

parameters{
    real<lower = 0.0> tau0;
    real<lower = 0.0> mu;
    real<lower = 0.0> sigma;
}

transformed parameters{
    real<lower = 0.0> tau0_resc;
    real<lower = 0.0> mu_resc;
    real<lower = 0.0> sigma_resc;

    tau0_resc = tau0 / tau0_scale;
    mu_resc = mu / mu_scale;
    sigma_resc = sigma / sigma_scale;

}

model{
    mu_resc ~ gamma(mu_shape, 1);
    sigma_resc ~ gamma(sigma_shape, 1);
    tau0_resc ~ gamma(tau0_shape, 1);

    for (i in 1:N){
        target += normal_lpdf(1/rt[i], mu / (tau0-logp[i]), sigma / (tau0-logp[i]))
    }
}