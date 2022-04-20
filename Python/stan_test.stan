data{
        int<lower=1> N;
        real x[N];
        
        real<lower=0.0> mu_mu;
        real<lower=0.0> mu_sigma;
        real<lower=0.0> sigma_mu;
        real<lower=0.0> sigma_sigma;
        
        real<lower=0.0> error_sigma;
        real<lower=0.0> x_max;
    }
    parameters{
        real<lower = 0.0> mu;
        real<lower = 0.0> sigma;
        real<lower = 0.0> r[N];
    }
    model{
        mu ~ normal(mu_mu, mu_sigma);
        sigma ~ normal(sigma_mu, sigma_sigma);
        for (i in 1:N){
            r[i] ~ normal(mu, sigma) T[0,1/x_max];
            x[i] ~ normal(1/r[i], error_sigma) T[0,x_max];
        }
    }