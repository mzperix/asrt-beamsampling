data{
    int<lower=1> K; // number of states currently represented
    simplex[K+1] pi[K]; // transition matrix
    
    // Prior parameters
    real<lower=0.0> alpha0_a;
    real<lower=0.0> alpha0_b;
    real<lower=0.0> gamma_a;
    real<lower=0.0> gamma_b;    
}

parameters{
    real<lower=0.0> alpha0;
    real<lower=0.0> gamma;
    simplex[K+1] beta;
}

model{
    real beta_sum;
    
    // Priors
    alpha0 ~ gamma(alpha0_a, alpha0_b);
    gamma ~ gamma(gamma_a, gamma_b);
    
    // Stick-breaking (backwards) for beta
    beta_sum = beta[K+1];
    for (k in K:1){
        target += beta_lpdf(beta[k]/(beta[k]+beta_sum) | 1 , gamma);
        beta_sum += beta[k];
    }
    
    // Likelihood of pi
    for (k in 1:K){
        pi[k] ~ dirichlet(alpha0*beta);
    }
}