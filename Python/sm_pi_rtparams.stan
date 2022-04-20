data{
    int<lower=1> T; // number of trials
    
    int<lower=1> D; // number of possible symbols
    int<lower=1> K; // number of represented states in current sample
    
    int Y[T];  // stimuli
    
    // matrix[K,D] M; // ~sufficient stats for earlier experiences
    // matrix[1,D] H; // dirichlet prior for emissions
    // simplex[D] phi[K]; // matrix of emissions

    real<lower=0.0> alpha0; // prior parameter for pi
    // simplex[K] beta;      // prior parameter for pi
    
    real rt[T];     // reaction times
    int<lower=0, upper=1> correct_response[T]; // whether responses are correct (1) or not (0)
    int<lower=0, upper=1> filters[T]; // other filters: 0: leave out from conditioning, 1: use
    
    // Prior for rt parameters
    real<lower=0.0> tau0_shape;  // prior parameters for tau
    real<lower=0.0> tau0_scale;  
    real<lower=0.0> mu_shape;    // prior parameters for mu
    real<lower=0.0> mu_scale;
    real<lower=0.0> sigma_shape; // prior parameters for sigma
    real<lower=0.0> sigma_scale;

    real<lower=0.0> error_sigma;
    real<lower=0.0> rt_max;
}

transformed data{
    real<lower = 1/rt_max> rt_reciprocal[T];
    for (i in 1:T){
        rt_reciprocal[i] = 1/rt[i];
    }   
}

parameters{
    real<lower = 0.0> tau0;
    real<lower = 0.0> mu;
    real<lower = 0.0> sigma;
    
    // emission "matrix" moved to data section
    simplex[K] pi[K];  // transition "matrix"
}

transformed parameters{
    real<lower = 0.0> tau0_resc;
    real<lower = 0.0> mu_resc;
    real<lower = 0.0> sigma_resc;
    
    matrix[K,K] pi_matrix;

    tau0_resc = tau0 / tau0_scale;
    mu_resc = mu / mu_scale;
    sigma_resc = sigma / sigma_scale;

    for (k in 1:K){
        for (j in 1:K){
            pi_matrix[k,j] = pi[k][j];
        }
        
    }
}

model{
    real p_pred;
    vector[K] pi_prior;
    
    // Prior for rt parameters
    mu_resc ~ gamma(mu_shape, 1);
    sigma_resc ~ gamma(sigma_shape, 1);
    tau0_resc ~ gamma(tau0_shape, 1);
    
    // Prior for pi
    pi_prior = rep_vector(alpha0/K, K);
    for (i in 1:K){
        pi[i] ~ dirichlet(pi_prior);
    }
    
    // Reaction times
    for (t in 2:T){
        p_pred = pi[Y[t-1],Y[t]];
        if (correct_response[t] == 1){
            if (filters[t] == 1){
                rt_reciprocal[t] ~ normal(mu/(tau0-log(p_pred)), sigma/(tau0-log(p_pred))) T[1/rt_max,]; 
            }
        }
    }

}