data{
    int<lower=1> T; // number of trials

    int<lower=1> D; // number of possible symbols
    int<lower=1> K; // number of represented states in current sample

    int Y[T];  // stimuli

    matrix[K,D] M; // ~sufficient stats for earlier experiences
    matrix[1,D] H; // dirichlet prior for emissions

    real<lower=0.0> alpha0; // prior parameter for pi
    // simplex[K] beta;      // prior parameter for pi

    real rt[T];     // reaction times
    int<lower=0, upper=1> correct_response[T]; // whether responses are correct (1) or not (0)
    int<lower=0, upper=1> filters[T]; // other filters: 0: leave out from conditioning, 1: use
    int<lower=0, upper=1> new_stream[T]; // 1: a new stream starts, that is, state posterior gets refreshed, 0: no new stream on this trial

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

    simplex[D] phi[K]; // emission "matrix"
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
    vector[K] s_hat[T];
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

    // Prior for phi
    for (i in 1:K){
        for (j in 1:D){
            target += phi[i][j]^(H[1,j]+M[i,j]-1);
        }
    }

    // Set up participant's posterior
    s_hat[1] = pi_matrix[1,:]';

    // Condition and propagate
    for (t in 2:T){
        if (new_stream[t] == 1){
            s_hat[t] = pi_matrix[1,:]';
        }
        else{
            // Condition on observation
            for (k in 1:K){
                s_hat[t][k] = s_hat[t-1][k]*phi[k][Y[t-1]];
            }

            // Use pi to propagate belief
            // Only values of pi that are above epsilon are used
            //s_hat[t] =  multiply(softplus_epsilon(pi_matrix', beta_softplus, epsilon),s_hat[t]);
            s_hat[t] =  multiply(pi_matrix',s_hat[t]);
            s_hat[t] /= sum(s_hat[t]);
        }
    }

    // Reaction times
    for (t in 1:T){
        p_pred = 0;
        for (k in 1:K){
            p_pred += s_hat[t][k]*phi[k][Y[t]];
        }
        if (correct_response[t] == 1){
            if (filters[t] == 1){
                rt_reciprocal[t] ~ normal(mu/(tau0-log(p_pred)), sigma/(tau0-log(p_pred))) T[1/rt_max,];
            }
        }
    }

}
