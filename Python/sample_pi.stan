functions {
    matrix softplus_epsilon(matrix x, real beta, real epsilon) {
        return (log1p_exp(beta*(x-epsilon))/beta+inv_logit((x-epsilon)*beta)*epsilon);
    }
    
    real log_max_epsilon(real x, real beta, real epsilon){
        // Soft 'enforcer' of x < epsilon. 
        // To be added to target
        return (-beta*(x-epsilon)-log1p_exp(-(beta*(x-epsilon))));
    }
}

data{
    int<lower=1> T; // number of trials
    
    int<lower=1> D; // number of possible symbols
    int<lower=1> K; // number of represented states in current sample
    
    int Y[T];  // stimuli
    real<lower=0.0> epsilon[T]; // slicing parameters
    real<lower=0.0> epsilon_min; // minimum of the slicing parameters

    simplex[D] phi[K]; // emission "matrix"
    
    matrix[K,D] M; // ~sufficient stats for earlier experiences
    
    real<lower=0.0> alpha0; // prior parameter for pi
    simplex[K+1] beta;      // prior parameter for pi
    
    real rt[T];     // reaction times

    real<lower=0.0> tau0;  // rt parameters
    real<lower=0.0> mu;   
    real<lower=0.0> sigma;
    
    real<lower=0.0> error_sigma;
    real<lower=0.0> rt_max;
    
    real<lower=0.0> beta_softmax;
}

transformed data{
    real<lower = 1/rt_max> rt_reciprocal[T];
    for (i in 1:T){
        rt_reciprocal[i] = 1/rt[i];
    }   
}

parameters{
    simplex[K+1] pi[K];  // transition "matrix"
}

transformed parameters{
    matrix[K,K] pi_matrix;

    for (k in 1:K){
        for (j in 1:K){
            pi_matrix[k,j] = pi[k][j];
        }
        
    }
}

model{
    real p_pred;
    vector[K] s_hat[T];
    
    // Set up posterior
    s_hat[1] = pi_matrix[1,:]';
    
    // Prior for pi
    for (i in 1:K){
        pi[i] ~ dirichlet(alpha0*beta); 
        
        // We condition on the 'unvisited' states to have goto probability below all the epsilons.
        target += log_max_epsilon(pi[i][K+1], beta_softmax, epsilon_min);
    }
    
    // Condition and propagate
    for (t in 2:T){
        // Condition on observation
        for (k in 1:K){
            s_hat[t][k] = s_hat[t-1][k]*phi[k][Y[t-1]]; 
        }
        
        // Use pi to propagate belief
        // Only values of pi that are above epsilon are used
        s_hat[t] =  multiply(softplus_epsilon(pi_matrix', beta_softmax, epsilon[t]),s_hat[t]); 
        s_hat[t] /= sum(s_hat[t]);
    }
    
    // Reaction times
    for (t in 1:T){
        p_pred = 0;
        for (k in 1:K){
            p_pred += s_hat[t][k]*phi[k][Y[t]];
        }
        rt_reciprocal[t] ~ normal(mu/(tau0-log(p_pred)), sigma/(tau0-log(p_pred))) T[1/rt_max,]; 
    }

}