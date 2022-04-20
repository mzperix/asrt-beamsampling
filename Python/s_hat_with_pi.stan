functions {
    vector condition_on_y(vector s, int K, matrix phi, int y) {
        vector[K] s_new;
        for (k in 1:K){
            s_new[k] = s[k]*phi[k][y];
        }
        return s_new / sum(s_new);
    }
}

data{
    int<lower=1> T; // number of trials
    
    int<lower=1> D; // number of possible symbols
    int<lower=1> K; // number of represented states in current sample
    int<lower=1> N; // number of particles
    
    int Y[T];  // stimuli
    real epsilon[T]; // slicing parameters
    
    matrix[1,D] H; // dirichlet prior for emissions
    matrix[K,D] M; // ~sufficient stats for earlier experiences
    
    real<lower=0.0> alpha0; // prior parameter for pi
    simplex[K+1] beta;      // prior parameter for pi
    
    real rt[T];     // reaction times

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
    vector[K+1]<upper=0.0> log_pi[K];  // transition "matrix"
}

transformed parameters{
    real<lower = 0.0> tau0_resc;
    real<lower = 0.0> mu_resc;
    real<lower = 0.0> sigma_resc;

    simplex[K+1] pi[K];

    tau0_resc = tau0 / tau0_scale;
    mu_resc = mu / mu_scale;
    sigma_resc = sigma / sigma_scale;
    
    for (k in 1:K){
        pi[k] = normed_exp(log_pi[k])
    }
}

model{
    real p_pred;
    vector[K] s_hat[T];
    vector[T] epsilon;
    mu_resc ~ gamma(mu_shape, 1);
    sigma_resc ~ gamma(sigma_shape, 1);
    tau0_resc ~ gamma(tau0_shape, 1);

    // Set up posterior
    s_hat[1] = Pi[1,:]';
    
    // Condition and propagate
    for (t in 2:T){
        for (k in 1:K){
            s_hat[t][k]=s_hat[t-1][k]*phi[k][Y[t-1]];
            s_hat[t][k]=
        }
        s_hat[t]=multiply(Pi',s_hat[t]);
        //epsilon ~ uniform(0,1);
        //s_hat[t]=s_hat[t]+epsilon;
        s_hat[t]/=sum(s_hat[t]);
    }
    
    // reaction times
    for (t in 1:T){
        p_pred = 0;
        for (k in 1:K){
            p_pred += s_hat[t][k]*phi[k][Y[t]];
        }
        rt_reciprocal[t] ~ normal(mu/(tau0-log(p_pred)), sigma/(tau0-log(p_pred))) T[1/rt_max,]; 
    }

    // prior for pi
    for (i in 1:K){
        pi[i] ~ dirichlet(alpha0*beta); 
    }

    // prior for phi
    for (i in 1:K){
        for (j in 1:D){
            target += phi[i][j]^(H[1,j]+M[i,j]-1);
        }
    }
}