data{
    int<lower=1> T; // number of trials
    
    int<lower=1> D; // number of possible symbols
    int<lower=1> K; // number of represented states in current sample
    int<lower=1> N; // number of particles
    
    int Y[T];  // stimuli
    real U[T]; // particle slicing parameters
    
    simplex[K] Pi[K];  // transition "matrix"
    
    matrix[1,D] H; // dirichlet prior for emissions
    matrix[K,D] M; // ~sufficient stats for earlier experiences
    
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
    
    simplex[D] phi[K];
    //matrix<lower=0.0, upper=1.0>[N,T] u_unnormed;
}

transformed parameters{
    real<lower = 0.0> tau0_resc;
    real<lower = 0.0> mu_resc;
    real<lower = 0.0> sigma_resc;
    //matrix[N,T] u;
        
    tau0_resc = tau0 / tau0_scale;
    mu_resc = mu / mu_scale;
    sigma_resc = sigma / sigma_scale;
    
    //for (t in 1:T){
    //    u[:,t] = u_unnormed[:,t] * U[t];
    //}
}

model{
    real p_pred;
    int<upper=K> s_hat[K,T];
    //matrix[N,T] w;
    mu_resc ~ gamma(mu_shape, 1);
    sigma_resc ~ gamma(sigma_shape, 1);
    tau0_resc ~ gamma(tau0_shape, 1);

    // u and U
    //for (t in 1:T){
    //    u_unnormed[:,t] ~ uniform(0,1);
    //    U[t] ~ uniform(max(u[:,T]),1);
    //}
    
    // PARTICLES
    // set up particles
    for (n in 1:N){
        //w[n,1] = 1/n;
        s[n,1] ~ categorical(Pi[1]);
    }
    
    // propagate particles
    for (t in 2:T){
        for (n in 1:N){ 
            s[n,t] ~ categorical(Pi[s[n,t-1]]); 
            
            // update weights using observation
            //w[n,t] = w[n,t-1]*Pi[s[n,t-1]][s[n,t]]*phi[s[n,t-1]][Y[t-1]]; 
        }
        for (n in 1:N){
            //w[n,t] /= sum(w[:,t]); // renormalize
        }
    }
    
    // reaction times
    for (t in 1:T){
        p_pred = 0;
        for (n in 1:N){
            //p_pred += phi[s[n,t]][Y[t]]*w[n,t];
        }
    //    rt_reciprocal[t] ~ normal(mu/(tau0-log(p_pred)), sigma/(tau0-log(p_pred))) T[1/rt_max,]; 
    }

    // prior for phi
    for (i in 1:K){
        for (j in 1:D){
            //target += phi[i][j]^(H[1,j]+M[i,j]-1);
        }
    }
}