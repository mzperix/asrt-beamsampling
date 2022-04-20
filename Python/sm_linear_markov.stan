data{
    int<lower=1> T; // number of trials
    
    int<lower=1> D; // number of possible symbols
    
    int Y[T];  // stimuli
    
    real rt[T];     // reaction times
    int<lower=0, upper=1> correct_response[T]; // whether responses are correct (1) or not (0)
    int<lower=0, upper=1> filters[T]; // other filters: 0: leave out from conditioning, 1: use
    
    // Prior for rt parameters
    real<lower=0.0> mu_mean;    // prior parameters for mu
    real<lower=0.0> mu_sd;
    real<lower=0.0> sigma_mean; // prior parameters for sigma
    real<lower=0.0> sigma_sd;
    
    real<lower=0.0> pi_sd;
}

transformed data{
     
}

parameters{
    real<lower = 0.0> mu;
    real<lower = 0.0> sigma;
    
    // emission "matrix" moved to data section
    matrix[D,D-1] pi;  // transition "matrix"
}

transformed parameters{
    matrix[D,D] pi_matrix;

    for (k in 1:D){
        for (j in 1:(D-1)){
            pi_matrix[k,j] = pi[k][j];
        }
        pi_matrix[k,D] = -sum(pi[k,:]);
    }
}

model{
    // Prior for rt parameters
    mu ~ normal(mu_mean, mu_sd) T[0,];
    sigma ~ normal(sigma_mean, sigma_sd) T[0,];
    
    // Prior for pi
    for (i in 1:D){
        for (j in 1:(D-1)){
            pi[i,j] ~ normal(0,pi_sd);
        }
    }
    
    // Reaction times
    for (t in 2:T){
        if (correct_response[t] == 1){
            if (filters[t] == 1){
                rt[t] ~ normal(mu+pi_matrix[Y[t-1],Y[t]], sigma);
            }
        }
    }

}