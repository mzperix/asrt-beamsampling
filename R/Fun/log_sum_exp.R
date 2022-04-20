stanmodelcode <- "
  data {
    int<lower=0> N;
    real rt_ct[N];
    real<lower=0.0> exp_rt_markov[N];
    real<lower=0.0> exp_rt_trigram[N];
    real m[N];
    //real<lower=0.1> sigma;
  }
  parameters {
    real a;
    real b;
    real sigma;
  }
  model {
    //target += normal_lpdf(a | 0, 10);
    //target += normal_lpdf(b | 0, 10);
    a ~ gamma(2,0.5);
    b ~ gamma(2,0.5);
    sigma ~ gamma(10,0.5);
    for (n in 1:N) {
      rt_ct[n] ~ normal(m[n] + log(a*exp_rt_markov[n]+b*exp_rt_trigram[n]), sigma);
    }
}
"

data <- df %>%
  subset(participant_test==119) %>%
  subset(session_test==8)

data <- df %>%
  #mutate(session_test = factor(session_test, levels=c(1,2,3,4,5,6,7,8))) %>%
  #filter_and_average_runs(sd_cutoff=100, apply_filters = TRUE) %>%
  spread(model, pred_prob) %>%
  subset(participant_test==119) %>%
  subset(session_test==8)

N <- length(data$Markov)
rt_markov <- data$Markov
rt_trigram <- data$Triplet
m <- pmax(rt_markov, rt_trigram)
m <- rep(0,N)
exp_rt_markov <- exp(rt_markov-m)
exp_rt_trigram <- exp(rt_trigram-m)
rt_ct <- data$iHMM
N <- length(rt_ct)

dat <- list(
  N = N, 
  rt_ct = rt_ct,
  exp_rt_markov = exp_rt_markov,
  exp_rt_trigram = exp_rt_trigram,
  m = m
);
fit <- stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 10012, chains = 3, verbose = TRUE,
            sample_file = file.path(tempdir(), 'norm.csv'))
print(fit)


get_pred <- function(a, b, exp_rt_markov, exp_rt_trigram, m){
  log(a*exp_rt_markov+b*exp_rt_trigram) + m
}
df <- data.frame(
  ct=rt_ct, 
  fitted=get_pred(0.32, 0.68, exp_rt_markov, exp_rt_trigram, m)) 
df <- df %>% mutate(resid = ct-fitted)
df %>% summarise(rmse = sqrt(mean(resid^2)))

df %>%
  ggplot(aes(x=fitted, y=ct)) +
  geom_point()

a <- 0.7
b <- 0.2
N <- 500
rt_markov <- rnorm(N)*0.5+0.5 #rnorm(N)*30+250
rt_trigram <- rnorm(N)*0.5+0.5 #rnorm(N)*34+400
m <- pmax(rt_markov, rt_trigram)
exp_rt_markov <- exp(rt_markov-m)
exp_rt_trigram <- exp(rt_trigram-m)
m <- rep(0,N)
rt_ct <- log(a*exp_rt_markov+b*exp_rt_trigram) + m + rnorm(N)*0.2
dat <- list(
  N = N, 
  rt_ct = rt_ct,
  exp_rt_markov = exp_rt_markov,
  exp_rt_trigram = exp_rt_trigram,
  m = m
  );
fit <- stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 10012, chains = 3, verbose = TRUE,
            sample_file = file.path(tempdir(), 'norm.csv'))
print(fit)
