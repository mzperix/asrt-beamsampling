library(readr)
library(ggplot2)

# params <- read_csv(here::here("Data/artificial_asrt_params.csv"))

rnormt <- function(n, range, mu, s = 1) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
  
}

sample_rt <- function(n, tau0, mu, sigma, p){
  r <- rnormt(n, c((tau0-log(p))/5000, 1000), mu, sigma)
  rt <- (tau0-log(p)) / r
  return(rt)
}

fig_prior_predictive_checks <- function(params){
  ps <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  n <- 1000
  samples <- data.frame()
  for (row in 1:nrow(params)){
    for (p in ps){
      tau0 <- params[[row, "tau0"]]
      mu <- params[[row, "mu"]]
      sigma <- params[[row, "sigma"]]
      samples <- rbind(samples, data.frame(id=row, rt=sample_rt(1000, tau0, mu, sigma, p), p=p, mu=mu, sigma=sigma, tau0=tau0))  
    }
  }

  samples %>%
    ggplot(aes(x=p, y=rt, group=p)) +
    geom_violin() +
    facet_wrap("id", nrow=9) +
    ylim(0,3000) +
    style
}