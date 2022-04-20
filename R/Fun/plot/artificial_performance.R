setwd('~/Science/asrt-beamsampling/')
source('R/Fun/plot/plot_styles.R')
library(ggplot2)
library(magrittr)
library(dplyr)

artificial_generation_data <- read.csv('Data/artificial_asrt_params.csv') %>%
  mutate(participant_test=0:80,
         model=factor(internal_model_steps, labels=c("Early","Middle","Late")),
         noise_level=tau0*sigma/mu)

performance_data <- read.csv('Python/Output/832a60f/84cabc2_PERFORMANCE_artificial_asrt.csv')
perf_data <- performance_data %>%
  group_by(participant_test) %>%
  summarise(rt_correlation = mean(correlation)) %>%
  ungroup()

pred_prob_data <- read.csv('Python/Output/832a60f/84cabc2_PRED_PROBS_artificial_asrt.csv')
log_pred_probs <- pred_prob_data %>%
  group_by(block, trial, participant_test) %>%
  summarise(mean_log_pred_prob = mean(LOG_PRED_PROB),
            y0 = mean(y0),
            y1 = mean(y1),
            y2 = mean(y2),
            y3 = mean(y3),
            rt = mean(rt)) %>%
  ungroup()

d <- data.frame()
for (i in 0:80){
  d <- rbind(d, read.csv(paste0('Data/artificial_asrt/artificial_asrt_ASRT_',i,'_pred_prob_all.csv')))
}
d <- d %>%
  group_by(participant_test) %>%
  mutate(n=1:n())

log_pred_probs <- log_pred_probs %>%
  group_by(participant_test) %>%
  mutate(n=1:n())

pred_corr <- merge(log_pred_probs, d, by=c("participant_test", "n", "block", "trial")) %>%
  group_by(participant_test) %>%
  summarise(#c = cor(mean_log_pred_prob, log_pred_prob),
            #cross_entropy = sum(-exp(log_pred_prob)*mean_log_pred_prob),
            cross_entropy=-sum(y0.x*log(y0.y+1e-16)+y1.x*log(y1.y+1e-16)+y2.x*log(y2.y+1e-16)+y3.x*log(y3.y+1e-16)),
            entropy=-sum(y0.x*log(y0.x+1e-16)+y1.x*log(y1.x+1e-16)+y2.x*log(y2.x+1e-16)+y3.x*log(y3.x+1e-16)),
            entropy_difference = (cross_entropy-entropy)/n(),
            original_rt_correlation = cor(rt.x, original_rt_predicted),
            rt_var = var(rt.x))

plot_data <- merge(pred_corr, perf_data) %>%
  merge(artificial_generation_data)

plot_data %>% 
  ggplot(aes(x=rt_correlation^2, y=entropy_difference, shape=model, color=sqrt(rt_var))) +
  geom_point(size=2.5) +
  xlab('Response time prediction performance') +
  ylab('Recovery performance using cross-entropy and entropy') +
  #xlim(0,1) +
  #ylim(0,1) +
  #scale_color_viridis_d() +
  scale_color_viridis_c() +
  style +
  #coord_fixed() +
  guides(shape=guide_legend(title='Model'),
         color=guide_colourbar(title='RT sd. (msec)'))

plot_data %>% 
  ggplot(aes(x=c^2, y=entropy_difference, shape=model, color=sqrt(rt_var))) +
  geom_point(size=2.5) +
  xlab('Correlation') +
  ylab('Cross-entropy evaluation') +
  #xlim(0,1) +
  #ylim(0,1) +
  #scale_color_viridis_d() +
  scale_color_viridis_c() +
  style +
  #coord_fixed() +
  guides(shape=guide_legend(title='Model'),
         color=guide_colourbar(title='RT sd. (msec)'))


plot_data %>%
  ggplot(aes(x=original_rt_correlation^2, y=rt_correlation^2, shape=model, color=sqrt(rt_var))) +
  geom_abline(slope=1.0, intercept=0, color="#cccccc") +
  geom_point(size=2.5) +
  xlab('R^2 using actual model') +
  ylab('R^2 with recovered model') +
  xlim(0,1) +
  ylim(0,1) +
  #scale_color_viridis_d() +
  scale_color_viridis_c() +
  style +
  coord_fixed(expand=FALSE) +
  guides(shape="none", color="none")
#ggsave('Figures/synthetic_performance.pdf', width=5, height=4.3)

