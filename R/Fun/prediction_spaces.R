## Are the CT predictions a linear combination of the Markov and Triplet predictions?

library(readr)
library(broom)
library(ggplot2)
source(here::here('R/Fun/stats.R'))
source(here::here('R/Fun/plot/plot_styles.R'))
#rt_prediction_data <- read_csv(here::here('Python/Output/832a60f/866066d_LEARNING_CURVES_elarasztas_pre_submission.csv'))

# error_prediction_data <- read_csv(here::here('Python/Output/832a60f/866066d_PREDICTED_PROBABILITIES_elarasztas_pre_submission.csv'))

pred_spaces_preprocess_error_prediction_data <- function(error_prediction_data){
  df <- error_prediction_data %>% 
    add_session_info(., dataset="elarasztas") %>%
    group_by(session_test, e_train, model, participant_test, correct_response, block, trial) %>%
    summarise(y0 = mean(y0), y1=mean(y1),y2=mean(y2),y3=mean(y3), pred_prob=mean(exp(LOG_PRED_PROB))) %>%
    select(session_test, e_train, model, participant_test, correct_response,block,trial,y0,y1,y2,y3) %>%
    gather(key="key", value="pred_prob", y0,y1,y2,y3)
  df[df$model=="Triplet", "pred_prob"] = exp(df[df$model=="Triplet", "pred_prob"])
  
  return(df)
}

# fit_ct_markov_trigram_log_exp_model <- function(pred_prob_individual, iter=5000){
#   stanmodelcode <- "
#     data {
#       int<lower=0> N;
#       real p_ct[N];
#       real<lower=0.0> p_markov[N];
#       real<lower=0.0> p_trigram[N];
#       //real<lower=0.1> sigma;
#     }
#     parameters {
#       real a;
#       real b;
#       real sigma;
#     }
#     model {
#       a ~ gamma(2,0.5);
#       b ~ gamma(2,0.5);
#       sigma ~ gamma(10,0.5);
#       for (n in 1:N) {
#         p_ct[n] ~ normal(a*p_markov[n]+b*p_trigram[n], sigma);
#       }
#     }
#   "
#   
#   data <- pred_prob_individual #%>%
#     #spread(model, pred_prob)
#     
#   dat <- list(
#     p_markov = data$Markov,
#     p_trigram = data$Triplet,
#     p_ct = data$iHMM,
#     N = length(data$Markov)
#   );
#   
#   print(min(dat[['p_ct']]))
#   
#   fit <- stan(model_code = stanmodelcode, model_name = "example",
#               data = dat, iter = iter, chains = 3, verbose = TRUE,
#               sample_file = file.path(tempdir(), 'norm.csv'))  
#   print(fit)
#   return(fit)
# }


# participants <- unique(error_prediction_data$participant_test)
# results <- data.frame()
# for (p in participants){
#   df <- error_prediction_data %>% 
#     pred_spaces_preprocess_error_prediction_data() %>%
#     subset(participant_test==p) %>%
#     subset(session_test==8) %>%
#     spread(model, pred_prob)
#   
#   fit <- df %>%
#     fit_ct_markov_trigram_log_exp_model(iter=5000)
#   samples <- extract(fit, c("a","b"), permuted = TRUE, inc_warmup = FALSE,
#                      include = TRUE)
#   p_result <- data.frame(session_test=8, participant_test=p, a=mean(samples$a), b=mean(samples$b))
#   
#   results <- rbind(results, p_result)
# }

# get_pred <- function(a, b, exp_rt_markov, exp_rt_trigram){
#  a*exp_rt_markov+b*exp_rt_trigram
#}

#p <- 110
#df <- error_prediction_data %>% 
#  pred_spaces_preprocess_error_prediction_data() %>%
#  subset(participant_test==p) %>%
#  subset(session_test==8) %>%
  #select(-c(y0,y1,y2,y3)) %>%
#  spread(model, pred_prob)

# fit <- df %>%
#   fit_ct_markov_trigram_log_exp_model()
# 
# samples <- extract(fit, c("a","b"), permuted = TRUE, inc_warmup = FALSE,
#                    include = TRUE)
# 
# exp_rt_markov <- df$Markov
# exp_rt_trigram <- df$Triplet
# rt_ct <- df$iHMM
# plot_data <- data.frame(
#   ct=rt_ct, 
#   fitted=get_pred(mean(samples$a), mean(samples$b), exp_rt_markov, exp_rt_trigram)) 
# plot_data <- plot_data %>% mutate(resid = ct-fitted)
# plot_data %>% summarise(rmse = sqrt(mean(resid^2)))
# 
# plot_data %>%
#   ggplot(aes(x=fitted, y=ct)) +
#   geom_point()

# RT_CT = const + log[ exp(RT_markov) + Alpha * exp(RT_trigram) ] + gauss error
# exp(RT) = (a * exp(RT_markov) + b * exp(RT_trigram))*exp(gauss error)
# RT_CT = const + log[ ]
#exp(RT_CT) = const * (exp(RT_markov) + Alpha * exp(RT_trigram))
#exp(RT_CT) = alpha * (exp(RT_markov)) + beta * exp(RT_trigram))

# errorf <- function(a,b, ct, markov, triplet){
#   ct - log(a*exp(markov) + b*exp(triplet))
# }
# 
# errorf(0.5, 0.5, 210, 200, 220)
# 


# df <- rt_prediction_data %>%
#   add_session_info(., dataset="elarasztas") %>%
#   filter_and_average_runs(sd_cutoff = 100) %>%
#   ungroup() %>%
#   subset(model %in% c("iHMM", "Markov", "Triplet")) %>%
#   relabel_and_order_models() %>%
#   #mutate(exp_rt_predicted = exp(rt_predicted)) %>%
#   select(-cutoff) %>%
#   spread(model, rt_predicted)
# 
# model_fits <- df %>%
#   mutate(session_test = factor(session_test, levels=c(1,2,3,4,5,6,7,8))) %>%
#   group_by(participant_test, session_test) %>%
#   do(model = lm(CT ~ Markov + Trigram - 1, data=., na.action = na.omit))
# 
# df <- error_prediction_data %>% 
#   add_session_info(., dataset="elarasztas") %>%
#   #subset(participant_test==101) %>%
#   subset(correct_response==1) %>%
#   group_by(session_test, e_train, model, participant_test, correct_response, block, trial) %>%
#   summarise(y0 = mean(y0), y1=mean(y1),y2=mean(y2),y3=mean(y3)) %>%
#   select(session_test, e_train, model, participant_test, correct_response,block,trial,y0,y1,y2,y3) %>%
#   gather(key="key", value="pred_prob", y0,y1,y2,y3)
# 
# df[df$model=="Triplet", "pred_prob"] = exp(df[df$model=="Triplet", "pred_prob"])

fig_predictability <- function(rt_prediction_data){
  model_fits <- rt_prediction_data %>%
    add_session_info("elarasztas") %>%
    #subset(session_test==8) %>%
    mutate(session_test = factor(session_test, levels=c(1,2,3,4,5,6,7,8))) %>%
    filter_and_average_runs(sd_cutoff=100, apply_filters = TRUE) %>%
    spread(model, rt_predicted) %>%
    group_by(session_test, participant_test) %>%
    do(model = lm(iHMM ~ Markov + GroundTruth, data=., na.action = na.omit))

  model_fits <- model_fits %>% mutate(rsq = summary(model)$r.squared)

  fig <- model_fits %>% 
    ggplot(aes(x=session_test, y=rsq)) +
    geom_boxplot() +
    xlab("Day") +
    ylab("R2 of model CT ~ Markov + Ideal Observer") +
    style

  return(fig)
}
