library(ggplot2)
library(magrittr)
library(dplyr)
library(cowplot)
library(lemon)

source(here::here('R/Fun/data_pipelines.R'))
source(here::here('R/Fun/plot/plot_styles.R'))

fig_learning_curves <- function(data){
  sl_plot_data <- data %>% 
    relabel_and_order_models() %>%
    compute_performance(., sd_cutoff=100)

  mean_data <- sl_plot_data %>%
    group_by(e_train, session_test, model) %>%
    summarize(group_performance = mean(performance),
              se = sd(performance)/sqrt(n()))

  p <- sl_plot_data %>%
    ggplot(aes(x=session_test, y = performance, fill = model)) +
    geom_line(data = mean_data,
              aes(x=session_test, y=group_performance, group=model, color=model),
              position = position_dodge(width=0.75),
              size=1.3) +
    geom_boxplot(data=sl_plot_data, aes(x=session_test, y = performance, color=model),
                 alpha=0,
                 fatten=NULL) +
    geom_point(data = mean_data,
               aes(x=session_test, y=group_performance, group=model, color=model),
               position = position_dodge(width=0.75),
               size=2.0) +
    xlab('Session (Day)') +
    ylab('Model performance (corr. squared)') +
    #ggtitle('Statistical learning') +
    model_colors +
    model_fills +
    #coord_cartesian(ylim=c(0,0.4), expand=FALSE) +
    style
  return(p)
}

fig_learning_curves_individual <- function(data){
  w <- 0.2
  sl_plot_data <- compute_performance(data, sd_cutoff=3)

  sl_plot_data %>%
    relabel_and_order_models() %>%
    ggplot(aes(x=session_test, y=performance, group=model, color=model)) +
    geom_line() +
    geom_point() +
    xlab("Session") +
    ylab("Model performance (var. expl.)") +
    guides(color=guide_legend(title="Model")) +
    facet_wrap("participant_test", ncol=5) +
    style +
    theme(legend.position = "bottom") +
    model_colors
}

fig_learning_curves_ct_vs_markov <- function(data){
  w <- 0.2
  plot_data <- data %>%
    compute_performance(sd_cutoff=3) %>%
    subset(model %in% c("iHMM", "Markov"))
  
  group1 <- c(101, 110, 131, 128, 129, 118, 109)
  
  grouped_data <- plot_data %>%
    spread(key=model, value=performance) %>%
    group_by(participant_test) %>%
    mutate(difference=iHMM-Markov,
           var_expl=1-(1-iHMM)/(1-Markov),
           difference=difference-mean(difference)) %>%
  #  mutate(difference=difference-mean(difference)) %>%
    mutate(group=ifelse(participant_test %in% group1, 1, 2))
  
  plot <- grouped_data %>%
    ggplot(aes(x=session_test, group=participant_test, color=group)) +
    # geom_point(aes(y=var_expl)) +
    # geom_line(aes(y=var_expl)) +
    geom_point(aes(y=difference)) +
    geom_line(aes(y=difference)) +
    # geom_point(data=subset(plot_data, model=="iHMM"), 
    #            aes(x=session_test, y=performance),
    #            color=model_colors_list[["CT"]],
    #            alpha=0.3) +
    # geom_line(data=subset(plot_data, model=="iHMM"),
    #           aes(x=session_test, y=performance),
    #           color=model_colors_list[["CT"]],
    #           alpha=0.3) +
    facet_wrap("participant_test", ncol=5) +
    #facet_wrap("group", ncol=2) +
    guides(color=FALSE) +
    xlab("Session (Day)") +
    ylab("CT - Markov performance") +
    style
  
  group1 <- c(101, 107, 109, 111, 118, 120, 122, 124, 126, 129)
  
  plot <- plot_data %>%
    subset(model == "Markov") %>%
    #subset(session_test %in% c(3,4,5,6,7,8)) %>%
    group_by(participant_test) %>%
    mutate(performance = performance-mean(performance)) %>%
    mutate(group=ifelse(participant_test %in% group1, 1, 2)) %>%
    ggplot(aes(x=session_test, group=participant_test, color=group)) +
    geom_point(aes(y=performance)) +
    geom_line(aes(y=performance)) +
    facet_wrap("participant_test", ncol=5) +
    #facet_wrap("group", ncol=2) +
    guides(color=FALSE) +
    xlab("Session (Day)") +
    ylab("CT - Markov performance") +
    style
  plot
  
  return(plot)
}

fig_individual_trials <- function(data){
  plot_data <- data %>%
    subset((filters == TRUE) & (correct_response == 1)) %>%
    group_by(participant_test,e_train,session_test,model,block,trial) %>%
    summarise(rt=mean(rt), rt_predicted=mean(rt_predicted))

  correlations <- data %>%
    compute_performance(., sd_cutoff=100)

  plot_data %>%
    group_by(participant_test, e_train, session_test, model) %>%
    mutate(rt_mean = mean(rt),
           rt_sd = sd(rt)) %>%
    subset(rt<rt_mean+2*rt_sd) %>%
    ggplot(aes(x=rt_predicted, y=rt, color=model)) +
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_point(size=0.7) +
    facet_rep_wrap(participant_test ~ model,
                   scales="free",
                   #repeat.tick.labels="bottom",
                   ncol=length(unique(plot_data$model))) +
    style +
    model_colors +
    guides(color=FALSE) +
    geom_text(data=correlations,
               aes(label=paste0("R^2 == " , format(performance, digits=3))),
               parse=TRUE,
               angle=90,
               x=Inf, y=-Inf,
               hjust=-0.2, vjust=0) +
    ylab('Response time') +
    xlab('Predicted Response Time')
}

fig_binned_predicted_and_actual_means <- function(data, participant){
  plot_data <- data %>%
    subset((filters == TRUE) & (correct_response == 1)) %>%
    subset(participant_test==participant) %>%
    subset(session_test!=1) %>%
    subset(session_test==8) %>%
    group_by(participant_test,e_train,session_test,model,block,trial) %>%
    summarise(rt=mean(rt), rt_predicted=mean(rt_predicted)) %>%
    group_by(participant_test, e_train, session_test, model) %>%
    mutate(rt_mean = mean(rt),
           rt_sd = sd(rt)) %>%
    subset(rt<rt_mean+3*rt_sd) %>%
    mutate(bin = cut(rt_predicted, breaks = seq(100,1600,5))) %>%
    group_by(participant_test, e_train, session_test, model, bin) %>%
    summarise(rt_mean = mean(rt),
              rt_se = sd(rt)/sqrt(n()),
              bin_mean = mean(rt_predicted))
  
  correlations <- plot_data %>%
    group_by(participant_test,model) %>%
    summarise(correlation = cor(rt_mean, bin_mean))
  
  plot_data %>%
    ggplot(aes(x=bin_mean, y=rt_mean, color=model)) +
    #geom_boxplot(aes(y=rt)) + 
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_point(size=0.7) +
    geom_errorbar(aes(ymax=rt_mean+rt_se, ymin=rt_mean-rt_se)) +
    facet_rep_wrap(participant_test ~ model,
                   scales="free",
                   #repeat.tick.labels="bottom",
                   ncol=length(unique(plot_data$model))) +
    style +
    model_colors +
    model_fills +
    guides(color=FALSE) +
    geom_text(data=correlations,
              aes(label=paste0("R^2 == " , format(correlation^2, digits=3))),
              parse=TRUE,
              angle=90,
              x=Inf, y=-Inf,
              hjust=-0.2, vjust=0) +
    ylab('Response time') +
    xlab('Predicted Response Time') +
    coord_fixed()
}

# library(readr)
# rt_prediction_data <- read_csv(here::here('Python/Output/832a60f/866066d_LEARNING_CURVES_elarasztas_pre_submission.csv'))
# data <- add_session_info(rt_prediction_data, "elarasztas")# %>%
#  filter_and_average_runs(100) %>%
#  ungroup() %>%
#  relabel_and_order_models() 
#fig_binned_predicted_and_actual_means(data,102)
