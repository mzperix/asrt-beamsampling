library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(lemon)
library(cowplot)

source(here::here('R/Fun/data_pipelines.R'))
source(here::here('R/Fun/plot/plot_styles.R'))
source(here::here('R/Fun/plot/learning_curves.R'))
source(here::here('R/Fun/plot/ihmm_markov_triplet.R'))
source(here::here('R/Fun/plot/across_sessions.R'))
source(here::here('R/Fun/plot/across_participants.R'))
source(here::here('R/Fun/plot/across_sequences.R'))
source(here::here('R/Fun/plot/error.R'))
source(here::here('R/Fun/plot/entropy.R'))
source(here::here('R/Fun/plot/fingerprint.R'))
source(here::here('R/Fun/plot/higher_order.R'))
source(here::here('R/Fun/plot/rt_cloud.R'))

fig_model_comparison <- function(data) {
  learning_curves_plot <- data %>%
    fig_learning_curves() +
    ylab("Model performance") +
    theme(
      legend.position = "right",
      legend.title = element_blank()
    )

  violin_plot <- data %>%
    relabel_and_order_models() %>%
    subset(session_test == 8) %>%
    compute_performance(., sd_cutoff = 100) %>%
    ggplot(aes(x = model, y = performance)) +
    geom_violin(aes(fill = model),
      scale = "count", trim = FALSE, bw = 0.025,
      color = pantone_2014[["dark_grey"]]
    ) +
    geom_point(color = pantone_2014[["grey"]]) +
    geom_line(aes(group = participant_test), color = pantone_2014[["grey"]]) +
    style +
    model_fills +
    guides(fill = "none") +
    ylab("Model performance") +
    theme(axis.title.x = element_blank())

  # trial_by_trial_plot <- data %>%
  #   relabel_and_order_models() %>%
  #   subset(participant_train == participant) %>%
  #   subset(e_train == "186_195") %>%
  #   fig_individual_trials() +
  #   facet_rep_wrap(model~.,
  #                  repeat.tick.labels=TRUE,
  #                  ncol=2
  #   ) +
  #   coord_fixed()+
  #   theme(
  #     strip.background = element_blank()
  #   )

  top_row <- plot_grid(learning_curves_plot, violin_plot, rel_widths = c(2, 1), nrow = 1, ncol = 2) +
    draw_label("A", 0, 1, hjust=-1, vjust=1) +
    draw_label("B", 0.66, 1, hjust=-1, vjust=1)
  # bottom_row <- plot_grid(violin_plot, trial_by_trial_plot, rel_widths=c(1,2), nrow=1, ncol=2)
  # fig <- plot_grid(top_row, trial_by_trial_plot, rel_heights=c(1,3), nrow=2, ncol=1)
  # draw_label(labels[1], 0, 1, hjust=-1, vjust=1) +
  # draw_label(labels[2], 0.66, 1, hjust=-1, vjust=1) +
  # draw_label(labels[3], 0.66, 0.2, hjust=-1, vjust=1)
  return(top_row)
}

fig2 <- function(data, participant, labels=c("A","B","C")){
  abc <- data %>%
    relabel_and_order_models() %>%
    subset(participant_train == participant) %>%
    subset(e_train == "186_195") %>%
    fig_individual_trials() +
    facet_rep_wrap(model~.,
               #scales="free",
               repeat.tick.labels=TRUE,
               ncol=2
               #nrow=length(unique(data$model))
               ) +
    #coord_capped_cart(bottom = capped_horizontal("both"),
    #                  left = capped_vertical("both"),
    #                  ratio=1) +
    coord_fixed()+
    theme(
      strip.background = element_blank()
      #strip.text.x = element_blank()
    )

  de_data <- data %>%
    relabel_and_order_models() %>%
    subset(e_train == "186_195") %>%
    subset((filters == TRUE) & (correct_response == 1)) %>%
    group_by(participant_test,e_train,session_test,model,block,trial) %>%
    summarise(rt=mean(rt), rt_predicted=mean(rt_predicted)) %>%
    group_by(participant_test, model) %>%
    summarise(performance = cor(rt, rt_predicted)^2)

  d <- de_data %>%
    ggplot(aes(x=performance, color=model, group=model, fill=model)) +
    geom_histogram() +
    facet_rep_wrap("model", ncol=1) +
    style +
    model_colors +
    model_fills +
    xlab("Model Performance") +
    ylab("Counts") +
    guides(fill=FALSE,
           group=FALSE,
           color=FALSE) +
    theme(legend.position = "bottom")

  de_means <- de_data %>%
    group_by(model) %>%
    summarise(performance=mean(performance))

  e <- de_data %>%
    ggplot(aes(x=model, y=performance, color=model, group=model, fill=model)) +
    geom_boxplot(show.legend = FALSE) +
    geom_point(size=-1) +
    geom_point(data=de_means, aes(x=model, y=performance),
      color="white", shape='square', size=1.6) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(de_data$model))) +
    style +
    model_colors +
    model_fills +
    xlab("") +
    ylab("Model Performance") +
    guides(color=FALSE,
           group=FALSE,
           fill=FALSE) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.y = element_blank()
          )

  right_col <- plot_grid(d, e, rel_heights=c(4,1), nrow=2, ncol=1)
  fig <- plot_grid(abc, right_col, rel_widths=c(2,1), nrow=1, ncol=2) +
    draw_label(labels[1], 0, 1, hjust=-1, vjust=1) +
    draw_label(labels[2], 0.66, 1, hjust=-1, vjust=1) +
    draw_label(labels[3], 0.66, 0.2, hjust=-1, vjust=1)
  return(fig)
}

fig2b <- function(data, participant, rt_cloud){
  fig <- plot_grid(rt_cloud_fig(relabel_and_order_models(rt_cloud)), fig2(data, participant, labels=c("C","D","E")),
                   rel_widths=c(1,1), nrow=1, ncol=2) +
         draw_label("A", 0, 1, hjust=0, vjust=1) +
         draw_label("B", 0, 0.33, hjust=0, vjust=1)
  return(fig)
}

fig3a <- function(data){
  data <- data %>%
    relabel_and_order_models() %>%
    compute_performance(sd_cutoff = 100)

  group_data <- data %>%
    group_by(model,session_test) %>%
    summarize(performance_mean = mean(performance),
              performance_se = sd(performance)/sqrt(n()))

  w <- 0.2
  a <- group_data %>%
    ggplot(aes(x=session_test, y=performance_mean, group=model, color=model)) +
    geom_errorbar(aes(ymin=performance_mean-2*performance_se,
                      ymax=performance_mean+2*performance_se),
                  position=position_dodge(), width=w) +
    geom_line(position=position_dodge(width=w)) +
    geom_point(position=position_dodge(width=w)) +
    ylab("Model Performance") +
    xlab("Session (Day)") +
    style +
    model_colors +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    guides(group=FALSE)
  
  return(a)
}

fig3 <- function(data, across_session_data){
  a <- fig3a(data)

  across_session_data <- across_session_data %>%
    relabel_and_order_models() %>%
    subset(session_train %in% c("Model Day 1", "Model Day 4", "Model Day 8"))

  across_session_data$session_train <- factor(
    across_session_data$e_train,
    levels=c("11_20", "36_45", "61_70", "86_95",
           "111_120", "136_145", "161_170", "186_195"),
    #levels=c("1_10", "26_35", "51_60", "76_85",
    #       "101_110", "126_135", "151_160", "176_185"),
    labels=c("Forward Prediction",
         "Backward-Forward\nPrediction","Backward-Forward\nPrediction","Backward-Forward\nPrediction",
         "Backward-Forward\nPrediction","Backward-Forward\nPrediction","Backward-Forward\nPrediction",
         "Backward Prediction"),
    ordered=TRUE)
  
  def <- across_session_data %>%
    subset(model %in% c("Markov","CT")) %>%
    fig_across_sessions(., ncol=min(4,length(unique(.$session_train)))) +
    ylim(0,0.3) +
    guides(color=FALSE)

  ggdraw() +
    draw_plot(a, x=0, y=0.4, width=1, height=0.6) +
    draw_plot(def, x=0, y=0, width=1, height=0.4) +
    draw_label("A", x=0, y=1, hjust=0, vjust=1) +
    draw_label("B", x=0, y=0.4, hjust=0, vjust=1) +
    draw_label("C", x=0.33, y=0.4, hjust=0, vjust=1) +
    draw_label("D", x=0.66, y=0.4, hjust=0, vjust=1) %>%
    return()
}

fig4 <- function(imt_data, lc_data){
  imt_data <- compute_performance(imt_data, sd_cutoff=100)

  a <- imt_data %>%
    subset(session_test == 8) %>%
    spread(., "model", "performance") %>%
    ggplot(aes(x=Triplet, y=iHMM-Markov, text=paste("participant:",participant_test))) +
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_point(color=pantone_2014[['dark_blue']], size=2) +
    coord_fixed() +
    xlab('Trigram performance') +
    ylab('CT perf. - Markov perf.') +
    style

  b_data <- imt_data %>%
    spread(., "model", "performance") %>%
    group_by(session_test, participant_test) %>%
    mutate(diff=iHMM-Markov-Triplet) %>%
    subset(session_test == 8)

  b <- b_data %>%
    ggplot(aes(x=diff)) +
    geom_histogram(bins=20, fill=pantone_2014[["violet"]]) +
    geom_vline(xintercept=mean(b_data$diff),
               color=pantone_2014[["purple"]],
               size=2) +
    coord_cartesian(expand = FALSE, xlim=c(-0.08,0.08), ylim=c(0,5)) +
    xlab("CT - Markov - Trigram") +
    ylab("Count") +
    style +
    theme(axis.text.x = element_text(angle=0))

  cd <- lc_data %>%
    subset(session_test %in% c(2,8)) %>%
    subset(model %in% c("iHMM","Triplet")) %>%
    fig_higher_order() +
    facet_wrap("session_test", ncol=2)

  ggdraw() +
    draw_plot(a, x=0, y=0.5, width=0.4, height=0.5) +
    draw_plot(b, x=0.4, y=0.5, width=0.55, height=0.5) +
    draw_plot(cd, x=0, y=0, width=1, height=0.5) #+
    #draw_label("A", x=0, y=1, hjust=-4, vjust=1) +
    #draw_label("B", x=0.4, y=1, hjust=-1, vjust=1) +
    #draw_label("C", x=0, y=0.66, hjust=-4, vjust=1) +
    #draw_label("D", x=0.5, y=0.66, hjust=-1, vjust=1) +
    #draw_label("E", x=0, y=0.33, hjust=-4, vjust=1) +
    #draw_label("F", x=0.5, y=0.33, hjust=-1, vjust=1)

}

fig4b <- function(imt_data, lc_data){ ## Change trigram to ideal observer
  imt_data <- compute_performance(imt_data, sd_cutoff=100)
  
  a <- imt_data %>%
    subset(session_test == 8) %>%
    spread(., "model", "performance") %>%
    ggplot(aes(x=GroundTruth, y=iHMM-Markov, text=paste("participant:",participant_test))) +
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_point(color=pantone_2014[['dark_blue']], size=2) +
    coord_fixed() +
    xlab('Ideal Observer performance') +
    ylab('CT perf. - Markov perf.') +
    style
  
  b_data <- imt_data %>%
    spread(., "model", "performance") %>%
    group_by(session_test, participant_test) %>%
    mutate(diff=iHMM-Markov-GroundTruth) %>%
    subset(session_test == 8)
  
  b <- b_data %>%
    ggplot(aes(x=diff)) +
    geom_histogram(bins=20, fill=pantone_2014[["violet"]]) +
    geom_vline(xintercept=mean(b_data$diff),
               color=pantone_2014[["purple"]],
               size=2) +
    coord_cartesian(expand = FALSE, xlim=c(-0.08,0.08), ylim=c(0,5)) +
    xlab("CT - Markov - Ideal Observer") +
    ylab("Count") +
    style +
    theme(axis.text.x = element_text(angle=0))
  
  cd <- lc_data %>%
    subset(session_test %in% c(2,8)) %>%
    subset(model %in% c("iHMM","GroundTruth","Triplet")) %>% # Needs triplet to identify triplet types
    fig_higher_order() #+
    #facet_wrap("session_test", ncol=2)
  
  ggdraw() +
    draw_plot(a, x=0, y=8/13, width=0.4, height=5/13) +
    draw_plot(b, x=0.4, y=8/13, width=0.55, height=5/13) +
    draw_plot(cd, x=0, y=0, width=1, height=8/13) #+
  #draw_label("A", x=0, y=1, hjust=-4, vjust=1) +
  #draw_label("B", x=0.4, y=1, hjust=-1, vjust=1) +
  #draw_label("C", x=0, y=0.66, hjust=-4, vjust=1) +
  #draw_label("D", x=0.5, y=0.66, hjust=-1, vjust=1) +
  #draw_label("E", x=0, y=0.33, hjust=-4, vjust=1) +
  #draw_label("F", x=0.5, y=0.33, hjust=-1, vjust=1)
  
}