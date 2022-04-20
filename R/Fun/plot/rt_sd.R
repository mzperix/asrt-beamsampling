library(ggplot2)
source(here::here('R/Fun/data_pipelines.R'))

#####################################
###   RT Descriptive stats (sd)   ###
#####################################
rt_sd <- function(rt_prediction_data){
  data <- rt_prediction_data %>%
    add_session_info(., dataset="elarasztas") %>%
    filter_and_average_runs(sd_cutoff=100) %>%
    group_by(participant_test,session_test) %>%
    summarise(mean_rt = mean(rt),
              sd_rt = sd(rt),
              min_rt = min(rt),
              max_rt = max(rt))
  
  data %>%
    ungroup() %>%
    mutate(session_test=as.character(session_test),
           participant_test=as.factor(participant_test)) %>%
    ggplot(aes(x=participant_test, y=sd_rt, color=session_test)) +
    geom_point(shape=1, size=1.5, stroke=1.0) +
    scale_color_brewer(palette="Spectral") +
    theme_bw() +
    xlab('Participant') +
    ylab('Standard Deviation of Response Times') +
    guides(color=guide_legend(title='Session')) +
    theme(axis.text.x=element_text(angle=90))
}

#ggsave('Figures/rt_sd.pdf',width=5, height=4.3)