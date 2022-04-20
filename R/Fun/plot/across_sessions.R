fig_across_sessions <- function(data, ncol=4){
  w_dodge <- 0.4
  p <- data %>%
    ungroup() %>%
    #relabel_and_order_models() %>%
    mutate(performance = correlation^2) %>%
    group_by(session_train, session_test, model) %>%
    summarize(group_performance = mean(performance),
              se = sd(performance)/sqrt(n())) %>%
    ggplot(aes(x=session_test, y = group_performance, group = model, color = model)) +
    style +
    geom_point(position = position_dodge(width=w_dodge)) +
    geom_line(position = position_dodge(width=w_dodge)) +
    geom_errorbar(aes(ymin = group_performance-se, ymax = group_performance+se),
                  width = w_dodge, position = "dodge") +
    xlab('Session (Day)') +
    ylab('Model performance') +
    #ggtitle('Learning curves') +
    model_colors +
    guides(colour = guide_legend(title="Model", label.position="right")) +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.5,"cm"),
          legend.spacing.x = unit(0.5,"cm"),
          legend.text.align = 0) +
    facet_rep_wrap('session_train', ncol = ncol, repeat.tick.labels=TRUE)
  return(p)
}

fig_across_sessions_individual <- function(data){
  w_dodge <- 0.4
  p <- data %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    mutate(performance = correlation^2) %>%
    ggplot(aes(x=session_test, y = performance, group = model, color = model)) +
    style +
    geom_point() +
    geom_line() +
    xlab('Session (Day)') +
    ylab('Model performance') +
    model_colors +
    guides(colour = guide_legend(title="Model", label.position="right")) +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.5,"cm"),
          legend.spacing.x = unit(0.5,"cm"),
          legend.text.align = 0) +
    facet_grid(participant_train ~ session_train)
  return(p)
}

fig_across_sessions_individual_sliced <- function(data, pages=3){
  participants <- unique(data$participant_test)
  n <- ceiling(length(participants)/pages)
  plots <- list()
  for (p in 1:pages){
    plots[[length(plots)+1]] <- data %>%
      subset(participant_test %in% participants[(n*(p-1)+1):(n*p)]) %>%
      fig_across_sessions_individual()
  }
  return(plots)
}
