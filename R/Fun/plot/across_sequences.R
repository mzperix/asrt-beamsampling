source(here::here('R/Fun/data_pipelines.R'))

across_sequence_linetypes <- scale_linetype_manual(values=c(
  "211_220"="dashed", "236_245"="solid", 
  "176_185"="solid", "201_225"="dashed",
  "186_195"="solid",
  "231_235"="dashed", "236_240"="solid", "241_245"="dashed",
  "Day 8"="solid", "Day 9"="dashed"))

fig_across_sequences <- function(data){
  w_dodge <- 0.3
  w_error <- 0.4
  
  plot_data <- across_sequences_performances(data)
  
  across_sequences_plot <- plot_data %>%
    subset(model=="iHMM") %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    group_by(epoch,e_train,model) %>%
    summarise(mean_performance = mean(performance), se_performance = sd(performance)/sqrt(n())) %>%
    mutate(model_e_train=paste(model,e_train)) %>%
    ggplot(aes(x=epoch, y=mean_performance, color=model, group=model_e_train, linetype=e_train, shape=e_train)) +
    geom_line(position = position_dodge(width=w_dodge),) +
    geom_errorbar(aes(ymin=mean_performance-se_performance, ymax=mean_performance+se_performance,
                      width=w_error),
                  position = position_dodge(width=w_dodge),) +
    geom_vline(aes(xintercept = 1.5), color=pantone_2014[['grey']], size=1.7) +
    geom_point(position = position_dodge(width=w_dodge), size=4) +
    geom_vline(aes(xintercept = 2.5), color=pantone_2014[['grey']], size=1.7) +
    style +
    model_colors +
    across_sequence_linetypes + 
    scale_shape_manual(values=list("186_195"=16,"211_220"=1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #ggtitle('Predicting different ASRT sequence \n (session 8 to 9) GROUP') +
    xlab('Blocks') +
    ylab('Model Performance') +
    guides(color=guide_legend(title="Model"))
  return(across_sequences_plot)
}

fig_across_sequences_individual <- function(data){
  plot_data <- across_sequences_performances(data)
  plot_data %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    subset(model=="CT") %>%
    mutate(model_e_train=paste(model,e_train),
           e_train_day=cut(
             as.integer(gsub('_[0-9]+','',e_train)),
             c(0,25,50,75,100,125,150,175,200,225),
             labels=paste('Day',c(1,2,3,4,5,6,7,8,9)))) %>% 
    ggplot(aes(x=epoch, y=performance, color=model, group=model_e_train, linetype=e_train_day, shape=e_train_day)) +
    geom_line() +
    geom_point() +
    facet_wrap("participant_test", ncol=5) +
    geom_vline(aes(xintercept = 1.5), color=pantone_2014[['grey']], size=1.7) +
    geom_vline(aes(xintercept = 2.5), color=pantone_2014[['grey']], size=1.7) +
    geom_vline(aes(xintercept = 6.5), color=pantone_2014[['grey']], size=1.7) +
    style +
    model_colors + 
    across_sequence_linetypes + 
    scale_shape_manual(values=list("Day 8"=16,"Day 9"=1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    xlab('Blocks') +
    ylab('Model Performance (var. expl.)') +
    guides(color=FALSE, #guide_legend(title="Model"),
           linetype=guide_legend(title="Model", ),
           shape=FALSE)
}

fig_across_sequences_ms <- function(data){
  indiv_data <- data %>%
    subset(model %in% c('iHMM','Markov')) %>%
    mutate(e_train = factor(.$e_train, 
      levels=c("186_195","211_220", "236_245", "236_240"),
      labels=c("Day 8",
               "Day 9",
               "Day 10",
               "Day 10"),
      ordered=TRUE)) %>%
    mutate(e_test = cut(
      block, 
      breaks=c(184,200,210,220,225,235,245,250), 
      labels = c('Day 8', 'Day 9', 'Day 9 Train','Day 9 End','Day 10','Day 10 Train','Day 10 Final'))) %>%
    subset(e_test %in% c('Day 8')) %>%
    subset(correct_response==1) %>%
    subset(filters == TRUE) %>%
    group_by(participant_test, e_train, block, model, trial, e_test) %>%
    summarise(rt_predicted=mean(rt_predicted), rt=mean(rt)) %>%
    group_by(participant_test, e_train, e_test, model) %>%
    mutate(cutoff = mean(rt) + 3 * sd(rt)) %>%
    subset(rt < cutoff) %>%
    summarise(performance = cor(rt,rt_predicted)) %>%
    mutate(performance = ifelse(model=="Triplet", -performance, performance)) %>%
    mutate(performance = sign(performance)*performance^2)
  
  d <- indiv_data %>%
    subset(model == 'iHMM') %>%
    spread(e_train, performance)
  
  #p_value <- t.test(x=d$`Day 9`, y=d$`Day 10`, paired = TRUE, alternative=c("less"))
  
  indiv_data <- indiv_data %>%
    mutate(model_e_train = factor(
      paste(model, e_train), 
      levels=c('iHMM Day 9','iHMM Day 10','Markov Day 9','Markov Day 10')))
  
  group_data <- indiv_data %>%
    group_by(e_train, e_test, model, model_e_train) %>%
    summarise(group_mean = mean(performance),
              group_se = sd(performance)/sqrt(n()))
  
  fig <- group_data %>%
    ggplot(aes(x=model_e_train, y=group_mean, 
               fill=model_e_train, color=model_e_train)) +
    geom_bar(stat="identity") +
    # geom_errorbar(aes(ymin=group_mean-2*group_se, ymax=group_mean+2*group_se),
    #               color=pantone_2014[['dark_grey']], width=0.5) +
    geom_point(data=indiv_data, aes(x=model_e_train, y=performance),
               color=pantone_2014[['light_grey']], show.legend = FALSE) +
    geom_line(data=subset(indiv_data,model=='iHMM'), 
              aes(x=model_e_train, y=performance, group=participant_test), 
              color=pantone_2014[['light_grey']]) +
    geom_line(data=subset(indiv_data,model=='Markov'), 
              aes(x=model_e_train, y=performance, group=participant_test), 
              color=pantone_2014[['light_grey']]) +
    facet_grid(~e_test) +
    coord_cartesian(expand = FALSE, xlim=c(0,5), ylim=c(0,0.45)) +
    ylab('Model Performance\n(Group mean)') +
    style +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_discrete(breaks=NULL) +
    model_fills +
    model_colors
  
  return(fig)
}
