library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)
library(Rtsne)

fig_across_participants <- function(data){
  data <- data %>%
    mutate(permuted=ifelse(permuted==TRUE,'sequence matched','not permuted')) %>%
    mutate(same=ifelse(participant_train==participant_test,'Within participant','Across participant')) %>%
    mutate(performance=correlation^2*sign(correlation))
  plot_across <- data %>%
    relabel_and_order_models() %>%
    ggplot(aes(x=performance, color=model, linetype=same)) + 
    #stat_bin(aes(y=..density..),position="dodge2", bins=20, color="white", size=0.3) +
    geom_line(stat="density", size=1.2) +
    #stat_density(adjust=0.5) +
    facet_grid(permuted~., scales="fixed") +
    style +
    model_fills +
    model_colors + 
    xlab('Model Performance') +
    ylab('Density') +
    theme(legend.position = "top",
          axis.text.y = element_blank()) +
    guides(fill=FALSE,
           linetype=guide_legend(title=""),
           color=guide_legend(title="Model"))
  return(plot_across)
}

fig_across_participants_bars <- function(data){
  w <- 0.95
  data <- data %>%
    relabel_and_order_models() %>%
    mutate(permuted=ifelse(permuted==TRUE,'sequence\nmatched','not\npermuted')) %>%
    mutate(same=ifelse(participant_train==participant_test,'Within participant','Across participant')) %>%
    mutate(performance=correlation^2*sign(correlation)) %>%
    group_by(model, same, permuted) %>%
    summarise(group_performance = mean(performance),
              se = sd(performance)/sqrt(n())) %>%
    mutate(same_permuted = paste(same, permuted, sep='\n')) %>%
    #subset((same=='Across participant') | (permuted=='sequence matched')) %>%
    subset(same=='Across participant')
  
  data[data$same_permuted == "Within participant\nsequence matched","same_permuted"] <- "Within participant"
  
  bar_plot <- data %>%
    ggplot(aes(x=permuted, y=group_performance, fill=model, color=model)) +
    geom_bar(stat="identity", position=position_dodge(width=w)) +
    geom_errorbar(
      aes(ymin=group_performance-se, ymax=group_performance+se), 
      position=position_dodge(width=w),
      width=0.5,
      color=pantone_2014[["dark_grey"]]) +
    style +
    model_fills +
    model_colors +
    xlab("") +
    ylab(expression(Mean~R^{2})) +
    coord_cartesian(expand = FALSE, xlim=c(0.3,2.7), ylim=c(0,0.13)) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle=0, vjust=0.5)) +
    guides(fill=guide_legend(title="Model"),
           color=FALSE)
  
  return(bar_plot)
}

fig_across_participants_all_days <- function(data){
  w <- 0.3
  data <- data %>%
    relabel_and_order_models() %>%
    mutate(permuted=ifelse(permuted==TRUE,'sequence matched','not permuted')) %>%
    subset(participant_test!=participant_train) %>%
    #mutate(same=ifelse(participant_train==participant_test,'Within participant','Across participant')) %>%
    mutate(performance=correlation^2*sign(correlation))
  
  p <- data %>%
    group_by(e_train, permuted, model) %>%
    summarise(performance_mean = mean(performance),
              performance_se = sd(performance)/sqrt(n()))
  
  p <- p %>%
    ggplot(aes(x=e_train, y=performance_mean, 
               color=model, group=model)) +
    geom_errorbar(aes(ymin=performance_mean-2*performance_se,
                      ymax=performance_mean+2*performance_se),
                  width=w,
                  position = position_dodge()) +
    geom_point(position = position_dodge(width=w)) +
    geom_line(position = position_dodge(width=w)) + 
    facet_grid(rows="permuted") +
    style +
    model_colors
  return(p)
}

cycle_sequence <- function(seq){
  while(seq[1] != 0){
    seq <- c(tail(seq,3),seq[1])
  }
  return(seq)
}

fig_across_participants_tsne <- function(data, sequence_data, model_filter, per, perplexity, theta, initial_dims){
  for (i in 1:nrow(sequence_data)){
    sequence_data[i,'seq'] <- paste0(cycle_sequence(sequence_data$sequence[[i]]), collapse='')
  }
  
  sequence_data <- sequence_data %>%
    mutate(participant_train = substr(participant_train,12,14),
           seq = as.factor(seq)) %>%
    select(participant_train, seq)
  
  df <- data %>%
    subset(permuted==per) %>%
    #subset(trial_type=="R") %>%
    #subset(participant_test == 101) %>%
    subset(model==model_filter) %>%
    group_by(trial,block,participant_test,participant_train,e_train,e_test,model) %>%
    summarise(rt_predicted=1/mean(rt_predicted)) %>%
    #subset(trial>5) %>%
    #subset(block<6) %>%
    # normalize
    group_by(participant_test, participant_train, e_train, e_test, model) %>%
    mutate(normalised_pred_rt = (rt_predicted - mean(rt_predicted)) / sd(rt_predicted)) %>%
    mutate(N = paste0('trial_',trial+(block-1)*85)) %>%
    select(participant_train, participant_test, e_train, e_test, model, N, normalised_pred_rt) %>%
    spread(key=N, value=normalised_pred_rt) %>%
    ungroup()
  
  df$e_train <- as.factor(plyr::mapvalues(
    df$e_train, 
    from=c("11_20", "36_45", "61_70", "86_95", 
           "111_120", "136_145", "161_170", "186_195"),
    to=c(1,2,3,4,5,6,7,8)))
  
  plot_data <- df %>%
    select(-c(participant_train,participant_test,e_train,e_test,model))
  plot_data <- Rtsne(plot_data, perplexity=perplexity, theta=theta, 
                     initial_dims=initial_dims, normalize = FALSE, max_iter=3000)[['Y']] %>%
    data.frame() %>%
    cbind(select(df,participant_train,participant_test,e_train,e_test,model)) %>%
    mutate(participant_train = factor(participant_train))
  
  plot_data <- merge(plot_data, sequence_data, by="participant_train")
  plot_data <- plot_data[order(plot_data$e_train),]
  #return(plot_data)
  
  p <- plot_data %>%
    ggplot(aes(x=X1,y=X2, color=participant_train, size=e_train, shape=seq)) + 
    #geom_path(aes(group=participant_train), size=1) +
    geom_point() +
    scale_shape_manual(values=c(0,1,2,5,6,9)) +
    coord_fixed() +
    style +
    theme(panel.grid.major = element_blank())
  
  return(p)
}