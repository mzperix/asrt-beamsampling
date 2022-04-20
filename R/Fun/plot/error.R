source(here::here('R/Fun/data_pipelines.R'))

fig_error <- function(data){
  data <- data %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    group_by(participant_test, e_train, model, block, trial, Y, correct_response, first_response) %>%
    summarise(y0=mean(y0),y1=mean(y1),y2=mean(y2),y3=mean(y3)) %>%
    ungroup()

  ranks <- t(apply(select(data,c(y0,y1,y2,y3)),1, function(x) rank(-x, ties.method = "random")))
  colnames(ranks) <- paste(colnames(ranks),"_rank",sep="")
  pred_stim <- data.frame(data,ranks) %>%
    select(-(y0:y3)) %>%
    gather(stimulus, rank, -(participant_test:first_response)) %>%
    subset(paste('y',Y,'_rank',sep='')==stimulus) %>%
    group_by(correct_response, model, rank)
  pred_stim[pred_stim$correct_response==0, 'correct_response'] <- 'Incorrect Response'
  pred_stim[pred_stim$correct_response==1, 'correct_response'] <- 'Correct Response'
  
  error_plot_stim <- pred_stim %>%
    #subset(correct_response == 0) %>%
    ggplot(aes(x=rank, fill = model, group=model)) +
    geom_bar(position="dodge", color='white', size=0.5) +
    facet_wrap("correct_response", scales="free") +
    style +
    xlab('Prediction Rank of Stimulus') +
    ylab('Count') +
    model_fills
  
  pred_resp <- data
  pred_resp$first_key <- plyr::mapvalues(pred_resp$first_response,
                                         from=c("z","c","b","m"),
                                         to=c(0,1,2,3))
  ranks <- t(apply(select(data,c(y0,y1,y2,y3)),1, function(x) rank(-x, ties.method='random')))
  colnames(ranks) <- paste(colnames(ranks),"_rank",sep="")
  pred_resp <- data.frame(pred_resp,ranks) %>%
    select(-(y0:y3)) %>%
    gather(stimulus, rank, -(participant_test:first_response), -first_key) %>%
    subset(paste('y',first_key,'_rank',sep='')==stimulus) %>%
    group_by(e_train, correct_response, model, rank)
  pred_resp[pred_resp$correct_response==0, 'correct_response'] <- 'Incorrect Response'
  pred_resp[pred_resp$correct_response==1, 'correct_response'] <- 'Correct Response'
  
  error_plot_resp <- pred_resp %>%
    subset(correct_response=='Incorrect Response') %>%
    ggplot(aes(x=rank, fill = model, group=model)) +
    geom_bar(position="dodge", color='white', size=0.5) +
    facet_wrap("correct_response", scales="free") +
    style +
    xlab('Prediction Rank of Response') +
    ylab('Count') +
    model_fills +
    guides(fill=FALSE)
  
  w <- 0.68
  correction <- 0.04
  fig <- ggdraw() +
    draw_plot(error_plot_stim, 0, 0, w-correction, 1) +
    draw_plot(error_plot_resp, w-correction, 0, 1-w-correction, 1) +
    draw_label('A',0,1, vjust=1, hjust=-0.5) +
    draw_label('B',w-correction,1, vjust=1, hjust=-0.5)
  
  return(fig)
}

fig_error_b <- function(data){
  data <- data %>%
    add_session_info(., dataset="elarasztas") %>%
    subset(session_test==8)
  
  ranks <- compute_prediction_ranks(data)
  stimulus_ranks <- ranks[["stimulus_ranks"]]
  response_ranks <- ranks[["response_ranks"]]
  p_stimulus_ranks <- stimulus_ranks %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    ggplot(aes(x=model, y=proportion,
               group=correct_response, fill=model, color=model,
               alpha=correct_response)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept=0.25) +
    facet_grid(.~response_type) +
    coord_cartesian(expand = FALSE, xlim=c(0,length(unique(stimulus_ranks$model))+1), ylim=c(0,0.7)) +
    scale_x_discrete(breaks=NULL) +
    ggtitle("Stimulus") +
    ylab("Proportion of top rank") +
    xlab("") +
    style +
    model_colors +
    model_fills +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(alpha=FALSE)
  
  p_response_ranks <- response_ranks %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    ggplot(aes(x=model, y=proportion,
               fill=model, color=model)) +
    geom_bar(stat="identity", position=position_dodge(), alpha=0.1) +
    geom_hline(yintercept=0.25) +
    facet_grid(.~response_type) +
    coord_cartesian(expand = FALSE, xlim=c(0,length(unique(response_ranks$model))+1), ylim=c(0,0.7)) +
    scale_x_discrete(breaks=NULL) +
    ggtitle("Response") +
    ylab("") +
    xlab("") +
    style +
    model_colors +
    model_fills +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(alpha=FALSE, fill=FALSE, color=FALSE)
  
  plot_grid(p_stimulus_ranks, p_response_ranks, ncol=2, rel_widths = c(2.6,1)) +
    draw_label('A', 0,1, hjust=0, vjust=1.5) +
    draw_label('B', 2.6/3.6,1, hjust=0, vjust=1.5)
}

fig_error_roc <- function(roc_data_all, auc, participants=NULL, ncol=5){
  if (is.null(participants)){
    participants <- unique(roc_data_all$participant_test)
  }
  
  roc_data_all <- roc_data_all %>%
    ungroup() %>%
    relabel_and_order_models()
  
  auc_wide <- auc %>%
    spread(model, area) %>%
    mutate(dauc = iHMM-GroundTruth)
  
  fig1 <- roc_data_all %>%
    subset(model != "Trigram") %>%
    subset(participant_test %in% participants) %>%
    ggplot(
      aes(x=FPR, y=TPR, color=model)) +
    geom_line(size=0.7) +
    # geom_point(data=subset(unique(roc_data_all), (model=="Trigram") &(participant_test %in% participants)),
    #            size=2.0)+
    style + 
    model_colors +
    coord_fixed() + 
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_text(angle=90)) +
    geom_text(data=subset(auc_wide, participant_test %in% participants),
              aes(label=paste0("AUC diff = ", format(dauc, digits=3))),
              color="black",
              x=Inf, y=-Inf,
              hjust=1, vjust=-1) +
    facet_wrap("participant_test", ncol=ncol)
  
  return(fig1)
}

fig_error_roc_all <- function(data, ncol=6){
  roc_auc <- compute_roc_auc(data)
  roc_data_all <- roc_auc[[1]]
  auc <- roc_auc[[2]]
  fig1 <- fig_error_roc(roc_data_all, auc, ncol=ncol)
  return(fig1)
}

fig_error_c <- function(data, participants, figs_as_list=FALSE){
  roc_auc <- compute_roc_auc(data)
  roc_data_all <- roc_auc[[1]]
  auc <- roc_auc[[2]]
  
  fig1 <- fig_error_roc(roc_data_all, auc, participants, ncol=1)
  
  auc <- auc %>%
    ungroup() %>%
    relabel_and_order_models() %>%
    mutate(model=factor(model, levels=c("CT","Ideal Observer","Trigram","Markov")))
  
  fig2 <- auc %>%
    ggplot(aes(x=model, y=area, group=participant_test)) +
    geom_bar(data=summarise(group_by(auc, model), mean_area=mean(area)), stat="identity",
             aes(x=model, y=mean_area, fill=model, group=model)) +
    geom_point(color=pantone_2014[['grey']]) +
    geom_line(color=pantone_2014[['grey']]) +
    coord_cartesian(expand = FALSE, xlim=c(0,length(unique(auc$model))+1), ylim=c(0.5,0.8)) +
    model_colors +
    model_fills +
    xlab("") +
    ylab("Area Under ROC Curve") +
    style +
    scale_x_discrete(breaks=NULL) +
    guides(fill=FALSE) +
    theme(legend.title = element_blank())
  
  if (figs_as_list){
    return(list(fig1, fig2))
  }
  return(plot_grid(fig1, fig2, rel_widths = c(1,1)))
  
}

fig_error_d <- function(data, participants){
  c_figs <- fig_error_c(data, participants, figs_as_list = TRUE)
  plot_grid(
    fig_error_b(data),
    plot_grid(c_figs[[1]] + facet_wrap("participant_test", ncol=2) + guides(color=FALSE),
              c_figs[[2]], rel_widths=c(2.2,1), nrow=1),
    rel_heights = c(1,1.2),
    nrow = 2) +
    draw_label('C', 0, 1.2/2.2, hjust=0, vjust=1.0) +
    draw_label('D', 2.2/3.2, 1.2/2.2, hjust=0, vjust=1.0)
}

fig_error_e <- function(data, participants){
  c_figs <- fig_error_c(data, participants, figs_as_list = TRUE)
  plot_grid(
      fig_error_b(data),
      plot_grid(
        c_figs[[1]] + facet_wrap("participant_test", ncol=2) + guides(color=FALSE),
        plot_grid(c_figs[[2]],ggplot(), rel_widths = c(1,2), ncol=2), 
        rel_heights=c(1.1,1.0), 
        nrow=2),
      rel_heights = c(1,2.2),
      nrow = 2)
}
