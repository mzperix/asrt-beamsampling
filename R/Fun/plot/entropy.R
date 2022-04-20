
fig_entropy <- function(data){
  data %>%
    relabel_and_order_models() %>%
    subset(filters == TRUE) %>%
    mutate(correct_response=as.factor(correct_response)) %>%
    select(-c(commit, model, e_test, participant_train,
              permuted, first_response, block, trial, ini, rt, filters,
              V1, Y)) %>%
    gather(key="entropy_type", value="value", -c(e_train, trial_type, correct_response, participant_test)) %>%
    group_by(participant_test, e_train, trial_type, entropy_type, correct_response) %>%
    summarise(value = mean(value)) %>%
    group_by(e_train, trial_type, entropy_type, correct_response) %>%
    summarise(entropy_mean = mean(value),
              entropy_se = sd(value)/sqrt(n())) %>%
    ggplot(
      aes(x=trial_type, y=entropy_mean, 
          color=trial_type, alpha=correct_response)) +
    geom_point() +
    geom_errorbar(
      aes(ymin=entropy_mean-entropy_se, ymax=entropy_mean+entropy_se),
      size=0.7,
      width=0.6) +
    facet_grid(entropy_type~e_train, scales="free_y") +
    style +
    xlab('') + 
    ylab('') +
    scale_size_manual(values=c(1,3)) +
    scale_alpha_manual(values=c(0.2,1)) +
    coord_capped_cart(expand=FALSE,
                      left=capped_vertical("both"),
                      bottom=capped_horizontal("both"),
                      gap=0.2) + 
    theme(axis.text = element_text(margin = c(10,10,10,10)))
}

fig_entropy_ms <- function(data){
  w = 0.5
  
  d <- data %>%
    subset(filters==TRUE) %>%
    subset(correct_response==1) %>%
    subset(e_train %in% c(1,8)) %>%
    mutate(correct_response=as.factor(correct_response)) %>%
    select(-c(commit, e_test, participant_train,
              permuted, first_response, block, trial, ini, rt, filters, V1, Y)) %>%
    gather(key="entropy_type", value="value", -c(e_train, trial_type, correct_response, participant_test, model)) %>%
    group_by(participant_test, e_train, trial_type, entropy_type, correct_response, model) %>%
    summarise(value = mean(value)) %>%
    spread(key=trial_type, value=value)
  
  d <- d %>%
    mutate(diff=P-R) %>%
    group_by(e_train, entropy_type, correct_response, model) %>%
    summarise(diff_mean = mean(diff),
              diff_se = sd(diff)/sqrt(n())
    ) %>%
    ungroup() %>%
    subset(entropy_type %in% c("PRED_ENTROPIES","STATE_PRIOR_ENTROPIES","STATE_POSTERIOR_ENTROPIES","PRED_PROBS")) %>%
    mutate(entropy_type = recode(
      .$entropy_type, 
      PRED_PROBS="A. Subjective\nprobability\nof stimulus",
      PRED_ENTROPIES="B. Entropy\nof predictions",
      STATE_POSTERIOR_ENTROPIES="D. Entropy of\nstate posterior",
      STATE_PRIOR_ENTROPIES="C. Entropy of\nstate prior")) %>%
    mutate(model_day = paste(model,e_train,sep='_')) %>%
    subset(model_day %in% c("GroundTruth_1","iHMM_1","iHMM_8")) %>%
    mutate(model_day = recode(.$model_day,
                              iHMM_8="CT day 8",
                              iHMM_1="CT day 1",
                              GroundTruth_1="Expert model"))
  
  d %>%
    relabel_and_order_models() %>%
    ggplot(
      aes(x=model_day, y=diff_mean, fill=model
          #alpha=correct_response,
          #color=model
      )) +
    geom_errorbar(aes(ymin=-abs(diff_mean)-2*diff_se,ymax=abs(diff_mean)+2*diff_se),
                  color = "white", alpha=0) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=diff_mean-2*diff_se, ymax=diff_mean+2*diff_se),
                  width=0.5, color=pantone_2014[['dark_grey']]) +
    facet_wrap(.~entropy_type, scales = "free", ncol=4) +
    model_colors +
    model_fills +
    xlab("") +
    ylab("Pattern-Random trials\n(Group mean)") +
    style +
    guides(fill=FALSE) + 
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle=90),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank())
}
