## Data pipelines for connecting the input data and plots / stats
##          --> plot
## input --X
##          --> stats

add_session_info <- function(data, dataset){
  if (dataset == "elarasztas"){
    data <- data %>%
      mutate(session_test = cut(
        block, c(0,25,50,75,100,125,150,175,200,225,250),
        labels = c(1,2,3,4,5,6,7,8,9,10)))
  }
  if (dataset == "ASRT_180"){
    data <- data %>%
      mutate(session_test = cut(
        block, c(0,5,10,15,20,25,30,35,40,45),
        labels = c(1,2,3,4,5,6,7,8,9)))
  }
  if (dataset == "TWM"){
    data <- data %>%
      mutate(session_test = cut(
        block, c(0,5,10,15,20,25),
        labels = c(1,2,3,4,5)))
  }
  if (dataset == "implicit_erp"){
    data <- data %>%
      mutate(session_test = cut(
        block, c(0,5,10,15,20,25,30),
        labels = c(1,2,3,4,5,6)))
  }
  if (dataset == "artificial_asrt_ASRT"){
    data <- data %>%
      mutate(session_test = cut(
        block, c(0,10),
        labels = c(1)))
  }
  return(data)
}

filter_and_average_runs <- function(rt_prediction_data, sd_cutoff, apply_filters=TRUE){
  if (apply_filters){
    rt_prediction_data <- rt_prediction_data %>%
      subset((filters == TRUE) & (correct_response == 1))
  }
  rt_prediction_data %>%
    group_by(participant_test,e_train,session_test,model,block,trial,trial_type,correct_response,filters) %>%
    summarise(rt=mean(rt), rt_predicted=mean(rt_predicted)) %>%
    group_by(participant_test,model) %>%
    mutate(cutoff = mean(rt)+sd_cutoff*sd(rt)) %>%
    subset(rt<cutoff) %>%
    return(.)
}

compute_performance <- function(rt_prediction_data, sd_cutoff=100){
  rt_prediction_data %>%
    filter_and_average_runs(., sd_cutoff, apply_filters=TRUE) %>%
    group_by(participant_test,e_train,session_test,model) %>%
    summarize(performance = cor(rt,rt_predicted)) %>%
    mutate(performance=performance^2*sign(performance))
}

compute_ihmm_markov_triplet_diff <- function(rt_prediction_data){
  compute_performance(rt_prediction_data, sd_cutoff=100) %>%
    spread(., "model","performance") %>%
    mutate(diff=iHMM-Markov-Triplet)
}

compute_ihmm_markov_io_diff <- function(rt_prediction_data){
  compute_performance(rt_prediction_data, sd_cutoff=100) %>%
    spread(., "model","performance") %>%
    mutate(diff=iHMM-Markov-GroundTruth)
}

get_triplet_types <- function(data){
  d <- data.table::setorder(data, participant_test,e_train,block,trial)

  # All pattern trials are high triplets.
  # All random trials that have the same prediction as pattern trials are also high
  # That is, if the prediction probability of the random trial is identical
  # to previous trial's prediction, we have high triplet.
  
  # Only look at trials after 7 after we filter
  d[(d$trial<=7),'filters'] <- FALSE

  d <- d %>%
    group_by(participant_test, e_train, session_test) %>%
    mutate(TT = ifelse(round(Triplet, digits=4)==round(min(Triplet, na.rm=TRUE),digits=4),'H','L'))
  return(d)
}

compute_higher_order_stats <- function(data){
  d <- data %>%
    add_session_info(dataset="elarasztas") %>%
    filter_and_average_runs(sd_cutoff=100, apply_filters=FALSE) %>%
    spread(.,"model","rt_predicted")
  
  d <- get_triplet_types(d) %>%
    subset(filters==TRUE) %>%
    subset(correct_response==1)
  
  d <- d %>%
    gather("model", "rt_predicted", 
           intersect(c("iHMM","Triplet","Markov","GroundTruth"), colnames(d))) %>%
    mutate(category = paste(TT,trial_type,sep='')) %>%
    subset(category %in% c('HR','HP')) %>%
    group_by(participant_test, model, session_test, category) %>%
    summarise(rt_predicted_mean = mean(rt_predicted),
              rt_predicted_se = sd(rt_predicted)/sqrt(n()),
              rt_mean = mean(rt),
              rt_se = sd(rt)/sqrt(n())) %>%
    subset(!is.na(category))
  
  il_mean <- d %>%
    subset(model!="Triplet") %>%
    gather("variable","value", -(participant_test:category)) %>%
    unite(cat_var,category,variable) %>%
    spread(cat_var, value) %>%
    mutate(rt_predicted_diff_mean=HR_rt_predicted_mean-HP_rt_predicted_mean,
           rt_predicted_diff_se=sqrt(HR_rt_predicted_se^2+HP_rt_predicted_se^2),
           rt_diff_mean = HR_rt_mean-HP_rt_mean,
           rt_diff_se = sqrt(HR_rt_se^2+HP_rt_se^2))
  
  correlations <- il_mean %>%
    na.omit() %>%
    filter(session_test %in% c(1,2,3,4,5,6,7,8)) %>%
    group_by(model,session_test) %>%
    summarise(correlation = cor(rt_predicted_diff_mean,rt_diff_mean, use="complete.obs"),
              RMSE = sqrt(mean((rt_predicted_diff_mean-rt_diff_mean)^2, na.rm=TRUE)),
              cor_test = list(cor.test(rt_predicted_diff_mean,rt_diff_mean, use="complete.obs")))
  
  return(list("il_mean"=il_mean, "correlations"=correlations))
}

filter_and_average_pred_probs <- function(pred_probs){
  pred_probs %>%
    #subset(filters==TRUE) %>%
    group_by(participant_test, e_train, model, block, trial, Y, correct_response, first_response, session_test) %>%
    summarise(y0=mean(y0),y1=mean(y1),y2=mean(y2),y3=mean(y3)) %>%
    ungroup() %>%
    mutate(first_response_id = recode(
      first_response,
      z=0,c=1,b=2,m=3)
    ) %>%
    gather(key, prediction, (y0:y3))
}

compute_prediction_ranks <- function(data){
  data <- filter_and_average_pred_probs(data)
  data[data$model=='Triplet', 'prediction'] <- exp(data[data$model=='Triplet', 'prediction'])
  data[data$correct_response==0, 'response_type'] <- 'Incorrect\nResponse'
  data[data$correct_response==1, 'response_type'] <- 'Correct\nResponse'
  
  rank_data <- data %>%
    group_by(participant_test, e_train, model, block, trial, Y, correct_response, response_type, first_response) %>%
    mutate(pred_rank = rank(-prediction),
           c_model = paste0('c',correct_response,model)) %>%
    ungroup()
  
  stimulus_ranks <- rank_data %>%
    subset(paste0('y',Y)==key) %>%
    group_by(model, correct_response, response_type, c_model) %>%
    summarise(proportion = sum(pred_rank==1) / n(),
              success = sum(pred_rank==1),
              n = n())
  
  response_ranks <- rank_data %>%
    subset(paste0('y',first_response_id)==key) %>%
    subset(correct_response==0) %>%
    group_by(model, correct_response, response_type, c_model) %>%
    summarise(proportion = sum(pred_rank==1) / n(),
              success = sum(pred_rank==1),
              n = n())
  
  return(list("stimulus_ranks"=stimulus_ranks,
              "response_ranks"=response_ranks))
}

compute_roc_auc <- function(data){
  data <- data %>%
    add_session_info(., dataset="elarasztas") %>%
    filter_and_average_pred_probs(.) %>%
    subset(paste0('y',Y)==key)
  
  data[data$model=='Triplet', 'prediction'] <- exp(data[data$model=='Triplet', 'prediction'])
  
  roc_data_all <- data.frame()
  for (p in seq(0.0,1.0,0.01)){
    roc_data <- data %>%
      group_by(participant_test, e_train, model, session_test) %>%
      summarise(TPR = sum((correct_response==1) & (prediction>p)) / sum(correct_response==1),
                FPR= sum((prediction>p) & correct_response==0) / sum(correct_response==0))
    roc_data_all <- bind_rows(roc_data,roc_data_all)
  }
  
  auc <- roc_data_all %>%
    mutate(
      dTPR = c(diff(TPR),0),
      dFPR = c(diff(FPR),0)) %>%
    group_by(participant_test, e_train, model, session_test) %>%
    summarise(area=sum(dFPR*TPR)+sum(dFPR*dTPR)/2) %>%
    ungroup() %>%
    mutate(model = ordered(model, c('iHMM','GroundTruth','Markov','Triplet')))
  
  return(list(roc_data_all, auc))
}

across_sequences_performances <- function(data){
  data <- data %>%
    subset(correct_response==1) %>%
    subset(filters == TRUE) %>%
    group_by(participant_test, e_train, block, model, trial) %>%
    summarise(rt_predicted=mean(rt_predicted), rt=mean(rt)) %>%
    #mutate(epoch = cut(block, breaks=c(184,200,208,216,225,233,241,250))) %>%
    mutate(epoch = cut(block, breaks=c(184,195,200,205,210,215,220,225,230,235,240,245,250)),
           day_train_test = cut(block, breaks=c(184,195,200,210,220,225,230,235,240,245),
                                labels=c("Day 8 Train", "Day 8 Test", "Day 9 Early" , "Day 9 Train", "Day 9 Test", "Day 10 Old", "Day 10 New", "Day 10 Old", "Day 10 New"))) %>%
    subset(day_train_test %in% c("Day 8 Test","Day 9 Test","Day 10 Old", "Day 10 New")) %>% 
    group_by(participant_test, e_train, epoch, model, day_train_test) %>%
    mutate(cutoff = mean(rt) + 3 * sd(rt)) %>%
    subset(rt < cutoff) %>%
    summarise(performance = cor(rt,rt_predicted)) %>%
    mutate(performance = ifelse(model=="Triplet", -performance, performance)) %>%
    mutate(performance = performance^2*sign(performance))
}

relabel_and_order_models <- function(data){
  data %>%
    mutate(model=ordered(
      model,
      levels=c('iHMM','GroundTruth','Markov','Triplet','Data'),
      labels=c('CT','Ideal Observer','Markov','Trigram','Data')
    ))
}
