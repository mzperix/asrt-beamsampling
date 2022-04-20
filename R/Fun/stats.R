## Stats generation for the manuscript

source(here::here('R/Fun/data_pipelines.R'))
library(magrittr)
library(readr)
library(dplyr)
library(tidyr)
library(lsr)

######################
## REPORT FUNCTIONS ##
######################
t_test_with_effect_size <- function(x,y,alternative,paired,mu=NULL){
  if (is.null(mu)){
    t_test <- t.test(x,y,alternative=alternative, paired=paired)  
    cohend <- cohensD(x,y)
  }
  else {
    t_test <- t.test(x, alternative=alternative, mu=mu)
    cohend <- cohensD(x,mu=mu)
  }
  return(list("t_test"=t_test, "cohend"=cohend))
}

convert_day <- function(day){
  if(day==1){return('One')}
  if(day==8){return('Eight')}
}

format_p <- function(p_value){
  if (p_value < 0.001){
    return('p<0.001')
  }
  paste0('p=',format(p_value, digits=3))
}

format_t_ci <- function(t_test){
  paste0("CI=[",
         format(t_test[['conf.int']][1],digits=3),
         ",",
         format(t_test[['conf.int']][2],digits=3),
         "]")
}

t_test_report <- function(t_test_and_cohend, pretext=NULL){
  t_test_result <- t_test_and_cohend[['t_test']]
  alternative <- ifelse(t_test_result[['alternative']]=='two.sided','two-sided ','one-sided ')
  cohend <- t_test_and_cohend[['cohend']]
  paste0(pretext,
         alternative,
         "$t(",t_test_result[['parameter']],")=",
         format(t_test_result[['statistic']],digits=4),
         ",\\;",
         format_p(t_test_result[['p.value']]),
         ",\\;", 
         format_t_ci(t_test_result),
         ",\\;",
         "d=",
         format(cohend,digits=3),
         '$'
         )
}

correlation_report <- function(corr_result, pretext=NULL){
  paste0(pretext,
         "$r(",corr_result[['parameter']],")=",
         format(corr_result[['estimate']],digits=3),
         ",\\;",
         format_p(corr_result[['p.value']]),
         "$")
}

r2_report <- function(corr_result, pretext=NULL){
  paste0(pretext,
         "$R^2(",corr_result[['parameter']],")=",
         format(corr_result[['estimate']]^2,digits=3),
         ",\\;",
         format_p(corr_result[['p.value']]),
         "$")
}

binomial_test_report_from_test <- function(binom_test, pretext=NULL){
  n <- binom_test[["parameter"]]
  p <- binom_test[["statistic"]] / n
  paste0(pretext,
         "$", format(p, digits=3),",\\;",
         "n=",n,",\\;",
         format_p(binom_test[["p.value"]]), "$")
}

binomial_test_report <- function(success, n, p, alternative, pretext=NULL){
  binomial_test_report_from_test(binom.test(success,n,p,alternative=alternative))
}

repeated_measures_one_way_anova_report <- function(aov_result, pretext=NULL){
  dim1 <- summary(aov_result)[[1]][["Df"]][[1]]
  dim2 <- aov_result[["df.residual"]]
  F_stat <- summary(aov_result)[[1]][["F value"]][[1]]
  p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
  eta_partial <- etaSquared(aov_result)[1,2]
  paste0(pretext,
         "$F(",dim1,",",dim2,")=",format(F_stat, digits=4),",\\;",
         format_p(p_value),
         ",\\;",
         "\\eta_{partial}^2=", format(eta_partial,digits=3),
         "$")
}

latex_format <- function(variable_name, report){
  paste0("\\newcommand{\\", variable_name,"}{",report,"}")
}

############################
###    TESTS / VALUES    ###
############################

### METHODS / RT CUTOFF
rtCutoff <- function(){
  dataset <- read_csv("Data/elarasztas_dataset.csv") %>%
    subset(TrialType != "Prac")
  
  accurate <- dataset %>%
    subset(firstACC == 1)
  
  pct_accurate <- nrow(accurate) / nrow(dataset) * 100
  
  ge_180 <- subset(accurate, firstRT >= 180)
  pct_ge_180 <- nrow(ge_180) / nrow(accurate) * 100
  
  ge_180 <- subset(dataset, firstRT >= 180)
  pct_ge_180_overall <- nrow(ge_180) / nrow(dataset) * 100
  
  
  # ONLY DAY 8
  dataset <- dataset %>%
    filter(Block >= 176) %>%
    filter(Block <= 200)
  ge_180 <- subset(dataset, firstRT >= 180)
  pct_ge_180_overall <- nrow(ge_180) / nrow(dataset) * 100
  
  
  paste(c(len(dataset[dataset.])), sep="\n")
}

### LEARNING CURVES / RT PREDICTION

singleParticipantIhmmPerformanceDayEight <- function(rt_prediction_data){
  if(length(unique(rt_prediction_data$participant_test))==1){
    df <- rt_prediction_data %>%
      add_session_info(.,dataset='elarasztas') %>%
      subset(model=='iHMM') %>%
      subset(session_test==8) %>%
      filter_and_average_runs(sd_cutoff=100)
    cor.test(x=df$rt,y=df$rt_predicted, alternative="two.sided") %>%
      r2_report(.) %>%
      latex_format('singleParticipantIhmmPerformanceDayEight',.)
  }
}

ihmmPerformanceDayEight <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(model=='iHMM') %>%
    subset(session_test==8) %>%
    filter_and_average_runs(sd_cutoff=100) %>%
    group_by(participant_test) %>%
    summarise(R2=cor(rt,rt_predicted)^2)
  
  paste0("$M=",format(mean(df$R2), digits=3),
         "$ ranging from $",
         format(min(df$R2), digits=2), "$ to $",
         format(max(df$R2), digits=2), "$") %>%
    latex_format(variable_name='ihmmPerformanceDayEight')
}

markovPerformanceDayOne <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(model=='Markov') %>%
    subset(session_test==1) %>%
    filter_and_average_runs(sd_cutoff=100) %>%
    group_by(participant_test) %>%
    summarise(R2=cor(rt,rt_predicted)^2)
  
  paste0("$M=",format(mean(df$R2), digits=3),
         "$ ranging from $",
         format(min(df$R2), digits=2), "$ to $",
         format(max(df$R2), digits=2), "$") %>%
    latex_format(variable_name='markovPerformanceDayOne')
}

singleParticipantMarkovPerformanceDayEight <- function(rt_prediction_data){
  if(length(unique(rt_prediction_data$participant_test))==1){
    df <- rt_prediction_data %>%
      add_session_info(.,dataset='elarasztas') %>%
      subset(model=='Markov') %>%
      subset(session_test==8) %>%
      filter_and_average_runs(sd_cutoff=100)
    cor.test(x=df$rt,y=df$rt_predicted, alternative="two.sided") %>%
      r2_report(.) %>%
      latex_format(variable_name='singleParticipantMarkovPerformanceDayEight')
  }
}

tripletPerformanceSignificanceDayOne <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(model=='Triplet') %>%
    subset(session_test==1) %>%
    filter_and_average_runs(., 100, apply_filters=TRUE) %>%
    group_by(participant_test,e_train,session_test,model) %>%
    summarize(performance = cor(rt,rt_predicted))
  t_test_with_effect_size(x=df$performance, mu=0, alternative = "greater") %>%
    t_test_report(.) %>%
    latex_format(variable_name='tripletPerformanceSignificanceDayOne')
}

idealobserverPerformanceSignificanceDayOne <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(model=='GroundTruth') %>%
    subset(session_test==1) %>%
    filter_and_average_runs(., 100, apply_filters=TRUE) %>%
    group_by(participant_test,e_train,session_test,model) %>%
    summarize(performance = cor(rt,rt_predicted))
  t_test_with_effect_size(x=df$performance, mu=0, alternative = "greater") %>%
    t_test_report(.) %>%
    latex_format(variable_name='idealobserverPerformanceSignificanceDayOne')
}

markovPerformanceSignificanceDayOne <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(model=='Markov') %>%
    subset(session_test==1) %>%
    filter_and_average_runs(., 100, apply_filters=TRUE) %>%
    group_by(participant_test,e_train,session_test,model) %>%
    summarize(performance = cor(rt,rt_predicted))
  t_test_with_effect_size(x=df$performance, mu=0, alternative = "greater") %>%
    t_test_report(.) %>%
    latex_format(variable_name='markovPerformanceSignificanceDayOne')
}

model_comparison_on_day <- function(rt_prediction_data, model1, model2, day){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(session_test==day) %>%
    filter_and_average_runs(sd_cutoff=100) %>%
    group_by(model, participant_test) %>%
    summarise(r = cor(rt, rt_predicted)) %>%
    spread(.,"model","r")
  binom.test(sum(df[[model1]]>df[[model2]]), n=nrow(df), p=0.5, alternative="two.sided")
  
  #t_test_with_effect_size(x=df[[model1]], y=df[[model2]], paired=TRUE, alternative="two.sided")
}

model_comparison_on_day_with_text <- function(rt_prediction_data, model1, model2, day){
  model_comparison_on_day(rt_prediction_data, model1, model2, day) %>%
    binomial_test_report_from_test(., pretext="binomial test on MSE values ") %>%
    #t_test_report(., pretext='paired t-test on correlation values ') %>%
    latex_format(variable_name=paste0(model1,model2,'ComparisonDay',convert_day(day)))
}

iHMMTripletComparisonDayOne <- function(rt_prediction_data){
  model_comparison_on_day_with_text(rt_prediction_data, 'iHMM', 'Triplet', 1)
}

iHMMTripletComparisonDayEight <- function(rt_prediction_data){
  model_comparison_on_day(rt_prediction_data, 'iHMM', 'Triplet', 8)
}

iHMMMarkovComparisonDayOne <- function(rt_prediction_data){
  model_comparison_on_day_with_text(rt_prediction_data, 'iHMM', 'Markov', 1)
}

iHMMMarkovComparisonDayEight <- function(rt_prediction_data){
  model_comparison_on_day(rt_prediction_data, 'iHMM', 'Markov', 8)
}

iHMMGroundTruthComparisonDayEight <- function(rt_prediction_data){
  model_comparison_on_day(rt_prediction_data, 'iHMM', 'GroundTruth', 8)
}

iHMMGroundTruthComparisonDayEight_text <- function(rt_prediction_data){
  model_comparison_on_day_with_text(rt_prediction_data, 'iHMM', 'GroundTruth', 8)
}


iHMMControlComparisonDayEight <- function(rt_prediction_data){
  # Three comparisons with each of the controls
  paste(
    iHMMTripletComparisonDayEight(rt_prediction_data) %>% binomial_test_report_from_test(.,pretext='CT vs Trigram '),
    iHMMMarkovComparisonDayEight(rt_prediction_data)  %>% binomial_test_report_from_test(.,pretext='CT vs Markov '),
    iHMMGroundTruthComparisonDayEight(rt_prediction_data) %>% binomial_test_report_from_test(.,pretext='CT vs Ground Truth '),
    sep=' ') %>%
    
    latex_format(variable_name = 'iHMMControlComparisonDayEight')
}

iHMMMeanPerformanceWithSEMDayEight <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(session_test==8) %>%
    subset(model=='iHMM') %>%
    compute_performance(.) %>%
    ungroup() %>%
    summarise(mean_r2 = mean(performance),
              se = sd(performance)/sqrt(n()))
  paste0("$R^2=",
         format(df$mean_r2, digits=3),
         "$ (s.e.m. $", format(df$se, digits=3),"$)") %>%
    latex_format(variable_name="iHMMMeanPerformanceWithSEMDayEight")
}

markovMeanPerformanceDayThree <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(session_test==3) %>%
    subset(model=='Markov') %>%
    compute_performance(., sd_cutoff=100) %>%
    ungroup() %>%
    summarise(mean_r2 = mean(performance),
              se = sd(performance)/sqrt(n()))
  paste0("$R^2=",
         format(df$mean_r2, digits=3),
         "$") %>%
    latex_format(variable_name="markovMeanPerformanceDayThree")
}

markovMeanPerformanceDayEight <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(.,dataset='elarasztas') %>%
    subset(session_test==8) %>%
    subset(model=='Markov') %>%
    compute_performance(., sd_cutoff=100) %>%
    ungroup() %>%
    summarise(mean_r2 = mean(performance),
              se = sd(performance)/sqrt(n()))
  paste0("$R^2=",
         format(df$mean_r2, digits=3),
         "$") %>%
    latex_format(variable_name="markovMeanPerformanceDayEight")
}

ihmmMarkovPlusGroundTruthCorrelationDayEight <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(., dataset="elarasztas") %>%
    subset(session_test==8) %>%
    compute_ihmm_markov_io_diff(.)
  cor.test(x=df$iHMM-df$Markov, y=df$GroundTruth, method="pearson") %>%
    correlation_report() %>%
    latex_format(variable_name="iHMMMarkovPlusGroundTruthCorrelationDayEight")
}

ihmmMarkovPlusGroundTruthCorrelationAllDays <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(., dataset="elarasztas") %>%
    compute_ihmm_markov_io_diff(.)
  cor.test(x=df$iHMM-df$Markov, y=df$GroundTruth, method="pearson") %>%
    correlation_report() %>%
    latex_format(variable_name="iHMMMarkovPlusGroundTruthCorrelationAllDays")
}

ihmmMarkovPlusGroundTruthComparisonDayEight <- function(rt_prediction_data){
  df <- rt_prediction_data %>%
    add_session_info(., dataset="elarasztas") %>%
    subset(session_test==8) %>%
    compute_ihmm_markov_io_diff(.)
  t_test_with_effect_size(x=df$diff, mu=0, alternative="greater") %>%
    t_test_report(.) %>%
    latex_format(variable_name="iHMMMarkovPlusGroundTruthComparisonDayEight")
}

markovPerformanceDecline <- function(rt_prediction_data){
  from <- 3
  to <- 8
  aov_data <- rt_prediction_data %>%
    add_session_info(., dataset="elarasztas") %>%
    subset(session_test %in% seq(from,to,1)) %>%
    subset(model=='Markov') %>%
    mutate(participant_test = as.factor(participant_test)) %>%
    compute_performance(., sd_cutoff=100) %>%
    ungroup()
  aov_result <- aov(performance ~ session_test + participant_test, data=aov_data)
  aov_result %>%
    repeated_measures_one_way_anova_report() %>%
    latex_format(variable_name="markovPerformanceDecline")
}

ihmmMarkovSpreadLinear <- function(rt_prediction_data, participant){
  aov_data <- rt_prediction_data %>%
    add_session_info(., dataset="elarasztas") %>%
    subset(model %in% c('iHMM','Markov')) %>%
    mutate(participant_test = as.factor(participant_test)) %>%
    compute_performance(., sd_cutoff=100) %>%
    ungroup() %>%
    spread(model, performance) %>%
    mutate(diff=iHMM-Markov)
  aov_result <- aov(diff ~ session_test + participant_test, data=aov_data)
  aov_result %>%
    repeated_measures_one_way_anova_report() %>%
    latex_format(variable_name="iHMMMarkovSpreadLinear")
}

learning_curve_stats <- function(rt_prediction_data, participant){
  paste(
    singleParticipantIhmmPerformanceDayEight(subset(rt_prediction_data,participant_test==participant)),
    ihmmPerformanceDayEight(rt_prediction_data),
    singleParticipantMarkovPerformanceDayEight(subset(rt_prediction_data,participant_test==participant)),
    tripletPerformanceSignificanceDayOne(rt_prediction_data),
    idealobserverPerformanceSignificanceDayOne(rt_prediction_data),
    markovPerformanceSignificanceDayOne(rt_prediction_data),
    iHMMTripletComparisonDayOne(rt_prediction_data),
    iHMMMarkovComparisonDayOne(rt_prediction_data),
    iHMMGroundTruthComparisonDayEight_text(rt_prediction_data),
    iHMMMeanPerformanceWithSEMDayEight(rt_prediction_data),
    markovPerformanceDayOne(rt_prediction_data),
    markovMeanPerformanceDayThree(rt_prediction_data),
    markovMeanPerformanceDayEight(rt_prediction_data),
    ihmmMarkovPlusGroundTruthCorrelationDayEight(rt_prediction_data),
    ihmmMarkovPlusGroundTruthCorrelationAllDays(rt_prediction_data),
    ihmmMarkovPlusGroundTruthComparisonDayEight(rt_prediction_data),
    iHMMControlComparisonDayEight(rt_prediction_data),
    ihmmMarkovSpreadLinear(rt_prediction_data),
    markovPerformanceDecline(rt_prediction_data),
    sep='\n'
  )
}

### ERROR PREDICTION
stimulus_rank_tests <- function(pred_probs){
  ranks <- pred_probs %>% 
    add_session_info(., dataset="elarasztas") %>% 
    subset(session_test==8) %>%
    #filter_and_average_pred_probs(.) %>% 
    compute_prediction_ranks(.)
  sr <- ranks[['stimulus_ranks']]
  paste(
    latex_format(variable_name='ihmmCorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="iHMM") & (sr$correct_response==1),]$success,
                                      sr[(sr$model=="iHMM") & (sr$correct_response==1),]$n, 
                                      0.25, "greater")),
    latex_format(variable_name='ihmmIncorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="iHMM") & (sr$correct_response==0),]$success,
                                      sr[(sr$model=="iHMM") & (sr$correct_response==0),]$n, 
                                      0.25, "less")),
    latex_format(variable_name='markovCorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="Markov") & (sr$correct_response==1),]$success,
                                      sr[(sr$model=="Markov") & (sr$correct_response==1),]$n, 
                                      0.25, "greater")),
    latex_format(variable_name='markovIncorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="Markov") & (sr$correct_response==0),]$success,
                                      sr[(sr$model=="Markov") & (sr$correct_response==0),]$n, 
                                      0.25, "less")),
    latex_format(variable_name='tripletCorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="Triplet") & (sr$correct_response==1),]$success,
                                      sr[(sr$model=="Triplet") & (sr$correct_response==1),]$n, 
                                      0.25, "greater")),
    latex_format(variable_name='tripletIncorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="Triplet") & (sr$correct_response==0),]$success,
                                      sr[(sr$model=="Triplet") & (sr$correct_response==0),]$n, 
                                      0.25, "less")),
    latex_format(variable_name='idealobserverCorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="GroundTruth") & (sr$correct_response==1),]$success,
                                      sr[(sr$model=="GroundTruth") & (sr$correct_response==1),]$n, 
                                      0.25, "greater")),
    latex_format(variable_name='idealobserverIncorrectStimulusRankTest', 
                 binomial_test_report(sr[(sr$model=="GroundTruth") & (sr$correct_response==0),]$success,
                                      sr[(sr$model=="GroundTruth") & (sr$correct_response==0),]$n, 
                                      0.25, "greater")),
    
    sep="\n")
}

auc_test <- function(pred_probs, model1, model2){
  roc_auc <- pred_probs %>%
    subset(model %in% c(model1, model2)) %>%
    add_session_info(., dataset="elarasztas") %>%
    #subset(session_test==8) %>%
    compute_roc_auc(.)
  auc <- roc_auc[[2]] %>%
    spread(model, area)
  
  t_test_with_effect_size(x=auc[[model1]], y=auc[[model2]], paired=TRUE, alternative="greater") %>%
    t_test_report(.,pretext='paired t-test on AUC values ') %>%
    latex_format(variable_name=paste0(model1,model2,'AUCTest'))
}

auc_tests <- function(pred_probs){
  paste(
    auc_test(pred_probs, "iHMM", "Markov"),
    auc_test(pred_probs, "iHMM", "Triplet"),
    auc_test(pred_probs, "iHMM", "GroundTruth"),
    sep="\n")
}

error_prediction_tests <- function(pred_probs){
  paste(auc_tests(pred_probs),
        stimulus_rank_tests(pred_probs),
        sep='\n')
}

### ACROSS PARTICIPANT TESTS
across_participant_test <- function(data_across_participants, group1, group2, variable_name){
  # Groups can be iHMM/Markov + Permuted/NotPermuted, e.g. iHMMPermuted
  data <- data_across_participants %>%
    subset(participant_test!=participant_train) %>%
    mutate(
      performance=correlation^2*sign(correlation),
      permuted=ifelse(permuted,'Permuted','NotPermuted'),
      group=paste0(model,permuted)) %>%
    select(participant_test,participant_train, e_test, e_train, group, performance) %>%
    spread(group, performance)
  
  t_test_with_effect_size(x=data[[group1]],y=data[[group2]], paired=TRUE, alternative="two.sided") %>%
    t_test_report() %>%
    latex_format(variable_name=variable_name)
}

across_participant_tests <- function(data_across_participants){
  paste(
    across_participant_test(data_across_participants, 'iHMMPermuted', 'iHMMNotPermuted', 'iHMMAcrossParticipantPermutedNotPermuted'),
    across_participant_test(data_across_participants, 'MarkovNotPermuted', 'iHMMNotPermuted', 'iHMMMarkovAcrossParticipantNotPermutedTest'),
    across_participant_test(data_across_participants, 'iHMMPermuted', 'MarkovPermuted', 'iHMMMarkovAcrossParticipantPermutedTest'),
    sep='\n')
}

### FINGERPRINT 

### Higher Order stats
higher_order_stats <- function(rt_prediction_data){
  ho_stats <- compute_higher_order_stats(rt_prediction_data)
  cors <- ho_stats[["correlations"]]
 
  paste(
    cors[(cors$model=='iHMM') & (cors$session_test==2),][['cor_test']][[1]] %>% correlation_report(.) %>% latex_format(variable_name="ihmmHigherOrderCorrelationDayTwo"),
    cors[(cors$model=='iHMM') & (cors$session_test==8),][['cor_test']][[1]] %>% correlation_report(.) %>% latex_format(variable_name="ihmmHigherOrderCorrelationDayEight"),
    cors[(cors$model=='GroundTruth') & (cors$session_test==2),][['cor_test']][[1]] %>% correlation_report(.) %>% latex_format(variable_name="ioHigherOrderCorrelationDayTwo"),
    cors[(cors$model=='GroundTruth') & (cors$session_test==8),][['cor_test']][[1]] %>% correlation_report(.) %>% latex_format(variable_name="ioHigherOrderCorrelationDayEight"),
    sep="\n"
  ) 
}

### Across sequences
ihmm_across_sequence_comparison <- function(across_sequences_data, e_train1, e_train2, session_test1, session_test2, variable_name){
  test_data <- across_sequences_data %>%
    add_session_info(dataset="elarasztas") %>%
    subset(model=="iHMM") %>%
    subset(session_test %in% c(session_test1, session_test2)) %>%
    subset(e_train %in% c(e_train1, e_train2)) %>%
    compute_performance(sd_cutoff=100) %>%
    mutate(group=paste(e_train, session_test, sep='_')) %>%
    ungroup() %>%
    select(participant_test, model, group, performance) %>%
    spread(group, performance)
  
  col1 <- paste(e_train1, session_test1, sep='_')
  col2 <- paste(e_train2, session_test2, sep='_')
  t_test_with_effect_size(test_data[[col1]],test_data[[col2]],alternative="greater", paired=TRUE) %>%
    t_test_report() %>%
    latex_format(variable_name=variable_name) %>%
    return()
}

day10_comparisons <- function(across_sequences_data, group1, group2, variable_name){
  test_data <- across_sequences_data %>%
    across_sequences_performances() %>%
    subset(model == 'iHMM') %>%
    group_by(participant_test, model, e_train, day_train_test) %>%
    summarise(performance = mean(performance)) %>%
    mutate(group = paste0(e_train, day_train_test)) %>%
    ungroup() %>%
    select(participant_test, group, performance) %>%
    spread(group, performance)
  t_test_with_effect_size(test_data[[group1]], test_data[[group2]], alternative="greater", paired=TRUE) %>%
    t_test_report() %>%
    latex_format(variable_name=variable_name)
}

across_sequence_comparisons <- function(across_sequences_data){
  paste(
    ihmm_across_sequence_comparison(across_sequences_data, '186_195','186_195',8,9,'ihmmDayEightToDayNineDifference'),
    ihmm_across_sequence_comparison(across_sequences_data, '211_220','211_220',9,8,'ihmmDayNineToDayEightDifference'),
    day10_comparisons(across_sequences_data, '186_195Day 10 Old','211_220Day 10 Old', 'ihmmDayTenOldSequenceComparison'),
    day10_comparisons(across_sequences_data, '211_220Day 10 New','186_195Day 10 New', 'ihmmDayTenNewSequenceComparison'),
    day10_comparisons(across_sequences_data, '186_195Day 10 Old','186_195Day 10 New', 'ihmmDayEightOldNewSequenceComparison'),
    day10_comparisons(across_sequences_data, '211_220Day 10 New','211_220Day 10 Old', 'ihmmDayNineOldNewSequenceComparison'),
    
    sep='\n'
  )
}
