adjust_fingerprint_data <- function(data){
  df <- data %>%
    mutate(sequence = as.factor(sequence)) %>%
    #subset(rt_count >= 5) %>%
    group_by(participant_test, e_train, sequence) %>%
    summarise(predicted_probability=
                exp(mean(log_predicted_probability, na.rm=TRUE)),
              rt_predicted = mean(rt_predicted),
              rt_mean = mean(rt_mean),
              rt_sem = mean(rt_sem),
              rt_count = mean(rt_count))
  
  return(df)
}

fig_fingerprint_across_sessions <- function(df){
  data <- adjust_fingerprint_data(df)
  f1 <- data %>%
    ggplot(aes(x=sequence, y=predicted_probability, fill=e_train)) +
    geom_bar(stat="identity") +
    facet_grid(rows=vars(participant_test)) +
    style +
    theme(axis.text.x = element_text(angle=90))
  data <- data %>%
    select(participant_test, sequence, e_train, predicted_probability) %>%
    spread(key=e_train, value=predicted_probability) %>%
    print(head(data))
  f2 <- data %>%
    ggplot(aes(x=Model1, y=Model2,
               text=paste("sequence:",sequence))) +
    geom_point(color=pantone_2014[["dark_blue"]], size=1) +
    facet_wrap(vars(participant_test), ncol=5) +
    coord_fixed() +
    style
  return(f2)
}

fig_fingerprint_predicted_rt <- function(df, min_rt_count=7){
  data <- adjust_fingerprint_data(df) %>%
    subset(!is.na(rt_sem)) %>%
    subset(rt_sem>0) %>%
    subset(rt_count >= min_rt_count)
  
  cors <- data %>%
    group_by(participant_test) %>%
    summarise(correlation = round(cor(rt_predicted, rt_mean), 3),
              x = min(rt_predicted),
              y = max(rt_mean),
              N = n())
  
  f3 <- data %>%
    ggplot(aes(x=rt_predicted, y=rt_mean)) +
    geom_abline(slope = 1.0, intercept = 0.0, color = pantone_2014[['grey']]) +
    geom_errorbar(aes(ymin=rt_mean-2*rt_sem, ymax=rt_mean+2*rt_sem),
                  width=0.07, color=pantone_2014[["light_grey"]], size=0.6) +
    geom_point(color=pantone_2014[["dark_blue"]], shape=20) +
    geom_text(data=cors, 
              aes(label=paste("r=", correlation, sep=""),
                  x=220, y=480),
              hjust=0) +
    geom_text(data=cors, 
              aes(label=paste("n=", N, sep=""),
                  x=220, y=460),
              hjust=0) +
    facet_wrap(vars(participant_test), ncol=5) +#, scales="free") +
    xlab("Predicted Response Time") +
    ylab("Mean Response Time") +
    coord_fixed() +
    ylim(180,500) +
    style
  return(f3)
}

fig_fingerprint_barplot <- function(df){
  data <- adjust_fingerprint_data(df)
  print(head(data))
  data <- data %>% 
    #subset(participant_test==101) %>%
    select(participant_test, sequence, e_train, predicted_probability) %>%
    mutate(predicted_probability=-log(predicted_probability)) %>%
    spread(key=e_train, value=predicted_probability)
  #data$sequence <- reorder(data$sequence, data$Model1)
  data <- data %>%
    arrange(participant_test, Model1) %>%
    # 3. Add order column of row numbers
    mutate(order = row_number())
  data %>%
    ggplot(aes(x=order, y=Model1)) +
    geom_bar(aes(y=Model2), stat="identity", 
             fill=pantone_2014[["light_blue"]]) +
    geom_bar(stat="identity", color = pantone_2014[["dark_blue"]], 
             fill=alpha("white",0.0)) +
    style +
    facet_wrap(vars(participant_test), ncol=2) +
    theme(axis.text.x = element_text(angle=90))
}