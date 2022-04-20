# Plot cross_entropy values

library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)

# here::set_here("~/Science/asrt-beamsampling")

source(here::here("R/Fun/plot/plot_styles.R"))

cross_entropy <- function(df, model1, model2){
  if (model1 == "Generative"){
    df1 <- subset(df, select=c("model", "participant_test", "block", "trial", "trial_type", "Y", "y0", "y1", "y2", "y3")) %>%
      generating_probabilities()
  }
  else {
    df1 <- subset(df, model==model1, select=c("model", "participant_test", "block", "trial", "y0", "y1", "y2", "y3"))
  }
  df2 <- subset(df, model==model2, select=c("model", "participant_test", "block", "trial", "y0", "y1", "y2", "y3"))
  
  result <- merge(df1, df2, by=c("participant_test", "block", "trial")) %>%
    mutate(cross_entropy=y0.x*log(y0.y+1e-16)+y1.x*log(y1.y+1e-16)+y2.x*log(y2.y+1e-16)+y3.x*log(y3.y+1e-16))
  return(-mean(result$cross_entropy))
}

generating_probabilities <- function(df) {
  df %>% mutate(
    y0 = ifelse(trial_type=="P", as.integer(Y == 0), 0.25),
    y1 = ifelse(trial_type=="P", as.integer(Y == 1), 0.25),
    y2 = ifelse(trial_type=="P", as.integer(Y == 2), 0.25),
    y3 = ifelse(trial_type=="P", as.integer(Y == 3), 0.25)
    )
}

prepare_cross_entropy_data <- function(df){
  df <- df %>%
    group_by(participant_test, block, trial, model, e_test, e_train, trial_type, Y) %>%
    summarise(y0=mean(y0), y1=mean(y1), y2=mean(y2), y3=mean(y3)) %>%
    ungroup()

  cross_entropy_with_ihmm <- data.frame()
  for (p in unique(df$participant_test)){
    for (epoch_test in unique(df$e_test)){
      for (model in c("Generative", "Markov", "GroundTruth")){
        cross_entropy_with_ihmm <- rbind(
          cross_entropy_with_ihmm, 
          data.frame(
            participant_test=p, 
            null_model=model, 
            e_test=epoch_test,
            cross_entropy=cross_entropy(
              subset(
                df, 
                (participant_test==p) & (e_test==epoch_test)
              ), 
              model, "iHMM"
              )
            )
          )
      } 
    }
  }

  return (cross_entropy_with_ihmm)
}

cross_entropy_plot <- function(cross_entropy_with_ihmm, e_test, participants_to_highlight){
  plot_data <- cross_entropy_with_ihmm %>%
  spread(key="null_model", value="cross_entropy") %>%
  subset(e_test==e_test) 

  plot_data %>%
    ggplot(aes(x=Generative, y=Markov)) +
    geom_point(color=pantone_2014[["darker_blue"]]) +
    geom_point(data=plot_data[which(plot_data$participant_test %in% participants_to_highlight),], color=pantone_2014[["light_blue"]]) +
    style + 
    theme(panel.border = element_rect(size = 1.0, fill=NA,
                                      colour = "#747474"), 
          axis.line = element_line(size = 0, color = "#747474"))+
    xlab(expression(paste(D[KL],"(Generative Probabilities || CT)"), parse=TRUE)) +
    ylab(expression(paste(D[KL],"(Markov || CT)"), parse=TRUE)) +
    xlim(1.25, 1.55) +
    ylim(0.85, 1.4)
}

fig6gh <- function(data, epoch_test, participants_to_highlight){
  cross_entropy_with_ihmm <- data %>%
    filter(e_test == epoch_test) %>%
    prepare_cross_entropy_data()

  return(cross_entropy_plot(cross_entropy_with_ihmm, e_test, participants_to_highlight))
}
