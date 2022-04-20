library(magrittr)
library(ggplot2)
library(dplyr)

setwd("/Users/balazstorok/SCIENCE/Projects/")
setwd("reactive-inhibition/Doc")
source("../R/Fun/getData.R")
Data.All <- getASRTData("180")

n_Subjects <- 80
Subjects <- unique(Data.All$Subject)[1:n_Subjects]

data <- Data.All %>%
  subset(Subject %in% Subjects) %>%
  group_by(Subject) %>%
  mutate(N = 1:n()) %>%
  subset(N < 1531) %>%
  subset(N > 1275)

data[data$TT == "H","logP"] <- -log(0.625)
data[data$TT == "L","logP"] <- -log(0.125)

data <- subset(data, logP > -Inf)
data <- subset(data, firstRT > 0)

data %>%
  group_by(.,Subject)%>%
  summarise(Corr = cor(-logP, 1/firstRT)) %>%
  ggplot(aes(x = Subject, y = Corr)) + 
  geom_point()

data %>%
  group_by(.,Subject)%>%
  summarise(Corr = cor(-logP, 1/firstRT)) %>%
  ungroup() %>%
  summarise(result = mean(Corr))
