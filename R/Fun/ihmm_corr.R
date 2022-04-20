library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(tidyr)
library(cowplot)
setwd("/Users/balazstorok/SCIENCE/Projects/asrt-beamsampling/R/Docs")
source("../Fun/getData.R")
data_ihmm <- getIHMMResult("mem-1-6")

setwd("/Users/balazstorok/SCIENCE/Projects/InverseLearning")
ASRT_data <- read_csv("ASRT_180_Matlab.csv", col_names = FALSE)
setwd("/Users/balazstorok/SCIENCE/Projects/asrt-beamsampling/R/Docs")

### Check nothing suspicious with results ###
samples <- sample(x = unique(data_ihmm$N), replace = FALSE, size = 500)
p1 <- data_ihmm %>%
  subset(N %in% samples) %>%
  subset(calcpoint == 2000) %>%
  gather(key = "memory", value = "PPred", -Subject, -calcpoint, -N) %>%
  ggplot(aes(x = N, y = PPred, color = memory)) +
  geom_point(size = 0.5) +
  facet_wrap("Subject", ncol = 10)

ggsave("Subject_101-353_iHMM_predictions.pdf", 
       plot   = p1,
       scale  = 1, 
       width  = 27, 
       height = 3*length(unique(data_ihmm$Subject))/10, 
       units = c("cm"), 
       dpi = 300)

### FILTER ERRONOUS SUBJECTS: 205, 228, 298, 321
Subject_error <- c(205, 228, 298, 321)
data_ihmm <- subset(data_ihmm, !(Subject %in% Subject_error))

### RESIDUALS ###

RT_data <- read_csv("../../Data/Asrt_180_Matlab.csv")
colnames(RT_data) <- c("Subject", "Session", "Block", "event", "RT", "firstACC", "button", "TrialType", "TT", "QT", "QiT")
RT_data <- RT_data %>%
  group_by(Subject) %>%
  mutate(N = 1:n())


source("../Fun/getData.R")
source("../Fun/ModelDefinitions.R")
source("../Fun/Residual.R")
fp <- getfitResult("allRT")

fl.testset <- list.files(path = "../../Data/", pattern = "ALL")
testset <- read_csv(paste("../../Data/", fl.testset, sep = ""), col_names = TRUE)
RT_data$firstRT <- RT_data$RT
residual <- calc.res.cv(RT_data, fp, testset, Model="1110110")
res_plot <- residual %>%
  ggplot(aes(x=N, y = resRT)) +
  geom_point(size = 0.2) +
  ylim(-100,100) +
  facet_grid(Subject~.)

Data <- residual %>%
  subset(select = c("Subject", "N", "resRT", "RT")) %>%
  merge(data_ihmm, ., by = c("Subject", "N"), all.y = FALSE)

p2 <- Data %>%
  subset(calcpoint == 2000) %>%
  gather(key = "memory", value = "PPred", -Subject, -calcpoint, -N, -resRT) %>%
  ggplot(aes(x = log(PPred), y = resRT, color = memory)) +
  geom_point(size = 0.5) +
  facet_wrap("Subject", ncol = 5)

ggsave("Subject_101-353_iHMM_resRT.pdf", 
       plot   = p2,
       scale  = 1, 
       width  = 10, 
       height = 49, 
       units = c("in"), 
       dpi = 300)


p4 <- Data %>%
  #mutate(Subject = as.factor(Subject)) %>%
  subset(calcpoint == 2000) %>%
  subset(N > 1275) %>%
  #subset(Subject < 120) %>%
  gather(., key = "memory", value = "PPred", -Subject, -calcpoint, -N, -resRT) %>%
  mutate(log_P    = log(PPred)) %>%
  mutate(P_bin    = cut(PPred, breaks = seq(0,1,0.1)),
         logP_bin = cut(log_P, breaks = seq(-10,10,0.2))) %>%
  group_by(Subject, calcpoint, memory, logP_bin) %>%
  summarise(medResRT  = median(resRT),
            meanResRT = mean(resRT),
            n = n()) %>%
  subset(n > 100) %>%
  ggplot(aes(x = logP_bin, y = medResRT, alpha = log(n), color = memory)) +
  geom_point() +
  facet_wrap("Subject", ncol = 8) +
  theme(axis.text.x = element_text(size = 9, angle = 90)) +
  ylim(-50, 50)
#p4

ggsave("Subject_101-353_iHMM_medianresRT.pdf", 
       plot   = p4,
       scale  = 1, 
       width  = 10, 
       height = 49, 
       units = c("in"), 
       dpi = 300)

p_correlations <- Data %>%
  subset(calcpoint == 2000) %>%
  subset(N > 1275) %>%
  group_by(Subject) %>%
  summarise(corr1 = cor(RT, log(mem1)),
         corr2 = cor(RT, log(mem2)),
         corr3 = cor(RT, log(mem3)),
         corr4 = cor(RT, log(mem4)),
         corr5 = cor(RT, log(mem5)),
         corr6 = cor(RT, log(mem6))) %>%
  gather(key = "mem", value = "corr", -Subject) %>%
  ggplot(aes(x = corr)) +
  geom_histogram() +
  facet_grid(mem~.)

Subjects <- Data %>%
  subset(calcpoint == 2000) %>%
  subset(N > 2550) %>%
  group_by(Subject) %>%
  summarise(corr6 = cor(resRT, log(mem6))) %>%
  subset(corr6 < -0.12)
Subjects <- unique(Subjects$Subject)

p5 <- Data %>%
  subset(calcpoint == 2000) %>%
  subset(N > 1275) %>%
  subset(Subject %in% Subjects) %>%
  ggplot(aes(y = log(mem6), x = 1/(300+resRT))) +
  geom_point(size = 0.2) +
  facet_grid(Subject ~.)


ggsave("iHMM_medianresRT_best_corr.pdf", 
       plot   = p5,
       scale  = 1, 
       width  = 10, 
       height = 30, 
       units = c("in"), 
       dpi = 300)


### CORRELATION FOR HIGH PROBABILITY EVENTS ###
p_corr_high <- Data %>%
  merge(., subset(residual, select = c("Subject", "N", "TrialType")), all.y = FALSE) %>%
  #subset(mem5 > 0.8) %>%
  subset(N > 1275) %>%
  subset(TrialType == "P") %>%
  group_by(Subject) %>%
  summarise(corr1 = cor(RT, log(mem1)),
            corr2 = cor(RT, log(mem2)),
            corr3 = cor(RT, log(mem3)),
            corr4 = cor(RT, log(mem4)),
            corr5 = cor(RT, log(mem5)),
            corr6 = cor(RT, log(mem6))) %>%
  gather(key = "mem", value = "corr", -Subject) %>%
  ggplot(aes(x = corr)) +
  geom_histogram() +
  facet_grid(mem~.)
  