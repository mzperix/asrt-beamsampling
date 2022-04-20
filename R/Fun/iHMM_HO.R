### iHMM predictions and Higher Order sequence learning ###
library(readr)
library(ggplot2)
library(dplyr)
library(magrittr)

setwd("/Users/balazstorok/SCIENCE/Projects/InverseLearning")
ASRT_data <- read_csv("ASRT_180_Matlab.csv", col_names = FALSE)
colnames(ASRT_data) <- c("Subject", "Session", "Block", "event", "RT", "ACC", "key", "TrialType","TT", "QT", "QiT")

Subjects <- unique(ASRT_data$Subject)
mem <- c(1,2,3,4,5,6)

Subjects <- c(101:108)

setwd("/Users/balazstorok/SCIENCE/Projects/asrt-beamsampling/Results")
iHMM_data = data.frame()
for (s in Subjects)
{
  for (m in mem)
  {
    dat <- data.frame(PPred = read_csv(paste("Subject-",s,"_distrMean_mem-",m,".csv",sep = ""), 
                                             col_names = FALSE,
                                             col_types = list(col_double())))
    colnames(dat) <- c("PPred")
    dat <- mutate(dat, Subject = s, mem = m, N = 1:n())
    iHMM_data <- rbind(iHMM_data, dat)
  }
  remove(dat)
}

Data <- ASRT_data %>%
  subset(Subject %in% Subjects) %>%
  group_by(Subject) %>%
  mutate(N = 1:n()) %>%
  merge(., iHMM_data, by = c("Subject", "N"))

p1 <- Data %>% 
  mutate(mem = as.factor(mem)) %>%
  ggplot(aes(x = N, y = PPred, color = mem)) +
  geom_point(size = 0.2) + 
  facet_wrap(c("Subject", "mem"), ncol = 6)

ggsave("Subject_101-108_iHMM_PPred.pdf", 
       plot   = p1,
       path = "../Figures",
       scale  = 1, 
       width  = 27, 
       height = 9*length(Subjects), 
       units = c("cm"), 
       dpi = 300)

### HO for iHMM ###
HO_iHMM <- Data %>%
  mutate(mem = as.factor(mem)) %>%
  group_by(Subject, mem) %>%
  mutate(Block = Block + (Session-1) *15) %>%
  subset(ACC == 1) %>%
  subset(TT %in% c("H","L")) %>%
  mutate(Type = paste(TrialType, TT, sep = "")) %>%
  group_by(Subject, Session, Block, Type, mem) %>%
  summarise(PPred.mean = mean(PPred))

p2 <- HO_iHMM %>%
  ggplot(aes(x = Block, y = PPred.mean, color = Type)) +
  geom_line() +
  facet_wrap(c("Subject", "mem"), ncol = 6)

ggsave("Subject_101-108_iHMM_HO.pdf", 
       plot   = p2,
       path = "../Figures",
       scale  = 1, 
       width  = 27, 
       height = 9*length(Subjects), 
       units = c("cm"), 
       dpi = 300)
