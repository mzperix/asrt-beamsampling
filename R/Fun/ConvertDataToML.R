library(readr)
library(magrittr)
library(dplyr)

setwd("/Users/balazstorok/SCIENCE/Projects/asrt-beamsampling/Data")
ASRT_data <- read_csv("ASRT_180_Matlab.csv", col_names = FALSE)
colnames(ASRT_data) <- c("Subject", "Session", "Block", "event", "RT", "ACC", "key", "TrialType","TT", "QT", "QiT")

files <- c("ASRT1_ML_text.csv",
           "ASRT2_ML_text.csv",
           "ASRT3_ML_text.csv")

i = 1
for (f in files)
{
  data <- subset(ASRT_data, Session == i)
  colnames(data) <- c("Subject", "Session", "Block", "event", "firstRT", "firstACC", "firstRESP", "TrialType")
  
  data <- data %>%
    subset(select = c("Subject", "Session", "Block", "event", "firstRT", "firstACC", "firstRESP", "TrialType")) %>%
    mutate(firstRESP = 0, TrialType = 0)
  
  write_csv(data, paste(substr(f,1,8),".csv",sep = ""), col_names = FALSE)
  i <- i+1
}
