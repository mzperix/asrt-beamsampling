### Functions for Loading in various sorts of data

library(readr)
library(magrittr)
library(dplyr)
library(tidyr)

Cols <- c("Subject", "Session", "Block", "event",
          "firstRT", "firstACC", "firstRESP", "TrialType",
          "TT", "QT", "QiT")#,
          #"filtered")

getASRTData <- function(File = "")
{
  datadir <- "../Data/ASRT/"
  data <- list.files(path = datadir, pattern=paste("*",File,"_*", sep = "")) %>%
    lapply(function(x) read_csv(paste(datadir,x,sep = ""), col_names = FALSE)) %>%
    rbind_all()
  colnames(data) <- Cols
  
  return(data)
}

getfitResult <- function(Subject = "")
{
  ### LOAD ###
  datadir <- "../../Data/"
  fl <- list.files(path = datadir, pattern = paste("*",Subject,"_*",sep = ""))
  fp <- fl %>%
    lapply(function(x) {
      cat(x, "\n")
      read_csv(paste(datadir,x,sep = ""), 
               col_names = TRUE) %>%
        subset(.,select = 1:25)
      }) %>%
    bind_rows()
  
  fp <- fp %>% mutate(model = as.character(A*1000000+B*100000+C*10000+D*1000+E*100+F*10+G))
  cat("Loaded ",length(fl)," files \n", length(fl))
  return(fp)
}

getIHMMResult <- function(title = "")
{
  datadir <- "../../Results/"
  fl <- list.files(path = datadir, pattern = paste("*",title,"*",sep = ""))
  ihmm_res <- fl %>%
    lapply(function(x) {
      cat(x, "\n")
      read_csv(paste(datadir,x,sep = ""), 
               col_names = TRUE,
               col_types = cols(
                 `1` = col_integer(),
                 `2` = col_integer(),
                 `3` = col_integer(),
                 `1_1` = col_double(),
                 `2_1` = col_double(),
                 `3_1` = col_double(),
                 `4` = col_double(),
                 `5` = col_double(),
                 `6` = col_double()
               ))
    }) %>%
    bind_rows()
  colnames(ihmm_res) <- c("Subject", "calcpoint", "N", "mem1", "mem2", "mem3", "mem4", "mem5", "mem6")
  
  cat("Loaded ",length(fl)," files \n", length(fl))
  return(ihmm_res)
}

read.tcsv <- function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

getTestSet <- function(Subjects = "")
{
  library(reshape2)
  datadir <- "../../Data/testset/"
  testset <- data.frame()
  if (is.null(Subjects)) Subjects <- unique(Data.All$Subject)
  for (s in Subjects)
  {
    fl <- list.files(path = datadir, pattern = s)
    testset.subj <- fl %>%
      lapply(function(x) {
        cat(x, "\n")
        y <- t(read.tcsv(paste(datadir,x,sep = ""),header = FALSE))
        return(y)
      })
    testset.subj <- data.frame(t(testset.subj[[1]]))
    colnames(testset.subj) <- c("incl")
    n <- nrow(testset.subj) /10
    cvi <- rep(1:10,n)
    cvi <- cvi[order(cvi)]
    testset.subj$cvi <- cvi
    testset.subj$N <- rep(1:n, 10)
    testset.subj <- spread(testset.subj, key="cvi", value = "incl")
    testset.subj$Subject <- s
    for (i in 1:nrow(testset.subj)) testset.subj[i,"rowsums"] <- sum(testset.subj[i,2:11])
    
    testset <- rbind(testset, testset.subj)
  }
  cat("Loaded ",length(fl)," files \n", length(fl))
  return(testset)
}

loadTestSet <- function(filename = "")
{
  datadir <- "../Data/ASRT/"
  fl <- list.files(path = datadir, pattern = paste("*trainsets_",filename,"*",sep = ""))
  trainset <- fl %>%
    lapply(function(x) {
      cat(x, "\n")
      read_csv(paste(datadir,x,sep = ""), 
               col_names = TRUE)
    }) %>%
    bind_rows()
}

getRSItiming <- function()
{
  return(read_delim("../Data/ASRT_RSI_csoportok.csv", delim = ","))
}

unifyFitResult <- function(Subject = "", fname)
{
  fitres <- Subject %>%
    lapply(., getfitResult) %>%
    rbind_all()
  datadir <- "../Matlab/Output/fitResult/"
  write_csv(fitres, path = paste(datadir, fname, sep = ""), col_names = TRUE)
}

unifyTestSet <- function(Subject = "", fname)
{
  fitres <- Subject %>%
    lapply(., getTestSet) %>%
    rbind_all()
  datadir <- "../Matlab/Output/fitResult/testset/"
  write_csv(fitres, path = paste(datadir, fname, sep = ""), col_names = TRUE)
}

selectBestFits <- function(fp)
{
  fp %>%
    group_by(Subject, model, cvi) %>%
    filter(Ressq == min(Ressq)) %>%
    return(.)
}