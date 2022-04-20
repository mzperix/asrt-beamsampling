### RESIDUALS ###

calc.res <- function(Data, fp)
  # Calculate residuals on Data for given fitted model parameters
{
  dat <- Data %>%
    group_by(Subject) %>%
    mutate(N = 1:n(), B = Block, S = Session) %>%
    group_by(Session, Block) %>%
    mutate(n = 1:n()) %>%
    ungroup()
  dat[dat$TT == "H", "logP"] = log(0.625)
  dat[dat$TT == "L", "logP"] = log(0.125)
  
  dat.res <- data.frame()
  for (r in 1:nrow(fp))
  {
    param <- fp[r,]
    dat.res <- transmute(dat, 
                         N = N, 
                         Subject = Subject,
                         model = param$model, 
                         predRT = ModelB(xdata = dat, x = param[,c("a","b","c","d","e","f","g")]),
                         resinvRT = 1/predRT-1/firstRT, 
                         resRT = predRT-firstRT,
                         cvi = as.factor(param$cvi)) %>%
      rbind(dat.res,.)
  }
  return(dat.res)
}

calc.res.vector <- function(xdata, x)
{
  return(xdata$firstRT-ModelB(xdata = xdata, x = x[,c("a","b","c","d","e","f","g")]))
}

# Data.res <- calc.res(Data %>% filter(Subject == 101), fp %>% filter(Subject == 101))
# Data.res %>% ggplot(aes(x = N, y = resinvRT, color = bootID)) + geom_point(size = 1, alpha = 0.2)

calc.res.cv <- function(Data.All, fp, testset, Model)
{
  ### Calculates residuals for each cross-validation testset
  ### Only for firstACC == 1 and TT  in c(H,L)
  
  fp.model <- subset(fp, model == Model)
  
  Data <- Data.All %>%
    group_by(Subject) %>%
    mutate(N = 1:n()) %>%
    subset(., firstACC == 1) %>%
    subset(., TT %in% c("H", "L")) %>%
    mutate(N2 = 1:n())
  
  df.cvi <- testset %>%
    gather(., key = "cvi", value = "train", -N, -Subject, -rowsums) %>%
    subset(., train != 1) %>%
    .[order(.$Subject, .$N),]

  Data <- merge(Data, df.cvi, by.x = c("Subject", "N2"), by.y = c("Subject", "N")) %>%
    .[order(.$Subject, .$N2),] %>%
    ungroup() %>%
    group_by(Subject) %>%
    mutate(B = Block, S = Session, n = (N-1) %% 85 + 1,
           cvi = as.integer(cvi))
    
  Data[Data$TT == "H", "logP"] = log(0.625)
  Data[Data$TT == "L", "logP"] = log(0.125)
  
  if(Model != "1000000")
  {
    Data <- Data %>%
      #subset(., Subject %in% head(unique(Data$Subject), 1)) %>%
      group_by(Subject, cvi) %>%
      mutate(predRT = ModelB(xdata = data.frame(N = N, B = B, S = S, logP = logP, n = n),
                             x = fp.model[fp.model$Subject == head(Subject,1) & 
                                            fp.model$cvi == head(cvi,1),]))
  } else 
  {  
    cv_means <- Data %>%
      group_by(Subject, cvi) %>%
      summarise(meanRT = mean(firstRT), 
                N = n(),
                sumRT = sum(firstRT))
    
    Data <- Data %>%
      group_by(Subject, cvi) %>%
      mutate(predRT = sum(cv_means[cv_means$Subject == head(Subject,1) &
                                      cv_means$cvi == head(cvi,1), "sumRT"]) /
               sum(cv_means[cv_means$Subject == head(Subject,1) &
                              cv_means$cvi == head(cvi,1) , "N"]))
  }
  
  Data$resRT <- Data$firstRT - Data$predRT
  
  return(Data)
}