### Model definitions

### Row format
Model_ABCDEFG <- function(xN,xP,xn,xB,xS,a,b,c,d,e,f,g)
{
  y = a + b*((xN)^(-c)) + d*(log(xP)) + (e+f*xB)*xn - g*xS
  return(y)
}

### Table format ###
ModelB <- function(xdata,x)
{
  y = x[["a"]] + x[["b"]]*(xdata$N)^(-x[["c"]]) - x[["d"]]*(xdata$logP) + (x[["e"]]+x[["f"]]*xdata$B)*xdata$n - x[["g"]]*xdata$S
  return(y)
}

predictRT <- function(x, BlockLength = 85, BlockInSession = 15, SessionNumber = 3)
{
  nn = rep(1:BlockLength,BlockInSession*SessionNumber)
  S = rep(1:SessionNumber, BlockLength*BlockInSession)
  S = S[order(S)]
  N = 1:(BlockLength*BlockInSession*SessionNumber)
  B = rep(1:BlockInSession,BlockLength)
  B = B[order(B)]
  
  xdata = data.frame(n = nn, N = N, B = B, logP = 0, S = S)
  
  return(ModelB(xdata, x))
}