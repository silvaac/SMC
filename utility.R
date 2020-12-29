require(tidyquant)
require(tidyverse)

cfun2 <- function(a,A){
  cf = mean(abs(a-A))
  return(cf)
}

pfun <- function(R){
  Rs = R
  ts = sample(1:length(R),2)
  Rs[ts[c(2,1)]] = R[ts]
  return(Rs)
}

efun6 <- function(a,b,c,b2,R,nLag){
  A  = acf(R,lag.max = nLag,plot=FALSE)
  B  = acf(abs(R),lag.max = 5*nLag,plot=FALSE)
  C  = ccf(abs(R),R,lag.max = nLag,plot=FALSE)
  B2 = acf(R*R,lag.max = 5*nLag,plot=FALSE)
  Z=cfun2(a,A$acf[-1])+cfun2(b,B$acf[-1])+cfun2(c,C$acf)
  Z = Z+cfun2(b2,B2$acf[-1])
  return(Z)  
}
