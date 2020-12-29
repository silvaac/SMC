source("utility.R")
mySeed = 37
fSave  = paste0("smcSPX_",mySeed)

temp = 1e-4
maxT = 20000
maxG = 2000
maxS = 200
cgoal= 5e-2  # 1e-3
wc   = 0.9
afac = 0.5
maxA = 120e6

#####
# Real mkt data: spy  
sp = readRDS("spxData")

######
# RV with empirical distribution
set.seed(mySeed)
a  = cdfR::cdfR.default(sp$R)
rs = cdfR::rCdf(a$cdf,length(sp$R))
rs[is.na(rs)] = 0

#######
nLag = 50
# Target stat & initialize
aR = acf(sp$R,lag.max = nLag)
aR = aR$acf[-1]
bR = acf(abs(sp$R),lag.max = 5*nLag)
bR = bR$acf[-1]
cR = ccf(abs(sp$R),sp$R,lag.max = nLag)
cR = cR$acf
bR2 = acf(sp$R^2,lag.max = 5*nLag)
bR2 = bR2$acf[-1]

# Error function: first guess
cmin = efun6(aR,bR,cR,bR2,rs,nLag = nLag)

# Loop
cost = cmin
foundMin = FALSE
nTot  = 0
nGood = 0
nAll  = 0
ra    = rs
tini  = temp
E     = matrix(NA,round(maxA/1e6),1)

while(nAll<maxA){
  nTot = nTot+1
  nAll = nAll+1
  if(nAll%%100000 == 0){
    print(nAll/1e6)
    E[round(nAll/1e6)] = tE
  }
  if(nAll%%500000 == 0) saveRDS(ra,file=fSave)
  
  cmax = cost-temp*log(runif(1))
  trs  = pfun(rs)
  tE   = efun6(aR,bR,cR,bR2,trs,nLag)
  
  if(tE<=cmax){
    nGood = nGood+1
    rs    = trs
    cost  = tE
  }
  if(cost<=cgoal){
    ra = rs
    foundMin = TRUE
    break
  }
  if(nTot<maxT && nGood<maxG){
    if(cost<(cmin*wc)){
      cmin = cost
      ra   = rs
    }
    next
  }
  if(temp==tini && nTot>1.5*nGood){
    tini = 10*temp
    temp = tini
  }else{
    if(nGood<=maxS){
      afac=sqrt(afac)
      maxT=maxT*sqrt(2)
      temp=tini
    }else{         
      temp=temp*afac
    }
  }
  nTot  = 0
  nGood = 0
}
saveRDS(ra,file=fSave)

A  = acf(ra,lag.max = nLag)
B  = acf(abs(ra),lag.max = 5*nLag)
C  = ccf(abs(ra),ra,lag.max = nLag)
D  = acf(ra^2,lag.max = 5*nLag)

# Quality of algo
plot(aR);points(A$acf[-1],col="blue")
plot(bR,ylim = c(-0.5,0.5));points(B$acf[-1],col="blue")
plot(cR,ylim = c(-0.5,0.5));points(C$acf,col="blue")
plot(bR2,ylim = c(-0.5,0.5));points(D$acf[-1],col="blue")
