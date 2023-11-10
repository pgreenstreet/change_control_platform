#traction example
# treatment 2 is the one we care about but treatmtn 1 we take forward first
#conditonal power
library(mvtnorm)


Olddata=function(theta0,sd,u,l,n)
{
  #Parts for each section of the equation
  Parts=rep(NA,5)
  
  muZ10_1=(theta0[1])*sqrt(n)/(sd*sqrt(2))
  muZ20_1=(theta0[2])*sqrt(n)/(sd*sqrt(2))
  muZ12_1=(theta0[1]-theta0[2])*sqrt(n)/(sd*sqrt(2))
  muZ13_1=(theta0[1]-theta0[3])*sqrt(n)/(sd*sqrt(2))
  muZ21_1=-muZ12_1
  muZ31_1=-muZ13_1
  
  muZ12_2=(theta0[1]-theta0[2])*sqrt(n)/(sd)
  muZ21_2=-muZ12_2
  
  Sigmapn=matrix(c(1,1/2,-1/2,-1/2,-.5*sqrt(.5),
                   1/2,1,1/2,0,.5*sqrt(.5),
                   -1/2,1/2,1,1/2,sqrt(.5),
                   -1/2,0,1/2,1,.5*sqrt(.5),
                   -.5*sqrt(.5),.5*sqrt(.5),sqrt(.5),.5*sqrt(.5),1),nrow = 5,byrow = T)
  ubn=c(Inf,Inf,0,0,Inf)
  lbn=c(u[1],l[1],-Inf,-Inf,u[2])
  meansn=c(muZ10_1,muZ20_1,muZ21_1,muZ31_1,muZ21_2)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]
  print(numerator)
  Sigmapd=matrix(c(1,1/2,-1/2,-1/2,
                   1/2,1,1/2,0,
                   -1/2,1/2,1,1/2,
                   -1/2,0,1/2,1),nrow = 4,byrow = T)
  ubd=c(Inf,Inf,0,0)
  lbd=c(u[1],l[1],-Inf,-Inf)
  meansd=c(muZ10_1,muZ20_1,muZ21_1,muZ31_1)
  set.seed(1)
  denominator=pmvnorm(lower=lbd,upper=ubd,mean = meansd,corr = Sigmapd)[[1]]
  print(denominator)
  #print(Parts)
  Oldpower=numerator/denominator
  return(Oldpower)
}
#Olddata(ttheta0,tsd,tu,tl,tn)


Newdata=function(theta0,sd,u,n)
{ 
  mu1=-(theta0[1]-theta0[2])*sqrt(n)/(sd*sqrt(2))
  set.seed(1)
  Total=pmvnorm(lower=c(u[2]),upper=c(Inf),mean = c(mu1),sigma = 1)[[1]]
  return(Total)
}
#Newdata(ttheta0,tsd,tu,tn)


#ttheta0=c(.4,1,-1000)
#tn=14
#tl=c(-0.5,1.5)
#tu=c(2,1.5)
#tsd=1
#Olddata(ttheta0,tsd,tu,tl,tn)
#Newdata(ttheta0,tsd,tu,tn)

ntestsT3=3
ntests=100
theta01= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntests),each=ntests),each=ntestsT3)
theta02= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntests),times = ntests),each=ntestsT3)
theta03= rep(rep(seq(-sqrt(2)*qnorm(0.55),sqrt(2)*qnorm(0.55),length.out = ntestsT3),times = ntests),times=ntests)
Mtheta0=cbind(theta01,theta02,theta03)
Mtheta0
dim(Mtheta0)[1]

#thinking do 40 cores



library(doParallel)
library(foreach) 
library(parallel)
# Useful commands
load("MAMStouse.Rdata")

cl <- makeCluster(4) 
registerDoParallel(cl) 


AlldataConpowerfocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% 
{
  
  library(mvtnorm)
  Olddata=Olddata
  Newdata=Newdata  
  ttheta0=Mtheta0[x,]
  tn=42
  tl=MAMStouse$l
  tu=MAMStouse$u
  tsd=1
  Olddatagtheta0=Olddata(ttheta0,tsd,tu,tl,tn)
  Newdatagtheta0=Newdata(ttheta0,tsd,tu,tn)
  return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
}


stopCluster(cl) 

save(AlldataConpowerfocus,file = "AlldataConpowerfocus.Rdata")
