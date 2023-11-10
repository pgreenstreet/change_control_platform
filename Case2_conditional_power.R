#traction example
# need not change anything with 2
# treatment 2 is the one we care about but treatment 1 we take forward first
#conditional power
library(mvtnorm)


OlddataT1T2=function(theta0,sd,u,l,n)
{
  #Parts for each section of the equation
  Parts=rep(NA,5)
  
  muZ10_1=(theta0[1])*sqrt(n)/(sd*sqrt(2))
  muZ20_1=(theta0[2])*sqrt(n)/(sd*sqrt(2))
  muZ12_1=(theta0[1]-theta0[2])*sqrt(n)/(sd*sqrt(2))
  muZ21_1=-muZ12_1
  
  muZ12_2=(theta0[1]-theta0[2])*sqrt(n)/(sd)
  muZ21_2=-muZ12_2
  
  SigmaT1T2=matrix(c(1,1/2,-1/2,-.5*sqrt(.5),
                   1/2,1,1/2,.5*sqrt(.5),
                   -1/2,1/2,1,sqrt(.5),
                   -.5*sqrt(.5),.5*sqrt(.5),sqrt(.5),1),nrow = 4,byrow = T)
  ubn=c(Inf,Inf,0,Inf)
  lbn=c(u[1],l[1],-Inf,u[2])
  meanT1T2=c(muZ10_1,muZ20_1,muZ21_1,muZ21_2)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meanT1T2,corr = SigmaT1T2)[[1]]
  print(numerator)
  
  ubd=c(Inf,Inf,0)
  lbd=c(u[1],l[1],-Inf)

  set.seed(1)
  denominator=pmvnorm(lower=lbd,upper=ubd,mean = meanT1T2[1:3],corr = SigmaT1T2[1:3,1:3])[[1]]
  print(denominator)
  #print(Parts)
  Oldpower=numerator/denominator
  return(Oldpower)
}
#Olddata(ttheta0,tsd,tu,tl,tn)




NewdataT1T2=function(theta0,sd,u,n)
{ 
  mu1=-(theta0[1]-theta0[2])*sqrt(n)/(sd*sqrt(2))
  set.seed(1)
  Total=pmvnorm(lower=c(u[2]),upper=c(Inf),mean = c(mu1),sigma = 1)[[1]]
  return(Total)
}
#Newdata(ttheta0,tsd,tu,tn)






#looking at treatment 3 now assuming treatment 1 is the control
#assuming 1 stage

OldorNewdataS1T3=function(theta0,sd,u,l,n)
{ 
  muS1T3=c((theta0[3]-theta0[1])*sqrt(n)/(sd*sqrt(2)),(theta0[3]-theta0[1])*sqrt(n)/(sd*sqrt(2)))
  SigmaS1T3=matrix(c(1,sqrt(1/2),sqrt(1/2),1),nrow = 2,byrow = T)
  set.seed(1)
  part1=pmvnorm(lower=c(u[1]),upper=c(Inf),mean = muS1T3[1],sigma = SigmaS1T3[1,1])[[1]]
  set.seed(1)
  part2=pmvnorm(lower=c(l[1],u[2]),upper=c(u[1],Inf),mean = muS1T3,corr = SigmaS1T3)[[1]]
  
  Total=part1+part2
  return(Total)
}




#looking at treatment 3 now assuming treatment 1 is the control
#assuming 2 stage

OlddataT3=function(theta0,sd,u,l,n)
{
  #Parts for each section of the equation
  Parts=rep(NA,5)
  
  muZ10_1=(theta0[1])*sqrt(n)/(sd*sqrt(2))
  muZ10_2=(theta0[1])*sqrt(n)/(sd)
  
  muZ30_1=(theta0[3])*sqrt(n)/(sd*sqrt(2))
  muZ31_1=(theta0[3]-theta0[1])*sqrt(n)/(sd*sqrt(2))
  
  muZ20_1=(theta0[2])*sqrt(n)/(sd*sqrt(2))
  
  muZ12_2=(theta0[1]-theta0[2])*sqrt(n)/(sd)
  muZ21_2=-muZ12_2
  
  muZ31_2=(theta0[3]-theta0[1])*sqrt(n)/(sd)
  
  
  
  SigmaT3=matrix(c(1,sqrt(.5),0,0,.5,-.5*sqrt(.5),0,
                   sqrt(.5),1,.5*sqrt(.5),-.5*sqrt(.5),.5*sqrt(.5),-.5,-.25,
                   0,.5*sqrt(.5),1,.5,0,0,.5*sqrt(.5),
                   0,-.5*sqrt(.5),.5,1,0,.5*sqrt(.5),sqrt(.5),
                   .5,.5*sqrt(.5),0,0,1,.5*sqrt(.5),0,
                   -.5*sqrt(.5),-.5,0,.5*sqrt(.5),.5*sqrt(.5),1,.25,
                   0,-.25,.5*sqrt(.5),sqrt(.5),0,.25,1),nrow = 7,byrow = T)
  
  meanT3=c(muZ10_1,muZ10_2,muZ30_1,muZ31_1,muZ20_1,muZ21_2,muZ31_2)
  
  ubnp1=c(u[1],Inf,u[1],l[1],Inf)
  lbnp1=c(l[1],u[2],l[1],-Inf,u[2])
  set.seed(1)
  numeratorp1=pmvnorm(lower=lbnp1,upper=ubnp1,mean = meanT3[c(1,2,3,5,7)],corr = SigmaT3[c(1,2,3,5,7),c(1,2,3,5,7)])[[1]]
  
  
  ubnp2=c(u[1],Inf,u[1],u[1],0,Inf)
  lbnp2=c(l[1],u[2],l[1],l[1],-Inf,u[2])
  set.seed(1)
  numeratorp2=pmvnorm(lower=lbnp2,upper=ubnp2,mean = meanT3[c(1,2,3,5,6,7)],corr = SigmaT3[c(1,2,3,5,6,7),c(1,2,3,5,6,7)])[[1]]
  
  
  ubnp3=c(u[1],Inf,Inf,0,l[1],Inf)
  lbnp3=c(l[1],u[2],u[1],-Inf,-Inf,u[2])
  set.seed(1)
  numeratorp3=pmvnorm(lower=lbnp3,upper=ubnp3,mean = meanT3[c(1,2,3,4,5,7)],corr = SigmaT3[c(1,2,3,4,5,7),c(1,2,3,4,5,7)])[[1]]
  
  ubnp4=c(u[1],Inf,Inf,0,u[1],0,Inf)
  lbnp4=c(l[1],u[2],u[1],-Inf,l[1],-Inf,u[2])
  set.seed(1)
  numeratorp4=pmvnorm(lower=lbnp4,upper=ubnp4,mean = meanT3[c(1,2,3,4,5,6,7)],corr = SigmaT3[c(1,2,3,4,5,6,7),c(1,2,3,4,5,6,7)])[[1]]
  
  numerator=numeratorp1+numeratorp2+numeratorp3+numeratorp4
  #print(numerator)
  
  ubdp1=c(u[1],Inf,u[1],l[1])
  lbdp1=c(l[1],u[2],l[1],-Inf)
  set.seed(1)
  denominatorp1=pmvnorm(lower=lbdp1,upper=ubdp1,mean = meanT3[c(1,2,3,5)],corr = SigmaT3[c(1,2,3,5),c(1,2,3,5)])[[1]]
  
  
  ubdp2=c(u[1],Inf,u[1],u[1],0)
  lbdp2=c(l[1],u[2],l[1],l[1],-Inf)
  set.seed(1)
  denominatorp2=pmvnorm(lower=lbdp2,upper=ubdp2,mean = meanT3[c(1,2,3,5,6)],corr = SigmaT3[c(1,2,3,5,6),c(1,2,3,5,6)])[[1]]
  
  
  ubdp3=c(u[1],Inf,Inf,0,l[1])
  lbdp3=c(l[1],u[2],u[1],-Inf,-Inf)
  set.seed(1)
  denominatorp3=pmvnorm(lower=lbdp3,upper=ubdp3,mean = meanT3[c(1,2,3,4,5)],corr = SigmaT3[c(1,2,3,4,5),c(1,2,3,4,5)])[[1]]
  
  ubdp4=c(u[1],Inf,Inf,0,u[1],0)
  lbdp4=c(l[1],u[2],u[1],-Inf,l[1],-Inf)
  set.seed(1)
  denominatorp4=pmvnorm(lower=lbdp4,upper=ubdp4,mean = meanT3[c(1,2,3,4,5,6)],corr = SigmaT3[c(1,2,3,4,5,6),c(1,2,3,4,5,6)])[[1]]
  
  denominator=denominatorp1+denominatorp2+denominatorp3+denominatorp4
  
 
  Oldpower=numerator/denominator
  return(Oldpower)
}

NewdataT3=function(theta0,sd,u,n)
{ 
  mu1=-(theta0[1]-theta0[3])*sqrt(n)/(sd*sqrt(2))
  set.seed(1)
  Total=pmvnorm(lower=c(u[2]),upper=c(Inf),mean = c(mu1),sigma = 1)[[1]]
  return(Total)
}


# for just treatment 1 vs 2 
library(doParallel)
library(foreach) 
library(parallel)

ntestsT2=100
ntestsT3=3
ntests=100
theta01= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntests),each=ntestsT2),each=ntestsT3)
theta02= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntestsT2),times = ntests),each=ntestsT3)
theta03= rep(rep(seq(-sqrt(2)*qnorm(0.55),sqrt(2)*qnorm(0.55),length.out = ntestsT3),times = ntests),times=ntestsT2)
Mtheta0=cbind(theta01,theta02,theta03)
Mtheta0
dim(Mtheta0)[1]


# Useful commands
load("triboundslater.Rdata")

cl <- makeCluster(4) 

registerDoParallel(cl) 


AlldataConpowerT12lateFocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    library(mvtnorm)
    OlddataT1T2=OlddataT1T2
    NewdataT1T2=NewdataT1T2
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=triboundslater[[3]][1,]
    tu=triboundslater[[2]][1,]
    tsd=1
    
    Olddatagtheta0=OlddataT1T2(ttheta0,tsd,tu,tl,tn)
    Newdatagtheta0=NewdataT1T2(ttheta0,tsd,tu,tn)
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataConpowerT12lateFocus,file = "AlldataConpowerT12lateFocus.Rdata")







#for treatment 3 against 1 
ntestsT2=3
ntests=100
theta01= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntests),each=ntests),each=ntestsT2)
theta02= rep(rep(seq(-sqrt(2)*qnorm(0.55),sqrt(2)*qnorm(0.55),length.out=ntestsT2),times = ntests),each=ntests)
theta03= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out = ntests),times = ntests),times=ntestsT2)
Mtheta0=cbind(theta01,theta02,theta03)
Mtheta0
dim(Mtheta0)[1]




cl <- makeCluster(4) 

registerDoParallel(cl) 


AlldataConpowerT3S1lateFocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    library(mvtnorm)
    OldorNewdataS1T3=OldorNewdataS1T3
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=triboundslater[[3]][1,]
    tu=triboundslater[[2]][1,]
    tsd=1
    
    Olddatagtheta0=Newdatagtheta0=OldorNewdataS1T3(ttheta0,tsd,tu,tl,tn)
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataConpowerT3S1lateFocus,file = "AlldataConpowerT3S1lateFocus.Rdata")

cl <- makeCluster(4) 

registerDoParallel(cl) 


AlldataConpowerT3S2lateFocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    library(mvtnorm)
    OlddataT3=OlddataT3
    NewdataT3=NewdataT3
    ttheta0=Mtheta0[x,]
    tn=42
    tl=triboundslater[[3]][1,]
    tu=triboundslater[[2]][1,]
    tsd=1
    Olddatagtheta0=OlddataT3(ttheta0,tsd,tu,tl,tn)
    Newdatagtheta0=NewdataT3(ttheta0,tsd,tu,tn)
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataConpowerT3S2lateFocus,file = "AlldataConpowerT3S2lateFocus.Rdata")



#still to look at 

library(doParallel)
library(foreach) 
library(parallel)
# Useful commands
load("infboundslater.Rdata")


ntestsT2=100
ntestsT3=3
ntests=100
theta01= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntests),each=ntestsT2),each=ntestsT3)
theta02= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntestsT2),times = ntests),each=ntestsT3)
theta03= rep(rep(seq(-sqrt(2)*qnorm(0.55),sqrt(2)*qnorm(0.55),length.out = ntestsT3),times = ntests),times=ntestsT2)
Mtheta0=cbind(theta01,theta02,theta03)
Mtheta0
dim(Mtheta0)[1]


cl <- makeCluster(4) 

registerDoParallel(cl) 


AlldataConpowerT12lateInfFocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    library(mvtnorm)
    OlddataT1T2=OlddataT1T2
    NewdataT1T2=NewdataT1T2
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=infboundslater[[3]][1,]
    tu=infboundslater[[2]][1,]
    tsd=1
    
    Olddatagtheta0=OlddataT1T2(ttheta0,tsd,tu,tl,tn)
    Newdatagtheta0=NewdataT1T2(ttheta0,tsd,tu,tn)
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataConpowerT12lateInfFocus,file = "AlldataConpowerT12lateInfFocus.Rdata")

ntestsT2=3
ntests=100
theta01= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out=ntests),each=ntests),each=ntestsT2)
theta02= rep(rep(seq(-sqrt(2)*qnorm(0.55),sqrt(2)*qnorm(0.55),length.out=ntestsT2),times = ntests),each=ntests)
theta03= rep(rep(seq(-sqrt(2)*qnorm(0.65),2*sqrt(2)*qnorm(0.65),length.out = ntests),times = ntests),times=ntestsT2)
Mtheta0=cbind(theta01,theta02,theta03)
Mtheta0
dim(Mtheta0)[1]

cl <- makeCluster(4) 

registerDoParallel(cl) 


AlldataConpowerT3S1lateInfFocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    library(mvtnorm)
    OldorNewdataS1T3=OldorNewdataS1T3
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=infboundslater[[3]][1,]
    tu=infboundslater[[2]][1,]
    tsd=1
    
    Olddatagtheta0=Newdatagtheta0=OldorNewdataS1T3(ttheta0,tsd,tu,tl,tn)
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataConpowerT3S1lateInfFocus,file = "AlldataConpowerT3S1lateInfFocus.Rdata")

cl <- makeCluster(4) 

registerDoParallel(cl) 


AlldataConpowerT3S2lateInfFocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    library(mvtnorm)
    OlddataT3=OlddataT3
    NewdataT3=NewdataT3
    ttheta0=Mtheta0[x,]
    tn=42
    tl=infboundslater[[3]][1,]
    tu=infboundslater[[2]][1,]
    tsd=1
    Olddatagtheta0=OlddataT3(ttheta0,tsd,tu,tl,tn)
    Newdatagtheta0=NewdataT3(ttheta0,tsd,tu,tn)
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataConpowerT3S2lateInfFocus,file = "AlldataConpowerT3S2lateInfFocus.Rdata")
