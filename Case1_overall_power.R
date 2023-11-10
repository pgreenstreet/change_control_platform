#traction example
#overall power looking at treatment 2 being the best

#theta0 is the effect against the control
#kprime is the chosen one
#K star is the one of interest
WrongfirstOD=function(kprime,kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial 
{

  kother=(1:3)[-c(kprime,kstar)]
  muZkprime0_1=theta0[kprime]*sqrt(n)/(sd*sqrt(2))
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  muZkotherkprime_1=(theta0[kother]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  muZkstarkprime_1=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  
  muZkstarkprime_2=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd) #let n be double the size 

  Sigmapn=matrix(c(1,1/2,-1/2,-1/2,-.5*sqrt(.5),
                   1/2,1,1/2,0,.5*sqrt(.5),
                   -1/2,1/2,1,1/2,sqrt(.5),
                   -1/2,0,1/2,1,.5*sqrt(.5),
                   -.5*sqrt(.5),.5*sqrt(.5),sqrt(.5),.5*sqrt(.5),1),nrow = 5,byrow = T)
  ubn=c(Inf,Inf,0,0,Inf)
  lbn=c(u[1],l[1],-Inf,-Inf,u[2])
  meansn=c(muZkprime0_1,muZkstar0_1,muZkstarkprime_1,muZkotherkprime_1,muZkstarkprime_2)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]
  
  return(numerator)
}

WrongfirstND=function(kprime,kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial 
{
  
  kother=(1:3)[-c(kprime,kstar)]
  muZkprime0_1=theta0[kprime]*sqrt(n)/(sd*sqrt(2))
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  muZkotherkprime_1=(theta0[kother]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  muZkstarkprime_1=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  
  muZkstarkprime_2=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd*sqrt(2)) #let n be double the size 
  
  Sigmapn=matrix(c(1,1/2,-1/2,-1/2,
                   1/2,1,1/2,0,
                   -1/2,1/2,1,1/2,
                   -1/2,0,1/2,1),nrow = 4,byrow = T)
  ubn=c(Inf,Inf,0,0)
  lbn=c(u[1],l[1],-Inf,-Inf)
  meansn=c(muZkprime0_1,muZkstar0_1,muZkstarkprime_1,muZkotherkprime_1)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]*pmvnorm(lower=u[2],upper=Inf,mean = muZkstarkprime_2,sigma = 1)[[1]]
  return(numerator)
}


Firstright=function(kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial (kstar now euals kprime)
{
  
  kother=(1:3)[-kstar]
  kother1=kother[1]
  kother2=kother[2]
  
 
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  muZkother1kstar_1=(theta0[kother1]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))
  muZkother2kstar_1=(theta0[kother2]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))
  
  
  Sigmapn=matrix(c(1,-1/2,-1/2,
                   -1/2,1,1/2,
                   -1/2,1/2,1),nrow = 3,byrow = T)
  ubn=c(Inf,0,0)
  lbn=c(u[1],-Inf,-Inf)
  meansn=c(muZkstar0_1,muZkother1kstar_1,muZkother2kstar_1)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]
  return(numerator)
}


Secondright=function(kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial (kstar now euals kprime)
{
  
  kother=(1:3)[-kstar]
  kother1=kother[1]
  kother2=kother[2]
  
  
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  
  muZkother10_1=(theta0[kother1])*sqrt(n)/(sd*sqrt(2))
  muZkother20_1=(theta0[kother2])*sqrt(n)/(sd*sqrt(2))
  
  muZkother1kstar_1=(theta0[kother1]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))
  muZkother2kstar_1=(theta0[kother2]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))
  
  muZkstar0_2=theta0[kstar]*sqrt(n)/(sd)
  
  muZkother1kstar_2=(theta0[kother1]-theta0[kstar])*sqrt(n)/(sd)
  muZkother2kstar_2=(theta0[kother2]-theta0[kstar])*sqrt(n)/(sd)
  
  
  #equation 4 (treatments other1 and other2  carrys on)
  
  Sigmapn4=matrix(c(1, sqrt(1/2), 1/2, -1/2*sqrt(1/2), 1/2, -1/2*sqrt(1/2),
                   sqrt(1/2), 1, 1/2*sqrt(1/2),-1/2, 1/2*sqrt(1/2), -1/2,
                   1/2, 1/2*sqrt(1/2) ,1, 1/2*sqrt(1/2),1/2, 0,
                   -1/2*sqrt(1/2),-1/2, 1/2*sqrt(1/2), 1, 0, 1/2,
                   1/2, 1/2*sqrt(1/2), 1/2, 0, 1, 1/2*sqrt(1/2),
                   -1/2*sqrt(1/2),-1/2, 0, 1/2, 1/2*sqrt(1/2), 1 ),nrow = 6,byrow = T)
  ubn4=c(u[1],Inf,u[1],0,u[1],0)
  lbn4=c(l[1],u[2],l[1],-Inf,l[1],-Inf)
  meansn4=c(muZkstar0_1,muZkstar0_2,muZkother10_1,muZkother1kstar_2,muZkother20_1,muZkother2kstar_2)
  set.seed(1)
  equation4=pmvnorm(lower=lbn4,upper=ubn4,mean = meansn4,corr = Sigmapn4)[[1]]

  
  #equation 1 (both other treatments stop early)
  
  
  Sigmapn1=Sigmapn4[c(1,2,3,5),c(1,2,3,5)]
  ubn1=c(u[1],Inf,l[1],l[1])
  lbn1=c(l[1],u[2],-Inf,-Inf)
  meansn1=meansn4[c(1,2,3,5)]
  set.seed(1)
  equation1=pmvnorm(lower=lbn1,upper=ubn1,mean = meansn1,corr = Sigmapn1)[[1]]
  
  
  #equation 2 (treatments other1 carrys on)
  
  Sigmapn2=Sigmapn4[c(1,2,3,4,5),c(1,2,3,4,5)]
  ubn2=c(u[1],Inf,u[1],0,l[1])
  lbn2=c(l[1],u[2],l[1],-Inf,-Inf)
  meansn2=meansn4[c(1,2,3,4,5)]
  set.seed(1)
  equation2=pmvnorm(lower=lbn2,upper=ubn2,mean = meansn2,corr = Sigmapn2)[[1]]
  
  
  #equation 3 (treatments other2 carrys on)
  
  
  Sigmapn3=Sigmapn4[c(1,2,3,5,6),c(1,2,3,5,6)]
  ubn3=c(u[1],Inf,l[1],u[1],0)
  lbn3=c(l[1],u[2],-Inf,l[1],-Inf)
  meansn3=meansn4[c(1,2,3,5,6)]
  set.seed(1)
  equation3=pmvnorm(lower=lbn3,upper=ubn3,mean = meansn3,corr = Sigmapn3)[[1]]
  
 
  
  result=equation1+equation2+equation3+equation4
  return(result)
}



#Firstright(kstar=2,theta0=c(0.1,0.5,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)

#Secondright(kstar=2,theta0=c(0.1,0.5,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)
#
#WrongfirstOD(kprime = 1,kstar=2,theta0=c(0.1,0.5,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)

#WrongfirstND(kprime = 1,kstar=2,theta0=c(0.1,0.5,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)


#Firstright(kstar=1,theta0=c(0.5,0.1,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)

#Secondright(kstar=1,theta0=c(0.5,0.1,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)

#WrongfirstOD(kprime = 2,kstar=1,theta0=c(0.5,0.1,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)

#WrongfirstND(kprime = 2,kstar=1,theta0=c(0.5,0.1,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)


OverallpowerOandN=function(theta0,sd,u,l,n) #Does it for whichever treatment is the best
{
  kstar=which.max(theta0)
  
  Nkstar=(1:3)[-kstar]
  Part1=Firstright(kstar,theta0,sd,u,l,n)
  
  Part2=Secondright(kstar,theta0,sd,u,l,n)
  
  OPart3a=WrongfirstOD(kprime = Nkstar[1],kstar,theta0,sd,u,l,n)
  
  OPart3b=WrongfirstOD(kprime = Nkstar[2],kstar,theta0,sd,u,l,n)
  
  NPart3a=WrongfirstND(kprime = Nkstar[1],kstar,theta0,sd,u,l,n)
  
  NPart3b=WrongfirstND(kprime = Nkstar[2],kstar,theta0,sd,u,l,n)

  OverallpowerN=Part1+Part2+NPart3a+NPart3b
    
  OverallpowerO=Part1+Part2+OPart3a+OPart3b

  return(list("OverallpowerN"=OverallpowerN,"OverallpowerO"=OverallpowerO))
}

#use overpower2 if needing to find a mistake.
#OverallpowerOandN(theta0=c(-100,5,5),sd=1,u=c(2,2),l=c(-2,2),n=10)
#OverallpowerOandN(theta0=c(0.1,-0.3,0.3),sd=1,u=c(2,2),l=c(-2,2),n=10)



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


AlldataOverpowerfocus <- foreach(x=1:dim(Mtheta0)[1]) %dopar% 
  {
    
    library(mvtnorm)
    OverallpowerOandN=OverallpowerOandN
    
    Firstright=Firstright
    
    Secondright=Secondright
    
    WrongfirstOD=WrongfirstOD
    
    WrongfirstND=WrongfirstND
     
    ttheta0=Mtheta0[x,]
    tn=42
    tl=MAMStouse$l
    tu=MAMStouse$u
    tsd=1
    ans=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)
    Olddatagtheta0=ans$OverallpowerO
    Newdatagtheta0=ans$OverallpowerN
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataOverpowerfocus,file = "AlldataOverpowerfocus.Rdata")


#Just looking at some points
ttheta0=c(0,0.6,0)
tn=42
tl=MAMStouse$l
tu=MAMStouse$u
somepoint1=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)

ttheta0=c(0.55,0.6,0)
tn=42
tl=MAMStouse$l
tu=MAMStouse$u
somepoint2=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)
somepoint1$OverallpowerO-somepoint2$OverallpowerO
