#traction example
library(mvtnorm)
#theta0 is the effect against the control
#kprime is the chosen one
#K star is the one of interest

#treatment 1 or 2 is the best
WrongfirstODT12=function(kprime,kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial 
{

  kother=(1:3)[-c(kprime,kstar)]
  muZkprime0_1=theta0[kprime]*sqrt(n)/(sd*sqrt(2))
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  muZkstarkprime_1=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  muZkstarkprime_2=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd)

  Sigmapn=matrix(c(1,1/2,-1/2,-.5*sqrt(.5),
                   1/2,1,1/2,.5*sqrt(.5),
                   -1/2,1/2,1,sqrt(.5),
                   -.5*sqrt(.5),.5*sqrt(.5),sqrt(.5),1),nrow = 4,byrow = T)
  ubn=c(Inf,Inf,0,Inf)
  lbn=c(u[1],l[1],-Inf,u[2])
  meansn=c(muZkprime0_1,muZkstar0_1,muZkstarkprime_1,muZkstarkprime_2)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]
  return(numerator)
}

WrongfirstNDT12=function(kprime,kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial 
{
  
  kother=(1:3)[-c(kprime,kstar)]
  muZkprime0_1=theta0[kprime]*sqrt(n)/(sd*sqrt(2))
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  muZkstarkprime_1=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  muZstarkstarkprime_2=(theta0[kstar]-theta0[kprime])*sqrt(n)/(sd*sqrt(2)) 
  
  Sigmapn=matrix(c(1,1/2,-1/2,
                   1/2,1,1/2,
                   -1/2,1/2,1),nrow = 3,byrow = T)
  ubn=c(Inf,Inf,0)
  lbn=c(u[1],l[1],-Inf)
  meansn=c(muZkprime0_1,muZkstar0_1,muZkstarkprime_1)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]*pmvnorm(lower=u[2],upper=Inf,mean = muZstarkstarkprime_2,sigma = 1)[[1]]
  return(numerator)
}


FirstrightT12=function(kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial (kstar now euals kprime)
{
  
  kother=(1:2)[-kstar]
  kother1=kother
  
 
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  muZkother1kstar_1=(theta0[kother1]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))

  
  Sigmapn=matrix(c(1,-1/2,
                   -1/2,1),nrow = 2,byrow = T)
  ubn=c(Inf,0)
  lbn=c(u[1],-Inf)
  meansn=c(muZkstar0_1,muZkother1kstar_1)
  set.seed(1)
  numerator=pmvnorm(lower=lbn,upper=ubn,mean = meansn,corr = Sigmapn)[[1]]
  return(numerator)
}


SecondrightT12=function(kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial (kstar now euals kprime)
{
  
  kother=(1:2)[-kstar]
  kother1=kother[1]
  
  k3=3
  
  
  muZkstar0_1=theta0[kstar]*sqrt(n)/(sd*sqrt(2))
  
  muZkother10_1=(theta0[kother1])*sqrt(n)/(sd*sqrt(2))
  muZk30_1=(theta0[3])*sqrt(n)/(sd*sqrt(2))
  
  muZkother1kstar_1=(theta0[kother1]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))
  muZk3kstar_1=(theta0[3]-theta0[kstar])*sqrt(n)/(sd*sqrt(2))
  
  muZkstar0_2=theta0[kstar]*sqrt(n)/(sd)
  
  muZkother1kstar_2=(theta0[kother1]-theta0[kstar])*sqrt(n)/(sd)
  
  
  
  #equation 4 (treatments other1 and other2  carrys on)
  
  Sigmapn4=matrix(c(1, sqrt(1/2), 1/2, -1/2*sqrt(1/2), 0,0,
                   sqrt(1/2), 1, 1/2*sqrt(1/2),-1/2, 1/2*sqrt(1/2), -1/2*sqrt(1/2),
                   1/2, 1/2*sqrt(1/2) ,1, 1/2*sqrt(1/2),1/2, 0,
                   -1/2*sqrt(1/2),-1/2, 1/2*sqrt(1/2), 1, 0, 1/2*sqrt(1/2),
                   0, 1/2*sqrt(1/2), 0, 0, 1, 1/2,
                   0,1/2*sqrt(1/2), 0, 1/2*sqrt(1/2), 1/2, 1 ),nrow = 6,byrow = T)
  
  ubn4=c(u[1],Inf,u[1],0,Inf,0)
  lbn4=c(l[1],u[2],l[1],-Inf,u[1],-Inf)
  meansn4=c(muZkstar0_1,muZkstar0_2,muZkother10_1,muZkother1kstar_2,muZk30_1,muZk3kstar_1)
  set.seed(1)
  equation4=pmvnorm(lower=lbn4,upper=ubn4,mean = meansn4,corr = Sigmapn4)[[1]]

  
  #equation 1 (both other treatments stop early)
  
  
  Sigmapn1=Sigmapn4[c(1,2,3,5),c(1,2,3,5)]
  ubn1=c(u[1],Inf,l[1],u[1])
  lbn1=c(l[1],u[2],-Inf,-Inf)
  meansn1=meansn4[c(1,2,3,5)]
  set.seed(1)
  equation1=pmvnorm(lower=lbn1,upper=ubn1,mean = meansn1,corr = Sigmapn1)[[1]]
  
  
  #equation 2 (treatments other1 carrys on)
  
  Sigmapn2=Sigmapn4[c(1,2,3,4,5),c(1,2,3,4,5)]
  ubn2=c(u[1],Inf,u[1],0,u[1])
  lbn2=c(l[1],u[2],l[1],-Inf,-Inf)
  meansn2=meansn4[c(1,2,3,4,5)]
  set.seed(1)
  equation2=pmvnorm(lower=lbn2,upper=ubn2,mean = meansn2,corr = Sigmapn2)[[1]]
  
  
  #equation 3 (treatments other2 carrys on)
  
  
  Sigmapn3=Sigmapn4[c(1,2,3,5,6),c(1,2,3,5,6)]
  ubn3=c(u[1],Inf,l[1],Inf,0)
  lbn3=c(l[1],u[2],-Inf,u[1],-Inf)
  meansn3=meansn4[c(1,2,3,5,6)]
  set.seed(1)
  equation3=pmvnorm(lower=lbn3,upper=ubn3,mean = meansn3,corr = Sigmapn3)[[1]]
  
  
  result=equation1+equation2+equation3+equation4
  return(result)
}


#looking at treatment 3
#stage 1
WrongfirstOaNDT3S2=function(kprime,kstar,theta0,sd,u,l,n)
{
  #again basically assuming 1 is the one taken forward but will not effect the code
  kother=(1:2)[-kprime]
  kother1=kother
  muE1_3=c((theta0[kprime])*sqrt(n)/(sd*sqrt(2)),(theta0[kother1]-theta0[kprime])*sqrt(n)/(sd*sqrt(2)))
  SigmaE1_3=matrix(c(1,-.5,-.5,1),nrow = 2,byrow = T)
  #print(muE1_3)
  #print(SigmaE1_3)
  #print(c(u[1],-Inf))
  #print(c(Inf,0))
  
  TotalE1_3=pmvnorm(lower=c(u[1],-Inf),upper=c(Inf,0),mean = muE1_3,corr = SigmaE1_3)[[1]]
  
  #part E_4
  muS1T3=c((theta0[3]-theta0[kprime])*sqrt(n)/(sd*sqrt(2)),(theta0[3]-theta0[kprime])*sqrt(n)/(sd))
  SigmaS1T3=matrix(c(1,sqrt(1/2),sqrt(1/2),1),nrow = 2,byrow = T)
  set.seed(1)
  part1=pmvnorm(lower=c(u[1]),upper=c(Inf),mean = muS1T3[1],sigma = SigmaS1T3[1,1])[[1]]
  set.seed(1)
  part2=pmvnorm(lower=c(l[1],u[2]),upper=c(u[1],Inf),mean = muS1T3,corr = SigmaS1T3)[[1]]
  
  Total=(part1+part2)*TotalE1_3
  return(Total)
  
}



#second stage
WrongfirstODT3S2=function(kprime,kstar,theta0,sd,u,l,n)
{
  kother=(1:3)[-c(kprime,kstar)]
  
  #Parts for each section of the equation
  Parts=rep(NA,5)
  
  muZ10_1=(theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  muZ10_2=(theta0[kprime])*sqrt(n)/(sd)
  
  muZ30_1=(theta0[3])*sqrt(n)/(sd*sqrt(2))
  muZ31_1=(theta0[3]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  
  muZ20_1=(theta0[kother])*sqrt(n)/(sd*sqrt(2))
  
  muZ12_2=(theta0[kprime]-theta0[kother])*sqrt(n)/(sd)
  muZ21_2=-muZ12_2
  
  muZ31_2=(theta0[3]-theta0[kprime])*sqrt(n)/(sd)
  
  #basically assuming that treatment 1 is taken forward to start
  
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
  
  return(numerator)
}

WrongfirstNDT3S2=function(kprime,kstar,theta0,sd,u,l,n)
{
  kother=(1:3)[-c(kprime,kstar)]
  
  #Parts for each section of the equation
  Parts=rep(NA,5)
  
  muZ10_1=(theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  muZ10_2=(theta0[kprime])*sqrt(n)/(sd)
  
  muZ30_1=(theta0[3])*sqrt(n)/(sd*sqrt(2))
  muZ31_1=(theta0[3]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  
  muZ20_1=(theta0[kother])*sqrt(n)/(sd*sqrt(2))
  
  muZ12_2=(theta0[kprime]-theta0[kother])*sqrt(n)/(sd)
  muZ21_2=-muZ12_2
  
  muZstar31_2=(theta0[3]-theta0[kprime])*sqrt(n)/(sd*sqrt(2))
  
  #basically assuming that treatment 1 is taken forward to start
  
  SigmaT3=matrix(c(1,sqrt(.5),0,0,.5,-.5*sqrt(.5),
                   sqrt(.5),1,.5*sqrt(.5),-.5*sqrt(.5),.5*sqrt(.5),-.5,
                   0,.5*sqrt(.5),1,.5,0,0,
                   0,-.5*sqrt(.5),.5,1,0,.5*sqrt(.5),
                   .5,.5*sqrt(.5),0,0,1,.5*sqrt(.5),
                   -.5*sqrt(.5),-.5,0,.5*sqrt(.5),.5*sqrt(.5),1),nrow = 6,byrow = T)
  
  meanT3=c(muZ10_1,muZ10_2,muZ30_1,muZ31_1,muZ20_1,muZ21_2)
  
  ubnp1=c(u[1],Inf,u[1],l[1])
  lbnp1=c(l[1],u[2],l[1],-Inf)
  set.seed(1)
  numeratorp1=pmvnorm(lower=lbnp1,upper=ubnp1,mean = meanT3[c(1,2,3,5)],corr = SigmaT3[c(1,2,3,5),c(1,2,3,5)])[[1]]
  
  
  ubnp2=c(u[1],Inf,u[1],u[1],0)
  lbnp2=c(l[1],u[2],l[1],l[1],-Inf)
  set.seed(1)
  numeratorp2=pmvnorm(lower=lbnp2,upper=ubnp2,mean = meanT3[c(1,2,3,5,6)],corr = SigmaT3[c(1,2,3,5,6),c(1,2,3,5,6)])[[1]]
  
  
  ubnp3=c(u[1],Inf,Inf,0,l[1])
  lbnp3=c(l[1],u[2],u[1],-Inf,-Inf)
  set.seed(1)
  numeratorp3=pmvnorm(lower=lbnp3,upper=ubnp3,mean = meanT3[c(1,2,3,4,5)],corr = SigmaT3[c(1,2,3,4,5),c(1,2,3,4,5)])[[1]]
  
  ubnp4=c(u[1],Inf,Inf,0,u[1],0)
  lbnp4=c(l[1],u[2],u[1],-Inf,l[1],-Inf)
  set.seed(1)
  numeratorp4=pmvnorm(lower=lbnp4,upper=ubnp4,mean = meanT3[c(1,2,3,4,5,6)],corr = SigmaT3[c(1,2,3,4,5,6),c(1,2,3,4,5,6)])[[1]]
  
  numerator=(numeratorp1+numeratorp2+numeratorp3+numeratorp4)*pmvnorm(lower=u[2],upper=Inf,mean=muZstar31_2,sigma=1)[[1]]
  
  return(numerator)
}


#stage 1 
FirstrightT3=function(kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial (kstar now euals kprime)
{
  muZ10_1=(theta0[1])*sqrt(n)/(sd*sqrt(2))
  muZ10_2=(theta0[1])*sqrt(n)/(sd)
  
  muZ20_1=(theta0[2])*sqrt(n)/(sd*sqrt(2))
  muZ20_2=(theta0[2])*sqrt(n)/(sd)
  
  muZ30_1=(theta0[3])*sqrt(n)/(sd*sqrt(2))
  
  muZ13_2=(theta0[1]-theta0[3])*sqrt(n)/(sd*sqrt(2))
  muZ23_2=(theta0[2]-theta0[3])*sqrt(n)/(sd*sqrt(2))
  

  #complete mu and sigma (as the ones for equation 9)
  meansn9=c(muZ30_1,muZ10_1,muZ10_2,muZ13_2,muZ20_1,muZ20_2,muZ23_2)
  
  Sigmapn9=matrix(c(1, 0, 1/2*sqrt(1/2), -1/2, 0,1/2*sqrt(1/2), -1/2,
                    0, 1, sqrt(1/2), 0, 1/2,1/2*sqrt(1/2), 0,
                    1/2*sqrt(1/2), sqrt(1/2), 1, 1/2*sqrt(1/2), 1/2*sqrt(1/2),1/2, 0,
                    -1/2, 0, 1/2*sqrt(1/2), 1, 0,0, 1/2,
                    0, 1/2, 1/2*sqrt(1/2), 0, 1,sqrt(1/2), 0,
                    1/2*sqrt(1/2), 1/2*sqrt(1/2), 1/2, 0, sqrt(1/2),1, 1/2*sqrt(1/2),
                    -1/2, 0, 0, 1/2, 0,1/2*sqrt(1/2), 1),nrow = 7,byrow = T)
  
  
  #equation 1 
  
  Sigmapn1=Sigmapn9[c(1,2,5),c(1,2,5)]
  ubn1=c(Inf,l[1],l[1])
  lbn1=c(u[1],-Inf,-Inf)
  meansn1=meansn9[c(1,2,5)]
  set.seed(1)
  equation1=pmvnorm(lower=lbn1,upper=ubn1,mean = meansn1,corr = Sigmapn1)[[1]]
  
  
  #equation 2 
  
  Sigmapn2=Sigmapn9[c(1,2,5,6),c(1,2,5,6)]
  ubn2=c(Inf,l[1],u[1],u[2])
  lbn2=c(u[1],-Inf,l[1],-Inf)
  meansn2=meansn9[c(1,2,5,6)]
  set.seed(1)
  equation2=pmvnorm(lower=lbn2,upper=ubn2,mean = meansn2,corr = Sigmapn2)[[1]]
  
  #equation 3 
  
  Sigmapn3=Sigmapn9[c(1,2,5,6,7),c(1,2,5,6,7)]
  ubn3=c(Inf,l[1],u[1],Inf,0)
  lbn3=c(u[1],-Inf,l[1],u[2],-Inf)
  meansn3=meansn9[c(1,2,5,6,7)]
  set.seed(1)
  equation3=pmvnorm(lower=lbn3,upper=ubn3,mean = meansn3,corr = Sigmapn3)[[1]]
  
  #equation 4 
  
  Sigmapn4=Sigmapn9[c(1,2,3,5),c(1,2,3,5)]
  ubn4=c(Inf,u[1],u[2],l[1])
  lbn4=c(u[1],l[1],-Inf,-Inf)
  meansn4=meansn9[c(1,2,3,5)]
  set.seed(1)
  equation4=pmvnorm(lower=lbn4,upper=ubn4,mean = meansn4,corr = Sigmapn4)[[1]]
  
  #equation 5 
  
  Sigmapn5=Sigmapn9[c(1,2,3,5,6),c(1,2,3,5,6)]
  ubn5=c(Inf,u[1],u[2],u[1],u[2])
  lbn5=c(u[1],l[1],-Inf,l[1],-Inf)
  meansn5=meansn9[c(1,2,3,5,6)]
  set.seed(1)
  equation5=pmvnorm(lower=lbn5,upper=ubn5,mean = meansn5,corr = Sigmapn5)[[1]]
  
  #equation 6 
  
  Sigmapn6=Sigmapn9[c(1,2,3,5,6,7),c(1,2,3,5,6,7)]
  ubn6=c(Inf,u[1],u[2],u[1],Inf,0)
  lbn6=c(u[1],l[1],-Inf,l[1],u[2],-Inf)
  meansn6=meansn9[c(1,2,3,5,6,7)]
  set.seed(1)
  equation6=pmvnorm(lower=lbn6,upper=ubn6,mean = meansn6,corr = Sigmapn6)[[1]]
  
  #equation 7 
  
  Sigmapn7=Sigmapn9[c(1,2,3,4,5),c(1,2,3,4,5)]
  ubn7=c(Inf,u[1],Inf,0,l[1])
  lbn7=c(u[1],l[1],u[2],-Inf,-Inf)
  meansn7=meansn9[c(1,2,3,4,5)]
  set.seed(1)
  equation7=pmvnorm(lower=lbn7,upper=ubn7,mean = meansn7,corr = Sigmapn7)[[1]]
  
  #equation 8
  
  Sigmapn8=Sigmapn9[c(1,2,3,4,5,6),c(1,2,3,4,5,6)]
  ubn8=c(Inf,u[1],Inf,0,u[1],u[2])
  lbn8=c(u[1],l[1],u[2],-Inf,l[1],-Inf)
  meansn8=meansn9[c(1,2,3,4,5,6)]
  set.seed(1)
  equation8=pmvnorm(lower=lbn8,upper=ubn8,mean = meansn8,corr = Sigmapn8)[[1]]
  
  #equation 9
  
  
  ubn9=c(Inf,u[1],Inf,0,u[1],Inf,0)
  lbn9=c(u[1],l[1],u[2],-Inf,l[1],u[2],-Inf)
  set.seed(1)
  equation9=pmvnorm(lower=lbn9,upper=ubn9,mean = meansn9,corr = Sigmapn9)[[1]]
  
  result=equation1+equation2+equation3+equation4+equation5+equation6+equation7+equation8+equation9
  return(result)
}

#stage 2 
SecondrightT3=function(kstar,theta0,sd,u,l,n) #all looking at a 2 stage trial (kstar now euals kprime)
{
  muZ10_1=(theta0[1])*sqrt(n)/(sd*sqrt(2))
  muZ10_2=(theta0[1])*sqrt(n)/(sd)
  
  muZ20_1=(theta0[2])*sqrt(n)/(sd*sqrt(2))
  muZ20_2=(theta0[2])*sqrt(n)/(sd)
  
  muZ30_1=(theta0[3])*sqrt(n)/(sd*sqrt(2))
  muZ30_2=(theta0[3])*sqrt(n)/(sd)
  
  
  #complete mu and sigma (as the ones for equation 9)
  meansn4=c(muZ30_1,muZ30_2,muZ10_1,muZ10_2,muZ20_1,muZ20_2)
  
  Sigmapn4=matrix(c(1, sqrt(1/2), 0, 1/2*sqrt(1/2), 0, 1/2*sqrt(1/2),
                    sqrt(1/2), 1, 0, 1/4, 0, 1/4,
                    0, 0, 1, sqrt(1/2), 1/2, 1/2*sqrt(1/2),
                    1/2*sqrt(1/2), 1/4, sqrt(1/2), 1, 1/2*sqrt(1/2), 1/2,
                    0, 0, 1/2, 1/2*sqrt(1/2), 1, sqrt(1/2),
                    1/2*sqrt(1/2), 1/4, 1/2*sqrt(1/2), 1/2, sqrt(1/2), 1),nrow = 6,byrow = T)
  
  
  #equation 1 
  
  Sigmapn1=Sigmapn4[c(1,2,3,5),c(1,2,3,5)]
  ubn1=c(u[1],Inf,l[1],l[1])
  lbn1=c(l[1],u[2],-Inf,-Inf)
  meansn1=meansn4[c(1,2,3,5)]
  set.seed(1)
  equation1=pmvnorm(lower=lbn1,upper=ubn1,mean = meansn1,corr = Sigmapn1)[[1]]
  
  
  #equation 2 
  
  Sigmapn2=Sigmapn4[c(1,2,3,5,6),c(1,2,3,5,6)]
  ubn2=c(u[1],Inf,l[1],u[1],u[2])
  lbn2=c(l[1],u[2],-Inf,l[1],-Inf)
  meansn2=meansn4[c(1,2,3,5,6)]
  set.seed(1)
  equation2=pmvnorm(lower=lbn2,upper=ubn2,mean = meansn2,corr = Sigmapn2)[[1]]
  
  #equation 3 
  
  Sigmapn3=Sigmapn4[c(1,2,3,4,5),c(1,2,3,4,5)]
  ubn3=c(u[1],Inf,u[1],u[2],l[1])
  lbn3=c(l[1],u[2],l[1],-Inf,-Inf)
  meansn3=meansn4[c(1,2,3,4,5)]
  set.seed(1)
  equation3=pmvnorm(lower=lbn3,upper=ubn3,mean = meansn3,corr = Sigmapn3)[[1]]
  
  #equation 4 
  
  ubn4=c(u[1],Inf,u[1],u[2],u[1],u[2])
  lbn4=c(l[1],u[2],l[1],-Inf,l[1],-Inf)
  set.seed(1)
  equation4=pmvnorm(lower=lbn4,upper=ubn4,mean = meansn4,corr = Sigmapn4)[[1]]
  
  result=equation1+equation2+equation3+equation4
  return(result)
}


OverallpowerOandN=function(theta0,sd,u,l,n) #Does it for whichever treatment is the best
{
  
  pkstar=which(theta0==max(theta0)) #possible kstars 
  OverallpowerNi=rep(NA,3)
  OverallpowerOi=rep(NA,3)
  if(any(pkstar==1))
  {
    kstar=1
    Part1=FirstrightT12(kstar,theta0,sd,u,l,n)
    Part2=SecondrightT12(kstar,theta0,sd,u,l,n)
    
    OPart3=WrongfirstODT12(kprime = 2,kstar,theta0,sd,u,l,n)
    NPart3=WrongfirstNDT12(kprime = 2,kstar,theta0,sd,u,l,n)
    
    OverallpowerNi[kstar]=Part1+Part2+NPart3
    
    OverallpowerOi[kstar]=Part1+Part2+OPart3
    

  }
  
  if(any(pkstar==2))
  {
    kstar=2
    Part1=FirstrightT12(kstar,theta0,sd,u,l,n)
    Part2=SecondrightT12(kstar,theta0,sd,u,l,n)
    
    OPart3=WrongfirstODT12(kprime = 1,kstar,theta0,sd,u,l,n)
    NPart3=WrongfirstNDT12(kprime = 1,kstar,theta0,sd,u,l,n)
    
    OverallpowerNi[kstar]=Part1+Part2+NPart3
    
    OverallpowerOi[kstar]=Part1+Part2+OPart3
    
  }
  
  if(any(pkstar==3))
  {
    kstar=3
    Part1=FirstrightT3(kstar,theta0,sd,u,l,n)
    #print("Part1")
    #print(Part1)
    
    Part2=SecondrightT3(kstar,theta0,sd,u,l,n)
    
    #print("Part2")
    #print(Part2)
    
    Part3a1=WrongfirstOaNDT3S2(kprime=1,kstar,theta0,sd,u,l,n)
    
    #print("Part3a1")
    #print(Part3a1)
    
    Part3b1=WrongfirstOaNDT3S2(kprime=2,kstar,theta0,sd,u,l,n)
    
    #print("Part3b1")
    #print(Part3b1)
    
    OPart3a2=WrongfirstODT3S2(kprime=1,kstar,theta0,sd,u,l,n)
    
    #print("OPart3a2")
    #print(OPart3a2)
    
    OPart3b2=WrongfirstODT3S2(kprime=2,kstar,theta0,sd,u,l,n)
    #print("OPart3b2")
    #print(OPart3b2)
    
    NPart3a2=WrongfirstNDT3S2(kprime=1,kstar,theta0,sd,u,l,n)
    NPart3b2=WrongfirstNDT3S2(kprime=2,kstar,theta0,sd,u,l,n)
    
    
    OverallpowerNi[kstar]=Part1+Part2+Part3a1+Part3b1+NPart3a2+NPart3b2
    
    OverallpowerOi[kstar]=Part1+Part2+Part3a1+Part3b1+OPart3a2+OPart3b2
    
    
  }  
  print(OverallpowerNi)
  print(OverallpowerOi)
  OverallpowerN=max(OverallpowerNi,na.rm = T)
  OverallpowerO=max(OverallpowerOi,na.rm = T)
  return(list("OverallpowerN"=OverallpowerN,"OverallpowerO"=OverallpowerO))
}


#thinking do 40 cores



library(doParallel)
library(foreach) 
library(parallel)
# Useful commands
load("triboundslater.Rdata")

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


AlldataOverpowerlateFocusT2agaistT1 <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    
    library(mvtnorm)
    OverallpowerOandN=OverallpowerOandN
    
    FirstrightT12=FirstrightT12
    SecondrightT12=SecondrightT12
    WrongfirstODT12=WrongfirstODT12
    WrongfirstNDT12=WrongfirstNDT12
    
    FirstrightT3=FirstrightT3
    SecondrightT3=SecondrightT3
    
    WrongfirstOaNDT3S2=WrongfirstOaNDT3S2
    WrongfirstODT3S2=WrongfirstODT3S2
    
    WrongfirstNDT3S2=WrongfirstNDT3S2
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=triboundslater[[3]][1,]
    tu=triboundslater[[2]][1,]
    tsd=1
    ans=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)
    Olddatagtheta0=ans$OverallpowerO
    Newdatagtheta0=ans$OverallpowerN
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataOverpowerlateFocusT2agaistT1,file = "AlldataOverpowerlateFocusT2agaistT1.Rdata")

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


AlldataOverpowerlateFocusT3agaistT1 <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    
    library(mvtnorm)
    OverallpowerOandN=OverallpowerOandN
    
    FirstrightT12=FirstrightT12
    SecondrightT12=SecondrightT12
    WrongfirstODT12=WrongfirstODT12
    WrongfirstNDT12=WrongfirstNDT12
    
    FirstrightT3=FirstrightT3
    SecondrightT3=SecondrightT3
    
    WrongfirstOaNDT3S2=WrongfirstOaNDT3S2
    WrongfirstODT3S2=WrongfirstODT3S2
    
    WrongfirstNDT3S2=WrongfirstNDT3S2
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=triboundslater[[3]][1,]
    tu=triboundslater[[2]][1,]
    tsd=1
    ans=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)
    Olddatagtheta0=ans$OverallpowerO
    Newdatagtheta0=ans$OverallpowerN
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataOverpowerlateFocusT3agaistT1,file = "AlldataOverpowerlateFocusT3agaistT1.Rdata")


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


AlldataOverpowerlateInfFocusT2agaistT1 <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    
    library(mvtnorm)
    OverallpowerOandN=OverallpowerOandN
    
    FirstrightT12=FirstrightT12
    SecondrightT12=SecondrightT12
    WrongfirstODT12=WrongfirstODT12
    WrongfirstNDT12=WrongfirstNDT12
    
    FirstrightT3=FirstrightT3
    SecondrightT3=SecondrightT3
    
    WrongfirstOaNDT3S2=WrongfirstOaNDT3S2
    WrongfirstODT3S2=WrongfirstODT3S2
    
    WrongfirstNDT3S2=WrongfirstNDT3S2
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=infboundslater[[3]][1,]
    tu=infboundslater[[2]][1,]
    tsd=1
    ans=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)
    Olddatagtheta0=ans$OverallpowerO
    Newdatagtheta0=ans$OverallpowerN
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataOverpowerlateInfFocusT2agaistT1,file = "AlldataOverpowerlateInfFocusT2agaistT1.Rdata")


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


AlldataOverpowerlateInfFocusT3agaistT1 <- foreach(x=1:dim(Mtheta0)[1]) %dopar% #for(x in 1:dim(Mtheta0)[1]) #
  {
    
    library(mvtnorm)
    OverallpowerOandN=OverallpowerOandN
    
    FirstrightT12=FirstrightT12
    SecondrightT12=SecondrightT12
    WrongfirstODT12=WrongfirstODT12
    WrongfirstNDT12=WrongfirstNDT12
    
    FirstrightT3=FirstrightT3
    SecondrightT3=SecondrightT3
    
    WrongfirstOaNDT3S2=WrongfirstOaNDT3S2
    WrongfirstODT3S2=WrongfirstODT3S2
    
    WrongfirstNDT3S2=WrongfirstNDT3S2
    
    ttheta0=Mtheta0[x,]
    tn=42
    tl=infboundslater[[3]][1,]
    tu=infboundslater[[2]][1,]
    tsd=1
    ans=OverallpowerOandN(ttheta0,tsd,tu,tl,tn)
    Olddatagtheta0=ans$OverallpowerO
    Newdatagtheta0=ans$OverallpowerN
    return(list("giventheta0"=ttheta0,"Olddatagtheta0"=Olddatagtheta0,"Newdatagtheta0"=Newdatagtheta0))
  }


stopCluster(cl) 

save(AlldataOverpowerlateInfFocusT3agaistT1,file = "AlldataOverpowerlateInfFocusT3agaistT1.Rdata")
