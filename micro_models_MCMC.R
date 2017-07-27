#This script is for estimating parameters for the models with 1 DOC pool using observations from East Long Lake.
#for questions and comments please contact Oleksandra "Sasha" Hararuk at ohararuk@gmail.com
#The logic of the code below as follows:
#     1. Load East Long Lake data and define additional input (entrainment flux, difference between stream and epilimnetic water densities)
#     2. Define initial parameter values and their prior ranges
#     3. Define elements of model structure
#     4. Calculate initial modeled epilimnetic CO2 and DOC pools
#     5. Estimate optimal model parameters
#     6. Calculate CO2 and DOC with best parameters and uncertainties around predictions
#     7. Plot CO2 and DOC pools and posterior parameter distributions

rm(list=ls())
library(MASS)
library(data.table)
library(R.utils)

#!!Important: set working directory here; code: setwd("C:/your_folder")
load('EnKF_LongData_20170223.RData') 
#variable description and units are located at https://github.com/jzwart/LakeCarbonEnKF/blob/master/Data/dataColumnDescriptions 

data<-data[1:134,] # discard data for 2015 (used later for validation)
obs<-cbind(data$dic,data$doc) #observed CO2 and DIC


par<-c(0.0048,40000,70000.0,7000000,6000)
# rdeath, Ea, Km, Vmax0, MIC_0
Min<-c(0.001,30000,30000,1000000,10)
Max<-c(0.01,48000,100000,100000000,7000)

run_type<-1 # 1=calibration, 2=validation
#estimate first guess pool values
source("MM_functions.R")
#for reverse MM kinetics: Xt1<-reverseMM(data,par,run_type)
#for MM kinetics: Xt1<-MM(data,par,run_type)
Xt1<-reverseMM(data,par,run_type)



#calculate initial value for likelihood function
n1<-as.numeric(length(na.omit(obs[,1])))
n2<-as.numeric(length(na.omit(obs[,2])))
error<-(Xt1-obs)^2
error[mapply(is.na,error)]<-0
error<-(colSums(error))
loglike<--0.5*(n1*log(error[1])+n2*log(error[2]))
J_old<- loglike

par_old<-par
par_keep<-par
J_keep<-J_old
simu<-1
par_rec<-par_old
n_par<-5 #calculating number of parameters for sampling and AIC calculation
sd<-2.38/sqrt(n_par) # from Gelman et al, 1996 Efficient metropolis jumping rules, in: Bernardo, J.M., Berger, J.O., Dawid, A.P., Smith, A.F.M. (eds.), Bayesian Statistics. Oxford University Press, pp. 599-607

#begin parameter estimation; details about the parameter proposal and acceptance criteria are in Haario et al., Bernoulli, 2001; Xu et al.GBC, 2006, 
for (simu in 1:500000) { # to get a representative parameter sample the routine can be stopped after 40,000 simulations
  
  while(TRUE){
    par_new<-par_old+runif(5,-0.5,0.5)*(Max-Min)/5 #propose new parameter set
    if(simu>10001){
      incr<-matrix(mvrnorm(n=1,matrix(0,1,5),as.matrix(Covars)*matrix(1,5,5),tol=1e-6),1,5)
      par_new<-par_old+incr
    }
    check1<-par_new < Max
    check2<-par_new > Min
    if(isTRUE(all(check1)==TRUE) && isTRUE(all(check2)==TRUE)){ #ensure that newly proposed parameters are within the prescribed boundaries
      break;
    }
  }
  par<-as.vector(par_new)
  
  #for reverse MM kinetics: Xt<-reverseMM(data,par,run_type)
  #for MM kinetics: Xt<-MM(data,par,run_type)
  Xt<-reverseMM(data,par,run_type)
  
  error<-(Xt-obs)^2
  error[mapply(is.na,error)]<-0
  error<-(colSums(error))
  
  loglike<--0.5*(n1*log(error[1])+n2*log(error[2]))
  J_new<- loglike
  
  delta<-J_new-J_old  #calculate the likelihood of the new parameters being accepted 
  
  if((min(1,exp(delta)))>runif(1,0,1)){
    
    par_old<-par
    J_old<-J_new
    J_keep<-cbind(J_keep,J_old)
    par_keep<-cbind(par_keep,par_old)
    
  }
  
  if(simu>10000){
    Covars<-cov(t(par_rec))
    Covars<-sd*Covars+diag(diag(sd*Covars),5,5)
  }
  print(simu)
  par_rec<-cbind(par_rec,par_old)
}

#calculating epilimnetic CO2 and DOC with best parameters and uncertainties around predictions

l<-length(par_keep[1,])
#discarding first half of the accepted parameters (burn-in period)
parf<-par_keep[,(l/2):l]
par<-rowMeans(parf)

load('EnKF_LongData_20170223.RData')
data<-data[135:256,] #2015 data
run_type<-2 # change the run type to "validation"

#for reverse MM kinetics: Xtf<-reverseMM(data,par,run_type)
#for MM kinetics: Xtf<-MM(data,par,run_type)
Xtf<-reverseMM(data,par,run_type)

obs<-cbind(data$dic,data$doc)
n1<-as.numeric(length(na.omit(obs[,1])))
n2<-as.numeric(length(na.omit(obs[,2])))
error<-(Xtf-obs)^2
error[mapply(is.na,error)]<-0
error<-(colSums(error))
loglike<--0.5*(n1*log(error[1])+n2*log(error[2]))

AIC<-2*n_par-2*loglike

#generate distribution of model predictions for each time point
chain<-500
Xt_multi<-array(0,dim=c(length(Xtf[,1]),length(Xtf[1,]),chain))


for (c in 1:chain){
  
  par<-parf[,c]
  #for reverse MM kinetics: Xt1<-reverseMM(data,par,run_type)
  #for MM kinetics: Xt1<-MM(data,par,run_type)
  Xt<-reverseMM(data,par,run_type)
  
  Xt_multi[,,c]<-Xt
  
}

library(matrixStats)
Xmean<-matrix(NA,length(Xtf[,1]),2)
Xsd<-matrix(NA,length(Xtf[,1]),2)
for (t in 1:length(Xtf[,1])){
  Xmean[t,]<-rowMeans(Xt_multi[t,,])
  Xsd[t,]<-rowSds(Xt_multi[t,,])
}


#Plot CO2 and DOC pools and posterior parameter distributions
par(mfrow=c(1,2))
plot(as.Date(data$datetime),Xmean[,1],ylab="CO2 (mol C)",type="l",col=rgb(0,0.5,0),lwd=2,xlab="date",ylim=c(0,4500))
polygon(c(as.Date(data$datetime),rev(as.Date(data$datetime))),c(Xmean[,1]+2*Xsd[,1],rev(Xmean[,1]-2*Xsd[,1])),col=rgb(0,1,0,0.5),border=FALSE)
lines(as.Date(data$datetime),data$dic,col="red",type="p",cex=1.2,pch=19)

plot(as.Date(data$datetime),Xmean[,2],ylab="DOC (mol C)",type="l",ylim=c(20000,90000),col=rgb(0,0.5,0),lwd=2,xlab="date")
polygon(c(as.Date(data$datetime),rev(as.Date(data$datetime))),c(Xmean[,2]+2*Xsd[,2],rev(Xmean[,2]-2*Xsd[,2])),col=rgb(0,1,0,0.5),border=FALSE)
lines(as.Date(data$datetime),data$doc,col="red",type="p",cex=1.2,pch=19)

par_names<-c("rdeath", "Ea", "Km", "Vmax0", "MIC_0")
units<-c("1/day","J/mol","mol","1/day","mol C")
plot.new()
par(mfrow=c(3,2))
for(i in 1:length(par)){
  hist(parf[i,],breaks=20,main=par_names[i],xlab=units[i],col="gray",xlim=c(Min[i],Max[i]))
}