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

#calculate entrainment flux
entr<-as.vector(matrix(0,length(data[,1]),1))
for (i in 2:length(data[,1])){
  entr[i-1]<-data$epiVol[i]-data$epiVol[i-1]-data$waterLoad[i-1]+data$waterLoss[i-1]
}
fromhypo<-entr
fromhypo<-replace(fromhypo,fromhypo<0,0)
tohypo<-entr
tohypo<-replace(tohypo,tohypo>0,0)

dens<-data$streamDens-data$epiDens

par<-c(0.007,1.047,0.007,0.99,0.99,0.0001,0)
# k20, q, kbase_photo, u1, u2, k_dens, k_photo
Min<-c(0.001,1,0.0001,0,0,0,-0.1)
Max<-c(0.02,1.2,0.02,1,1,1,1)

#define elements of model structure
run_type<-1 # 1=calibration, 2=validation
entr_switch<-1 # entrainment effect on C dynamics; 1=present, 0=absent
nu_hDIC_switch<-1 # assumption about depth gradient in hypolimnetic CO2; 1=present, 0=absent (can activate only if entrainment is included)
nu_hDOC_switch<-1 # assumption about depth gradient in hypolimnetic DOC; 1=present, 0=absent (can activate only if entrainment is included)
dens_switch<-1 # density effect on C load partitioning between epi- and hypolimnion; 1=present, 0=absent
photo_switch<-1 # photodegradation; 1=present, 0=absent

#estimate first guess pool values
source("1_pool_function.R")
Xt1<-OnePool(data,par,run_type,entr_switch,nu_hDIC_switch,nu_hDOC_switch,dens_switch,photo_switch)

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
n_par<-2+nu_hDIC_switch+nu_hDOC_switch+dens_switch+2*photo_switch #calculating number of parameters for sampling and AIC calculation
sd<-2.38/sqrt(n_par) # from Gelman et al, 1996 Efficient metropolis jumping rules, in: Bernardo, J.M., Berger, J.O., Dawid, A.P., Smith, A.F.M. (eds.), Bayesian Statistics. Oxford University Press, pp. 599-607

#begin parameter estimation; details about the parameter proposal and acceptance criteria are in Haario et al., Bernoulli, 2001; Xu et al.GBC, 2006, 
for (simu in 1:500000) { # to get a representative parameter sample the routine can be stopped after 40,000 simulations
  
  while(TRUE){
    par_new<-par_old+runif(7,-0.5,0.5)*(Max-Min)/5 #propose new parameter set
    if(simu>10001){
      incr<-matrix(mvrnorm(n=1,matrix(0,1,7),as.matrix(Covars)*matrix(1,7,7),tol=1e-6),1,7)
      par_new<-par_old+incr
    }
    check1<-par_new < Max
    check2<-par_new > Min
    if(isTRUE(all(check1)==TRUE) && isTRUE(all(check2)==TRUE)){ #ensure that newly proposed parameters are within the prescribed boundaries
      break;
    }
  }
  par<-as.vector(par_new)
  
  Xt<-OnePool(data,par,run_type,entr_switch,nu_hDIC_switch,nu_hDOC_switch,dens_switch,photo_switch)
  
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
    Covars<-sd*Covars+diag(diag(sd*Covars),7,7)
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
entr<-as.vector(matrix(0,length(data[,1]),1))
for (i in 2:length(data[,1])){
  entr[i-1]<-data$epiVol[i]-data$epiVol[i-1]-data$waterLoad[i-1]+data$waterLoss[i-1]
}
fromhypo<-entr
fromhypo<-replace(fromhypo,fromhypo<0,0)
tohypo<-entr
tohypo<-replace(tohypo,tohypo>0,0)
dens<-data$streamDens-data$epiDens

run_type<-2 # change the run type to "validation"

Xtf<-OnePool(data,par,run_type,entr_switch,nu_hDIC_switch,nu_hDOC_switch,dens_switch,photo_switch)

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
 Xt<-OnePool(data,par,run_type,entr_switch,nu_hDIC_switch,nu_hDOC_switch,dens_switch,photo_switch)
  
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

par_names<-c("k20","q","kbase_photo","u1","u2","k_dens","k_photo")
units<-c("1/day","unitless","1/day","unitless","unitless","unitless","unitless")
plot.new()
par(mfrow=c(4,2))
for(i in 1:length(par)){
  hist(parf[i,],breaks=20,main=par_names[i],xlab=units[i],col="gray",xlim=c(Min[i],Max[i]))
}