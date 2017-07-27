#The functions here are used in the script micro_models.R
MM<-function(data,par,run_type){
  
  Xt<-matrix(NA,length(data[,1]),3)
  
  Ea<-par[2]
  Vmax<-par[4]*exp(-Ea/(8.31*(data$epiTemp+273)))
  CUE<-0.5
  rdeath<-par[1]
  Km<-par[3]
 
  if(run_type==1){ #this is the first time point index observations for DOC and CO2 were available in 2014
    s<-6
  }
  if(run_type==2){ #this is the first time point index observations for DOC and CO2 were available in 2015
    s<-10
  }
  
  for (t in s:length(data[,1])){
    if (t==s){
      Xt[t-1,1]<-data$dic[t-1]
      Xt[t-1,2]<-data$doc[t-1]
      Xt[t-1,3]<-par[5]
      
    }
    
    
    change_dic<-0
    change_doc<-0
    change_mic<-0
    X1<-Xt[t-1,1]
    X2<-Xt[t-1,2]
    X3<-Xt[t-1,3]
    t2<-48
    for (t2 in 1:t2){
      change1<-(-(data$QoutInt[t-1]/data$epiVol[t-1])*X1-(data$kCO2[t-1]/data$thermo.depth[t-1])*X1+
                  ((1-CUE)*Vmax[t-1]*X2*X3/(Km+X2))+((data$DICeq[t-1]*data$epiVol[t-1]/data$thermo.depth[t])*data$kCO2[t-1])+data$dicIn[t-1]-
                  data$GPP[t-1]*(1-GPPrespired))/t2
      change2<-(-(Vmax[t-1]*X2*X3/(Km+X2))-(data$QoutInt[t-1]/data$epiVol[t-1])*X2+exudeTotal*data$GPP[t-1]+
                  data$docIn[t-1]+rdeath*X3)/t2
      change3<-((CUE*Vmax[t-1]*X2*X3/(Km+X2))-rdeath*X3)/t2
      X1<-X1+change1
      X2<-X2+change2
      X3<-X3+change3
    }
    
    Xt[t,1]<-X1
    Xt[t,2]<-X2
    Xt[t,3]<-X3
    
  }
  Xt<-Xt[,1:2]
  return(Xt)
}

reverseMM<-function(data,par,run_type){
  
  Xt<-matrix(NA,length(data[,1]),3)
  
  Ea<-par[2]
  Vmax<-par[4]*exp(-Ea/(8.31*(data$epiTemp+273)))
  CUE<-0.5
  rdeath<-par[1]
  Km<-par[3]
  
  if(run_type==1){ #this is the first time point index observations for DOC and CO2 were available in 2014
    s<-6
  }
  if(run_type==2){ #this is the first time point index observations for DOC and CO2 were available in 2015
    s<-10
  }
  
  for (t in s:length(data[,1])){
    if (t==s){
      Xt[t-1,1]<-data$dic[t-1]
      Xt[t-1,2]<-data$doc[t-1]
      Xt[t-1,3]<-par[5]
      
    }
    
    
    change_dic<-0
    change_doc<-0
    change_mic<-0
    X1<-Xt[t-1,1]
    X2<-Xt[t-1,2]
    X3<-Xt[t-1,3]
    t2<-48
    for (t2 in 1:t2){
      change1<-(-(data$QoutInt[t-1]/data$epiVol[t-1])*X1-(data$kCO2[t-1]/data$thermo.depth[t-1])*X1+
                  ((1-CUE)*Vmax[t-1]*X2*X3/(Km+X3))+((data$DICeq[t-1]*data$epiVol[t-1]/data$thermo.depth[t])*data$kCO2[t-1])+data$dicIn[t-1]-
                  data$GPP[t-1]*(1-GPPrespired))/t2
      change2<-(-(Vmax[t-1]*X2*X3/(Km+X3))-(data$QoutInt[t-1]/data$epiVol[t-1])*X2+exudeTotal*data$GPP[t-1]+
                  data$docIn[t-1]+rdeath*X3)/t2
      change3<-((CUE*Vmax[t-1]*X2*X3/(Km+X3))-rdeath*X3)/t2
      X1<-X1+change1
      X2<-X2+change2
      X3<-X3+change3
    }
    
    Xt[t,1]<-X1
    Xt[t,2]<-X2
    Xt[t,3]<-X3
    
  }
  Xt<-Xt[,1:2]
  return(Xt)
}